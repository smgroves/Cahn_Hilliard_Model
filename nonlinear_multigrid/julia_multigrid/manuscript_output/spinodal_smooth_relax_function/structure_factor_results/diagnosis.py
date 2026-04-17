"""
Diagnostic script: Help identify why structure factor analysis is failing.

Run this on your actual data to generate detailed diagnostics.
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2
from scipy.ndimage import gaussian_filter


def diagnose_phase_field(phi, title="Phase Field Diagnostics"):
    """
    Run basic sanity checks on a single phase field snapshot.

    Parameters
    ----------
    phi : ndarray (N, N)
        Phase field
    title : str
        Label for diagnostic output

    Returns
    -------
    dict with diagnostic metrics
    """

    print(f"\n{'='*70}")
    print(f"{title}")
    print(f"{'='*70}")

    # Basic statistics
    phi_mean = np.mean(phi)
    phi_std = np.std(phi)
    phi_min = np.min(phi)
    phi_max = np.max(phi)

    print(f"Shape: {phi.shape}")
    print(f"Mean:  {phi_mean:+.6f}")
    print(f"Std:   {phi_std:.6f}")
    print(f"Min:   {phi_min:+.6f}")
    print(f"Max:   {phi_max:+.6f}")
    print(f"Range: [{phi_min:+.4f}, {phi_max:+.4f}]")

    # Sanity checks
    checks = {}

    # Check 1: Phase field is finite
    if not np.all(np.isfinite(phi)):
        print("❌ ERROR: Phase field contains NaN or Inf!")
        checks['finite'] = False
    else:
        print("✓ Phase field is finite")
        checks['finite'] = True

    # Check 2: Phase field has variation
    if phi_std < 1e-6:
        print("❌ WARNING: Phase field is nearly constant (no variation)")
        checks['variation'] = False
    else:
        print(f"✓ Phase field has variation (std = {phi_std:.4f})")
        checks['variation'] = True

    # Check 3: Phase field is approximately normalized
    if abs(phi_mean) > 0.1:
        print(f"⚠ WARNING: Mean is {phi_mean:.4f} (expected ~0)")
    else:
        print(f"✓ Mean is close to zero ({phi_mean:+.6f})")

    if phi_std < 0.3 or phi_std > 2.0:
        print(f"⚠ WARNING: Std is {phi_std:.4f} (expected ~1)")
    else:
        print(f"✓ Std is reasonable ({phi_std:.4f})")

    # Check 4: Range is sensible
    if phi_max - phi_min < 0.5:
        print(f"❌ WARNING: Range is tiny ({phi_max - phi_min:.4f})")
        checks['range'] = False
    else:
        print(f"✓ Range is adequate ({phi_max - phi_min:.4f})")
        checks['range'] = True

    # Check 5: Structure is present
    grad_mag = np.sqrt(np.gradient(phi, axis=0)**2 +
                       np.gradient(phi, axis=1)**2)
    grad_mean = np.mean(grad_mag)

    if grad_mean < 1e-6:
        print(f"❌ ERROR: No spatial structure (gradients ~0)")
        checks['structure'] = False
    else:
        print(f"✓ Spatial structure present (avg |∇φ| = {grad_mean:.6f})")
        checks['structure'] = True

    return checks


def diagnose_time_series(phi_list, times=None, n_sample=3):
    """
    Run diagnostics across a time series.

    Parameters
    ----------
    phi_list : list of ndarray
        Phase fields at successive times
    times : array-like, optional
        Time values
    n_sample : int
        Number of snapshots to diagnose
    """

    n_total = len(phi_list)
    if times is None:
        times = np.arange(n_total)

    print(f"\n{'='*70}")
    print("TIME SERIES DIAGNOSTICS")
    print(f"{'='*70}")
    print(f"Total snapshots: {n_total}")
    print(f"Time range: {times[0]:.3f} to {times[-1]:.3f}")

    # Sample early, middle, late
    if n_total <= 3:
        indices = list(range(n_total))
    else:
        indices = [0, n_total // 2, n_total - 1]

    for i in indices:
        t = times[i]
        phi = phi_list[i]
        diagnose_phase_field(phi, f"Snapshot {i} (t = {t:.3f})")


def diagnose_coarsening(phi_list, times=None, L=None, N=None):
    """
    Check if domains are actually coarsening.

    Computes simple size metrics (no FFT) to see if l(t) is growing.

    Parameters
    ----------
    phi_list : list of ndarray
        Phase fields
    times : array-like, optional
    L, N : float, int
        Domain size and resolution (optional)
    """

    if times is None:
        times = np.arange(len(phi_list))

    print(f"\n{'='*70}")
    print("COARSENING DIAGNOSTICS (No FFT)")
    print(f"{'='*70}")

    # Simple metric: fraction of interface pixels
    interface_fracs = []

    for phi in phi_list:
        # Threshold gradient to find interfaces
        grad_mag = np.sqrt(np.gradient(phi, axis=0)**2 +
                           np.gradient(phi, axis=1)**2)
        threshold = np.mean(grad_mag) + np.std(grad_mag)
        interface_pixels = np.sum(grad_mag > threshold)
        total_pixels = phi.size
        frac = interface_pixels / total_pixels
        interface_fracs.append(frac)

    interface_fracs = np.array(interface_fracs)

    print(f"Interface pixel fraction over time:")
    for i, (t, frac) in enumerate(zip(times, interface_fracs)):
        print(f"  t={t:8.3f}: {frac:.4f}")

    # Check if decreasing
    early_frac = np.mean(interface_fracs[:max(1, len(interface_fracs)//4)])
    late_frac = np.mean(interface_fracs[int(0.75*len(interface_fracs)):])

    print(f"\nAverage interface fraction:")
    print(f"  Early times: {early_frac:.4f}")
    print(f"  Late times:  {late_frac:.4f}")

    if late_frac < early_frac * 0.8:
        print(f"✓ Interface fraction decreasing → domains coarsening")
        return True
    elif late_frac < early_frac * 1.2:
        print(f"⚠ Interface fraction roughly constant → slow/no coarsening")
        return False
    else:
        print(f"❌ Interface fraction increasing → unphysical (check dynamics?)")
        return False


def check_energy_conservation(phi_list, times=None, alpha=0.01):
    """
    Check if free energy is monotonically decreasing.

    For Cahn-Hilliard, should have dF/dt ≤ 0 everywhere.

    Parameters
    ----------
    phi_list : list of ndarray
    times : array-like, optional
    alpha : float
        Parameter in F = ∫(φ²/4 + α|∇φ|²/2)
    """

    if times is None:
        times = np.arange(len(phi_list))

    print(f"\n{'='*70}")
    print("ENERGY CONSERVATION CHECK")
    print(f"{'='*70}")

    energies = []

    for phi in phi_list:
        # F = ∫(φ⁴/4 - φ²/2 + α|∇φ|²/2) (roughly)
        # Simplified: F ~ ∫(φ² + α|∇φ|²)
        bulk_energy = np.sum(phi**2)
        grad_x = np.gradient(phi, axis=0)
        grad_y = np.gradient(phi, axis=1)
        grad_energy = alpha * np.sum(grad_x**2 + grad_y**2)
        F = bulk_energy + grad_energy
        energies.append(F)

    energies = np.array(energies)

    print(f"Free energy over time:")
    for t, F in zip(times, energies):
        print(f"  t={t:8.3f}: F={F:12.6e}")

    # Check monotonicity
    dF = np.diff(energies)
    increasing = np.sum(dF > 0)
    decreasing = np.sum(dF < 0)

    print(f"\nEnergy changes:")
    print(f"  Decreasing: {decreasing}/{len(dF)}")
    print(f"  Increasing: {increasing}/{len(dF)}")

    if increasing > len(dF) * 0.1:
        print(
            f"❌ WARNING: Energy increases at {increasing} steps (not physical)")
        return False
    else:
        print(f"✓ Energy is monotonically decreasing")
        return True


def check_mass_conservation(phi_list, times=None):
    """
    Check if total mass ∫φ dA is conserved.

    For Cahn-Hilliard, should be constant.
    """

    if times is None:
        times = np.arange(len(phi_list))

    print(f"\n{'='*70}")
    print("MASS CONSERVATION CHECK")
    print(f"{'='*70}")

    masses = [np.sum(phi) for phi in phi_list]
    masses = np.array(masses)

    mass_ref = masses[0]
    mass_changes = np.abs(masses - mass_ref) / (np.abs(mass_ref) + 1e-10)

    print(f"Total mass over time:")
    for t, m, dm in zip(times, masses, mass_changes):
        print(f"  t={t:8.3f}: M={m:+.6e}  (Δ={dm*100:6.3f}%)")

    max_drift = np.max(mass_changes)

    if max_drift > 0.01:
        print(
            f"❌ WARNING: Mass drifts by {max_drift*100:.2f}% (should be <1%)")
        return False
    else:
        print(f"✓ Mass is conserved (drift = {max_drift*100:.4f}%)")
        return True


def diagnose_structure_factor(phi, L, N, title="Structure Factor"):
    """
    Check if structure factor computation makes sense.

    Parameters
    ----------
    phi : ndarray (N, N)
    L : float
        Domain size
    N : int
        Grid points
    """

    print(f"\n{'='*70}")
    print(f"{title}")
    print(f"{'='*70}")

    # Compute S(k) manually
    phi_centered = phi - np.mean(phi)
    phi_hat = fft2(phi_centered)
    power = np.abs(phi_hat) ** 2
    power_shifted = np.fft.fftshift(power)

    # Find peak
    center = N // 2
    yy, xx = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    rr = np.sqrt((xx - center)**2 + (yy - center)**2)

    # Radial average
    radii = np.arange(1, N // 2)
    S_k = np.zeros(len(radii))
    for i, r in enumerate(radii):
        mask = (rr >= r - 0.5) & (rr < r + 0.5)
        if np.any(mask):
            S_k[i] = np.mean(power_shifted[mask])

    k_grid = 2 * np.pi * radii / L

    # Find peak
    skip = max(2, int(0.05 * len(radii)))
    idx_max = np.argmax(S_k[skip:]) + skip
    k_max = k_grid[idx_max]
    l_char = 2 * np.pi / k_max if k_max > 0 else np.inf

    print(f"Peak wavenumber: k_max = {k_max:.6f}")
    print(f"Characteristic length: l = 2π/k = {l_char:.6f}")

    # Check for reasonable peak
    peak_height = S_k[idx_max]
    avg_height = np.mean(S_k[skip:])

    if peak_height < 2 * avg_height:
        print(
            f"⚠ WARNING: Peak is not very pronounced (ratio = {peak_height/avg_height:.2f}x)")
        print(f"  → S(k) may be too noisy; try smoothing or higher resolution")
    else:
        print(
            f"✓ Peak is well-defined (ratio = {peak_height/avg_height:.2f}x)")

    # Check against domain size
    if l_char > L * 0.3:
        print(f"⚠ WARNING: Length scale is large relative to domain")
        print(f"  → May be in saturation regime; exclude late times from fit")
    elif l_char < L / 100:
        print(f"⚠ WARNING: Length scale is very small")
        print(f"  → May need finer grid; currently limited by resolution")
    else:
        print(f"✓ Length scale is reasonable (l/L = {l_char/L:.3f})")

    return {'k_max': k_max, 'l_char': l_char, 'S_k': S_k, 'k_grid': k_grid}

# %%
# ============================================================================
# MAIN: Run all diagnostics
# ============================================================================


if __name__ == '__main__':

    print("""
    ╔════════════════════════════════════════════════════════════════════════╗
    ║            CAHN-HILLIARD SIMULATION DIAGNOSTICS                       ║
    ║          Help identify why structure factor is failing                 ║
    ╚════════════════════════════════════════════════════════════════════════╝
    """)

    # EDIT THESE FOR YOUR DATA:
    # phi_list = [...]  # Your phase fields
    # times = np.array([...])
    # L = 64.0
    # N = 256

    print("\n" + "="*70)
    print("USAGE:")
    print("="*70)
    print("""
    1. Edit this file: set phi_list, times, L, N at the bottom
    2. Run: python diagnose.py
    3. Review the output for warnings and errors
    4. Fix identified issues
    5. Re-run structure factor analysis
    """)

    print("\nExample code to add:\n")
    print("""
    import numpy as np
    from structure_factor_simple import StructureFactorAnalyzer
    
    # Load your data
    phi_data = np.load('your_simulation.npy')
    phi_list = list(phi_data)
    times = np.logspace(0, 3, len(phi_list))  # Adjust time range
    L = 64.0
    N = 256
    
    # Run diagnostics
    diagnose_time_series(phi_list, times, n_sample=5)
    diagnose_coarsening(phi_list, times, L, N)
    check_energy_conservation(phi_list, times)
    check_mass_conservation(phi_list, times)
    
    if len(phi_list) > 0:
        # Check first and last snapshots
        diagnose_structure_factor(phi_list[0], L, N, "First snapshot")
        diagnose_structure_factor(phi_list[-1], L, N, "Last snapshot")
    """)
