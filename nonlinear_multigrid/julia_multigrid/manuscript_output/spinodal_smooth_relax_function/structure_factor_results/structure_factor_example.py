"""
SPINODAL DECOMPOSITION: S(k) Evolution with Visual Phase Field Correspondence

Shows:
1. Phase field φ(x,y) at different times
2. Structure factor S(k) at those same times
3. How peak wavenumber k_max relates to visual domain size
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2
from matplotlib.gridspec import GridSpec


def generate_synthetic_coarsening(N=128, n_times=8, random_seed=42):
    """
    Generate synthetic coarsening data for visualization.

    Simulates spinodal decomposition:
    - Start: random noise (small domains everywhere)
    - End: large bicontinuous structure
    """

    np.random.seed(random_seed)

    times = np.logspace(0, 3, n_times)  # t = 1, 3.16, 10, 31.6, ..., 1000
    phi_list = []

    for t in times:
        # Synthetic coarsening: domain size grows as t^(1/3)
        l_t = 2 * (t ** (1/3))  # Characteristic length

        # Generate phase field with given characteristic length
        # Use combination of sine waves at appropriate scale
        x = np.linspace(0, 2*np.pi, N)
        y = np.linspace(0, 2*np.pi, N)
        xx, yy = np.meshgrid(x, y)

        # Wave vector corresponding to characteristic length
        k_dominant = 2*np.pi / l_t

        # Sum of sine waves at the dominant wavenumber + noise
        phi = (np.sin(k_dominant * xx) * np.cos(k_dominant * yy) +
               0.3 * np.random.randn(N, N))

        # Apply smooth threshold to make it look binary-ish
        phi = np.tanh(phi)

        phi_list.append(phi)

    return phi_list, times


def compute_structure_factor(phi, L=2*np.pi):
    """
    Compute structure factor S(k) from phase field.
    """

    # Remove mean
    phi_centered = phi - np.mean(phi)

    # FFT
    phi_hat = fft2(phi_centered)
    power = np.abs(phi_hat) ** 2

    # Shift to center
    power_shifted = np.fft.fftshift(power)

    # Radial averaging
    N = phi.shape[0]
    center = N // 2
    coords = np.arange(-center, center)
    yy, xx = np.meshgrid(coords, coords, indexing='ij')
    r_pix = np.sqrt(xx**2 + yy**2)

    n_bins = N // 2
    radii = np.arange(n_bins, dtype=float)
    S_k = np.zeros(n_bins)

    for i in range(1, n_bins):
        mask = (r_pix >= i - 0.5) & (r_pix < i + 0.5)
        if np.any(mask):
            S_k[i] = np.mean(power_shifted[mask])

    # Convert to wavenumber
    k = 2 * np.pi * radii / L

    # Find peak
    skip_idx = max(2, int(0.05 * n_bins))
    idx_max = np.argmax(S_k[skip_idx:]) + skip_idx
    k_max = k[idx_max]
    l_char = 2 * np.pi / k_max if k_max > 0 else np.inf

    return k, S_k, k_max, l_char


# Generate data
print("Generating synthetic spinodal decomposition data...")
phi_list, times = generate_synthetic_coarsening(N=128, n_times=8)

# Compute structure factors
S_k_list = []
k_values = None
l_values = []
k_max_values = []

for phi in phi_list:
    k, S_k, k_max, l_char = compute_structure_factor(phi)
    S_k_list.append(S_k)
    l_values.append(l_char)
    k_max_values.append(k_max)
    if k_values is None:
        k_values = k

# Create figure with phase fields and S(k) side by side
fig = plt.figure(figsize=(16, 10))
gs = GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.3)

# Select 4 times to show (early, 2 middle, late)
if len(times) >= 4:
    indices = [0, len(times)//3, 2*len(times)//3, len(times)-1]
else:
    indices = list(range(len(times)))

print(f"\nShowing snapshots at times: {times[indices]}")
print(f"Characteristic lengths: {[l_values[i] for i in indices]}")
print(f"Peak wavenumbers k_max: {[k_max_values[i] for i in indices]}")

for col, idx in enumerate(indices):
    t = times[idx]
    phi = phi_list[idx]
    S_k = S_k_list[idx]
    k_max = k_max_values[idx]
    l_char = l_values[idx]

    # ========== PHASE FIELD (top row) ==========
    ax_phi = fig.add_subplot(gs[0, col])

    im = ax_phi.imshow(phi, cmap='RdBu_r', origin='lower', vmin=-1, vmax=1)
    ax_phi.set_title(f't = {t:.2f}\nl = {l_char:.3f}',
                     fontsize=11, fontweight='bold')
    ax_phi.set_xlabel('x')
    ax_phi.set_ylabel('y')
    plt.colorbar(im, ax=ax_phi, label='φ', fraction=0.046, pad=0.04)

    # ========== STRUCTURE FACTOR S(k) (middle row, log scale) ==========
    ax_sk = fig.add_subplot(gs[1, col])

    # Plot S(k) with peak highlighted
    ax_sk.semilogy(k_values[1:], S_k[1:], 'o-',
                   linewidth=1.5, markersize=4, color='steelblue')
    ax_sk.axvline(k_max, color='red', linestyle='--',
                  linewidth=2, label=f'k_max = {k_max:.3f}')
    ax_sk.set_xlabel('Wavenumber k', fontsize=10)
    ax_sk.set_ylabel('S(k)', fontsize=10)
    ax_sk.set_title(f'Structure Factor', fontsize=11, fontweight='bold')
    ax_sk.grid(True, alpha=0.3, which='both')
    ax_sk.legend(fontsize=9)
    ax_sk.set_xlim(left=0)

    # ========== INTERPRETATION (bottom row) ==========
    ax_text = fig.add_subplot(gs[2, col])
    ax_text.axis('off')

    # Interpretation text
    if idx == indices[0]:
        status = "EARLY TIME\nSmall domains"
        meaning = f"High k_max\n→ Many small features"
    elif idx == indices[-1]:
        status = "LATE TIME\nLarge domains"
        meaning = f"Low k_max\n→ Few large features"
    else:
        status = "INTERMEDIATE"
        meaning = f"k_max decreasing\n→ Coarsening"

    text = f"""{status}

k_max = {k_max:.4f}
l = 2π/k_max = {l_char:.3f}

{meaning}

Domain size ∝ l(t)
Peak moves LEFT over time!"""

    ax_text.text(0.5, 0.5, text,
                 ha='center', va='center',
                 fontsize=10,
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                 family='monospace')

# Overall title
fig.suptitle('Spinodal Decomposition: Phase Field → Structure Factor S(k) Evolution\n' +
             'Peak k_max moves LEFT (decreasing wavenumber) as domains GROW',
             fontsize=14, fontweight='bold', y=0.995)

plt.savefig('spinodal_sk_evolution.png', dpi=150, bbox_inches='tight')
print("\n✓ Saved: spinodal_sk_evolution.png")
plt.show()


# ============================================================================
# ADDITIONAL FIGURE: k_max vs time (coarsening law)
# ============================================================================

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: k_max vs t
ax1.semilogy(times, k_max_values, 'o-', linewidth=2.5,
             markersize=8, color='blue')
ax1.set_xlabel('Time t', fontsize=12, fontweight='bold')
ax1.set_ylabel('Peak wavenumber k_max', fontsize=12, fontweight='bold')
ax1.set_title('Peak Wavenumber Decreases (Coarsening)',
              fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3, which='both')
ax1.invert_yaxis()  # Invert so "large k" is at top (intuitive)

# Annotate
ax1.text(times[0], k_max_values[0]*1.3, 'Small domains\n(high k)', fontsize=10, ha='center',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
ax1.text(times[-1], k_max_values[-1]*0.7, 'Large domains\n(low k)', fontsize=10, ha='center',
         bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))

# Right: l(t) vs t
ax2.loglog(times, l_values, 's-', linewidth=2.5, markersize=8, color='green')
ax2.set_xlabel('Time t', fontsize=12, fontweight='bold')
ax2.set_ylabel('Characteristic length l(t)', fontsize=12, fontweight='bold')
ax2.set_title('Domain Size Grows (l ∝ t^(1/3))',
              fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3, which='both')

# Theory line
l_theory = l_values[0] * (times / times[0]) ** (1/3)
ax2.loglog(times, l_theory, '--', linewidth=2, color='gray',
           alpha=0.6, label='Theory: l ∝ t^(1/3)')
ax2.legend(fontsize=11)

fig2.suptitle('Coarsening Kinetics: How Peak Position Reflects Domain Growth',
              fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('coarsening_kinetics.png', dpi=150, bbox_inches='tight')
print("✓ Saved: coarsening_kinetics.png")
plt.show()


# ============================================================================
# REFERENCE: What different k_max values mean visually
# ============================================================================

print("\n" + "="*70)
print("REFERENCE: WHAT DIFFERENT k_max VALUES MEAN VISUALLY")
print("="*70)
print("""
k_max (Peak Wavenumber)  |  l = 2π/k_max  |  What You See Visually
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

HIGH k_max (e.g., 5.0)   |  SMALL l~1.3   |  Lots of small domains
                         |                |  Fine-scale oscillations
                         |                |  EARLY TIME
                         |
MEDIUM k_max (e.g., 1.0) |  MEDIUM l~6.3  |  Medium-sized domains
                         |                |  Mix of features
                         |                |  INTERMEDIATE TIME
                         |
LOW k_max (e.g., 0.3)    |  LARGE l~21    |  Few large domains
                         |                |  Bicontinuous structure
                         |                |  LATE TIME

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

KEY INSIGHT:
    
    As time goes on: k_max DECREASES
                     ↓
                     l increases
                     ↓
                     Domains grow larger
                     ↓
                     COARSENING!

The peak of S(k) tells you the DOMINANT DOMAIN SIZE at that moment.
As domains merge, the peak moves LEFT (to lower k).
""")

print("\nPeak wavenumber progression:")
for i, idx in enumerate(indices):
    print(
        f"  Step {i+1}: t={times[idx]:7.2f}  →  k_max={k_max_values[idx]:.4f}  →  l={l_values[idx]:.4f}")

print("\n" + "="*70)
