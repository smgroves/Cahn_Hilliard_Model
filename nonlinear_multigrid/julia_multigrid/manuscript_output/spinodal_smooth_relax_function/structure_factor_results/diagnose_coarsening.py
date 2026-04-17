"""
INVESTIGATION: Why does the interface metric not detect coarsening 
               even though you see it visually?

This is actually GOOD NEWS - it usually means the coarsening is real,
but the metric is wrong (not that your simulation is broken).
"""
# %%
import numpy as np
import matplotlib.pyplot as plt


def investigate_interface_metric(phi_list, times, verbose=True):
    """
    Deep dive into interface detection. Figure out why the metric might fail.
    """

    print("\n" + "="*70)
    print("INVESTIGATING INTERFACE METRIC MISMATCH")
    print("="*70)

    # Collect metrics with MULTIPLE thresholds
    thresholds_to_try = [
        ('Very loose (mean only)', lambda grad: np.mean(grad)),
        ('Loose (mean + 0.5*std)', lambda grad: np.mean(grad) + 0.5*np.std(grad)),
        ('Standard (mean + std)', lambda grad: np.mean(grad) + np.std(grad)),
        ('Tight (mean + 2*std)', lambda grad: np.mean(grad) + 2*np.std(grad)),
        ('Very tight (95th percentile)', lambda grad: np.percentile(grad, 95)),
    ]

    all_fracs = {}

    for label, thresh_fn in thresholds_to_try:
        fracs = []
        for phi in phi_list:
            grad_x = np.gradient(phi, axis=0)
            grad_y = np.gradient(phi, axis=1)
            grad_mag = np.sqrt(grad_x**2 + grad_y**2)
            thresh = thresh_fn(grad_mag)
            frac = np.sum(grad_mag > thresh) / phi.size
            fracs.append(frac)

        all_fracs[label] = np.array(fracs)

        # Check if this threshold shows coarsening
        early = np.mean(fracs[:max(1, len(fracs)//4)])
        late = np.mean(fracs[int(0.75*len(fracs)):])

        trend = "declining" if late < 0.95 * \
            early else ("flat" if late < 1.05*early else "increasing")

        if verbose:
            print(f"\n{label}:")
            print(f"  Early: {early:.6f}, Late: {late:.6f}")
            print(f"  Trend: {trend}")

    # Which threshold shows coarsening?
    print("\n" + "-"*70)
    print("SUMMARY: Which threshold detects coarsening?")
    print("-"*70)

    coarsening_detected = False
    for label, fracs in all_fracs.items():
        early = np.mean(fracs[:max(1, len(fracs)//4)])
        late = np.mean(fracs[int(0.75*len(fracs)):])
        if late < 0.95*early:
            print(f"✓ {label}")
            coarsening_detected = True
        else:
            print(f"✗ {label}")

    if not coarsening_detected:
        print("\n⚠ NO THRESHOLD DETECTS COARSENING")
        print("This means:")
        print("  → Domain morphology is changing in confusing ways")
        print("  → Metric is fundamentally wrong for your system")
        print("  → Check spectrum width instead (Method 2)")

    # Plot all thresholds
    fig, ax = plt.subplots(figsize=(12, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(all_fracs)))

    for (label, fracs), color in zip(all_fracs.items(), colors):
        ax.plot(times, fracs, 'o-', label=label,
                color=color, linewidth=2, markersize=5)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Interface pixel fraction', fontsize=12)
    ax.set_title('Interface Metric with Different Thresholds',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    fig.savefig('interface_threshold_comparison.png',
                dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved: interface_threshold_comparison.png")

    return all_fracs


def check_domain_morphology(phi_list, times):
    """
    Check if domain SHAPE is changing in ways that confuse the interface metric.

    Example: If domains go from "large connected blobs" to "small disconnected dots",
    interface fraction might INCREASE even though domains are coarsening.
    """

    print("\n" + "="*70)
    print("CHECKING DOMAIN MORPHOLOGY (Shape Changes)")
    print("="*70)

    # Sample early, middle, late phases
    indices = [0, len(phi_list)//2, len(phi_list)-1]
    labels = ['Early', 'Middle', 'Late']

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    for ax_idx, (i, label) in enumerate(zip(indices, labels)):
        phi = phi_list[i]

        # Threshold at +/- 0.5 (roughly the phase boundary)
        boundary_pixels = np.abs(phi) < 0.1  # Very thin interface
        domain_plus = phi > 0.5
        domain_minus = phi < -0.5

        # Create visualization
        viz = np.zeros((phi.shape[0], phi.shape[1], 3))
        viz[domain_plus] = [1, 0.7, 0]   # Orange: + domain
        viz[domain_minus] = [0, 0.7, 1]  # Blue: - domain
        viz[boundary_pixels] = [1, 1, 1]  # White: interface

        axes[ax_idx].imshow(viz, origin='lower')
        axes[ax_idx].set_title(
            f'{label} (t={times[i]:.3f})', fontsize=11, fontweight='bold')
        axes[ax_idx].set_xlabel('x')
        axes[ax_idx].set_ylabel('y')

    plt.tight_layout()
    fig.savefig('domain_morphology_evolution.png',
                dpi=150, bbox_inches='tight')
    print(f"✓ Saved: domain_morphology_evolution.png")

    print("""
    VISUAL CHECK:
    - Orange regions = + phase
    - Blue regions = - phase
    - White pixels = interface
    
    If you see:
    ✓ Fewer, larger orange/blue blobs → Real coarsening
    ✗ Many small fragments → Fragmentation (not coarsening)
    ~ Changing from blobs to foam-like → Topology change
    """)


def check_gradient_histogram(phi_list, times):
    """
    Check the distribution of gradients over time.

    If |∇φ| distribution changes shape (not just shifts),
    a fixed threshold will give confusing results.
    """

    print("\n" + "="*70)
    print("CHECKING GRADIENT DISTRIBUTIONS")
    print("="*70)

    # Sample at early, middle, late times
    indices = [0, len(phi_list)//2, len(phi_list)-1]
    labels = ['Early', 'Middle', 'Late']

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    for ax_idx, (i, label) in enumerate(zip(indices, labels)):
        phi = phi_list[i]
        grad_x = np.gradient(phi, axis=0)
        grad_y = np.gradient(phi, axis=1)
        grad_mag = np.sqrt(grad_x**2 + grad_y**2)

        ax = axes[ax_idx]
        ax.hist(grad_mag.flatten(), bins=50, color='steelblue',
                alpha=0.7, edgecolor='black')
        ax.set_xlabel('|∇φ|', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(f'{label} (t={times[i]:.3f})',
                     fontsize=11, fontweight='bold')
        ax.set_yscale('log')

        # Mark some reference thresholds
        mean_val = np.mean(grad_mag)
        std_val = np.std(grad_mag)
        ax.axvline(mean_val, color='green', linestyle='--',
                   linewidth=2, label=f'Mean={mean_val:.3f}')
        ax.axvline(mean_val + std_val, color='red', linestyle='--',
                   linewidth=2, label=f'Mean+std={mean_val+std_val:.3f}')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig('gradient_distribution_evolution.png',
                dpi=150, bbox_inches='tight')
    print(f"✓ Saved: gradient_distribution_evolution.png")

    print("""
    INTERPRETATION:
    - If histograms shift LEFT over time → Interfaces are SHARPENING
      (coarsening means fewer, sharper interfaces)
    - If histograms stay SAME → Coarsening not detected
    - If histograms shift RIGHT → Interfaces are BLURRING
      (could indicate: saturation, numerical dissipation, or saturation)
    """)


# ============================================================================
# MAIN
# ============================================================================

if __name__ == '__main__':

    print("""
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                  VISUAL vs METRIC INVESTIGATION                     ║
    ║                                                                      ║
    ║  If you see coarsening visually but metrics say "no coarsening"    ║
    ║  this script will tell you WHY and which method to trust.          ║
    ╚══════════════════════════════════════════════════════════════════════╝
    """)

    print("""
    USAGE:
    
    import numpy as np
    from visual_vs_metric_diagnostic import (
        investigate_interface_metric,
        check_domain_morphology,
        check_gradient_histogram
    )
    
    # Your data
    phi_list = [...]
    times = np.array([...])
    
    # Run investigations
    all_fracs = investigate_interface_metric(phi_list, times)
    check_domain_morphology(phi_list, times)
    check_gradient_histogram(phi_list, times)
    
    # This creates 3 diagnostic plots:
    # - interface_threshold_comparison.png
    # - domain_morphology_evolution.png
    # - gradient_distribution_evolution.png
    """)
# %%
