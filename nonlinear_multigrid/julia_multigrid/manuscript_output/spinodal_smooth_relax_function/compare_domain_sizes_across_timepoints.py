# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function"

indir_small = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output_large_domain"
phi_name = "MG_2000_dt_0.1_Nx_128_n_relax_4_eps_1.9211993488463508_large_domain_phi.txt"

phi_MG = np.genfromtxt(
    f"{indir_small}/{phi_name}",
)

phi_MG = phi_MG.reshape(-1, 128, 128).transpose(1, 2, 0)


# %%
indir_FD = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function"
phi_name_FD = "FD_20000_dt_6.10e-08_Nx_128_n_relax_4_gam_3.69e+00_D_16384_final_phi.csv"
phi_final_FD = np.genfromtxt(
    f"{indir_FD}/{phi_name_FD}",
    delimiter=",",
)

# %%
sses = []
for t in range(phi_MG.shape[2]):
    tmp = phi_MG[:, :, t]
    sse = np.sum((tmp - phi_final_FD) ** 2)
    sses.append(sse)

print("Matching timepoint =", np.argmin(sses), ", SSE =", sses[np.argmin(sses)])


# %%
def plot_differences(
    phi_MG,
    phi_final_FD,
    timepoint=0,
    title="",
    ax0_title="",
    ax1_title="",
    outdir="",
    save=False,
    suffix="",
):
    difference = phi_MG[:, :, timepoint] - phi_final_FD

    vcenter = 0
    vmin, vmax = difference.min(), difference.max()
    normalize_phis = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=-1, vmax=1)
    normalize_diff = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(1, 3, figsize=(24, 6))
    s = sns.heatmap(
        phi_MG[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        ax=axs[0],
        xticklabels=20,
        yticklabels=20,
        norm=normalize_phis,
    )
    axs[0].set_title(ax0_title)

    s = sns.heatmap(
        phi_final_FD,
        square=True,
        cmap=cm.RdBu_r,
        ax=axs[1],
        xticklabels=20,
        yticklabels=20,
        norm=normalize_phis,
    )
    axs[1].set_title(ax1_title)

    s = sns.heatmap(
        difference,
        square=True,
        norm=normalize_diff,
        cmap=cm.PiYG,
        ax=axs[2],
        xticklabels=20,
        yticklabels=20,
    )
    axs[2].set_title("Difference")
    plt.suptitle(title)
    plt.tight_layout()
    if save:
        plt.savefig(f"{outdir}/method_comparison{suffix}.png")
        plt.close()
    else:
        plt.show()


# %%
timepoint = np.argmin(sses)  # from MG
dt = 0.1
ns = 10
# note that rescaled t is the rescaled version due to the large domain size with M = 1 instead of 128^2.
title = "SSE =" + str(round(sses[timepoint], 3))
plot_differences(
    np.argmin(sses),
    title=title,
    ax0_title="NMG domain = [0,128], rescaled_t="
    + str(round((timepoint * dt * ns / 128**2), 6)),
    ax1_title="FD domain = [0,128], t=" + str(round(20000 * 6.10e-08, 6)),
)


# %%
indir_small = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"
phi_name_small = (
    "MG_2000_dt_6.103515625e-6_Nx_128_n_relax_4_eps_0.015009369912862116_phi.txt"
)
phi_M_small = np.genfromtxt(
    f"{indir_small}/{phi_name_small}",
)

phi_M_small = phi_M_small.reshape(-1, 128, 128).transpose(1, 2, 0)
# %%
s = sns.heatmap(
    phi_M_small[:, :, 2],
    square=True,
    cmap=cm.RdBu_r,
    xticklabels=20,
    yticklabels=20,
)

# %%
sses = []
for t in range(phi_M_small.shape[2]):
    tmp = phi_M_small[:, :, t]
    sse = np.sum((tmp - phi_final_FD) ** 2)
    sses.append(sse)

print("Matching timepoint =", np.argmin(sses), ", SSE =", sses[np.argmin(sses)])

# %%
timepoint = np.argmin(sses)  # from MG
dt = 6.103515625e-7
ns = 1000
# note that rescaled t is the rescaled version due to the large domain size with M = 1 instead of 128^2.
title = f"{phi_name_small} vs {phi_name_FD}; SSE =" + str(round(sses[timepoint], 3))
plot_differences(
    phi_M_small,
    phi_final_FD,
    np.argmin(sses),
    title=title,
    ax0_title="NMG domain = [0,1], t=" + str(round((timepoint * dt * ns), 6)),
    ax1_title="FD domain = [0,128], t=" + str(round(20000 * 6.10e-08, 6)),
    save=False,
    # outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/",
    # suffix="_large_FD_vs_small_MG",
)

# %%
