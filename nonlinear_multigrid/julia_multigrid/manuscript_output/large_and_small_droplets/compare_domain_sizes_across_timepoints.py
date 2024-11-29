# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets"

indir_small = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/output"
phi_name = (
    "MG_5000_dt_0.001_Nx_128_eps_0.015009369912862116_r1_20_r2_30_space_10_phi.txt"
)

phi_MG = np.genfromtxt(
    f"{indir_small}/{phi_name}",
)

phi_MG = phi_MG.reshape(-1, 128, 128).transpose(1, 2, 0)


# %%
indir_FD = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/large_and_small_droplets"
phi_name_FD = "FD_3000000_dt_1.00e-06_Nx_128_gam_1.48e+01_D_16384_r1_20_r2_30_space_10_final_phi.csv"
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
def plot_differences(timepoint=0, title=""):
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
    axs[0].set_title("NMG")

    s = sns.heatmap(
        phi_final_FD,
        square=True,
        cmap=cm.RdBu_r,
        ax=axs[1],
        xticklabels=20,
        yticklabels=20,
        norm=normalize_phis,
    )
    axs[1].set_title("FD")

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
    plt.show()


# plt.savefig(f"{outdir}/method_comparison_NMG_vs_FD_0_128.png")


# %%
timepoint = np.argmin(sses)
title = (
    "Matching timepoint =" + str(timepoint) + ", SSE =" + str(round(sses[timepoint], 3))
)
plot_differences(np.argmin(sses), title=title)

timepoint = 90
title = (
    "Matching timepoint =" + str(timepoint) + ", SSE =" + str(round(sses[timepoint], 3))
)

plot_differences(90, title=title)

# %%
plt.plot(phi_MG[:, 64, np.argmin(sses)], label="MG")
plt.plot(phi_final_FD[:, 64], label="FD")
plt.legend()
plt.title("Vertical profile at y = 64")
plt.show()
plt.plot(phi_MG[64, :, np.argmin(sses)], label="MG")
plt.plot(phi_final_FD[64, :], label="FD")
plt.legend()
plt.title("Horizontal profile at x = 64")
plt.show()

# %%
