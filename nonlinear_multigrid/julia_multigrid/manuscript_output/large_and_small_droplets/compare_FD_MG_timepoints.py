# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function"

indir_small = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"
phi_name = "MG_200_dt_5.5e-7_Nx_128_n_relax_4_eps_0.015009369912862116_phi.txt"

phi_MG = np.genfromtxt(
    f"{indir_small}/{phi_name}",
)

phi_MG = phi_MG.reshape(-1, 128, 128).transpose(1, 2, 0)

# %%
indir_FD = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function"
phi_name_FD = "FD_200_dt_5.50e-07_Nx_128_n_relax_4_gam_3.69e+00_D_16384_phi.csv"

phi_FD = np.genfromtxt(f"{indir_FD}/{phi_name_FD}", delimiter=",")

phi_FD = phi_FD.reshape(-1, 128, 128).transpose(1, 2, 0)


# %%
save = True
dt = 5.5e-7
timepoints = [0, 25, 50, 75, 100, 200]
fig, axs = plt.subplots(2, len(timepoints), figsize=(9 * len(timepoints), 15))
for i, timepoint in enumerate(timepoints):
    # for timepoint in [100]:
    normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
    s = sns.heatmap(
        phi_FD[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        ax=axs[0, i],
        xticklabels=20,
        yticklabels=20,
        norm=normalize_phis,
    )
    axs[0, i].set_title(f"Time= {timepoint*dt}")

    s = sns.heatmap(
        phi_MG[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        ax=axs[1, i],
        xticklabels=20,
        yticklabels=20,
        norm=normalize_phis,
    )
    # axs[1, i].set_title(f"Time= {timepoint*dt}")
    # axs[1].collections[0].cmap.set_bad(".5")

plt.suptitle(f"FD (Top) vs NMG (Bottom); dt = {dt}")
plt.tight_layout()

if save:
    plt.savefig(f"{outdir}/FD_vs_NMG_dt5.5e-7_v2.png")
    plt.close()
else:
    plt.show()

# %%
dt = 5.5e-6

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"
phi_name = "MG_2000_dt_5.5e-6_Nx_128_n_relax_4_eps_0.015009369912862116_phi.txt"

phi_MG_large_dt = np.genfromtxt(
    f"{indir}/{phi_name}",
)

phi_MG_large_dt = phi_MG_large_dt.reshape(-1, 128, 128).transpose(1, 2, 0)
timepoints = [10, 20, 100, 1000, 2000]
fig, axs = plt.subplots(1, len(timepoints), figsize=(9 * len(timepoints), 8))
for i, timepoint in enumerate(timepoints):
    # for timepoint in [100]:
    normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
    s = sns.heatmap(
        phi_MG_large_dt[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        ax=axs[i],
        xticklabels=20,
        yticklabels=20,
        norm=normalize_phis,
    )
    axs[i].set_title("Time={:.2e}".format(timepoint * dt))
plt.tight_layout()
plt.savefig(f"{outdir}/NMG_dt5.5e-6_v2.png")
# plt.show()
# %%
