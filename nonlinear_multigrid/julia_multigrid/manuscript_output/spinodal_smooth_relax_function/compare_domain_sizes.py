import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function"

indir_small = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"
indir_large = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output_large_domain"
phi_name_small = (
    "MG_2000_dt_6.103515625e-6_Nx_128_n_relax_4_eps_0.015009369912862116_final_phi.csv"
)
phi_name_large = (
    "MG_2000_dt_0.1_Nx_128_n_relax_4_eps_1.9211993488463508_large_domain_final_phi.csv"
)
phi_final_small = np.genfromtxt(
    f"{indir_small}/{phi_name_small}",
    delimiter=",",
)
phi_final_large = np.genfromtxt(
    f"{indir_large}/{phi_name_large}",
    delimiter=",",
)

difference = phi_final_small - phi_final_large
vcenter = 0
vmin, vmax = difference.min(), difference.max()
normalize_phis = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=-1, vmax=1)
normalize_diff = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

fig, axs = plt.subplots(1, 3, figsize=(24, 6))
s = sns.heatmap(
    phi_final_large,
    square=True,
    cmap=cm.RdBu_r,
    ax=axs[0],
    xticklabels=20,
    yticklabels=20,
    norm=normalize_phis,
)
axs[0].set_title("NMG (domain = 0,128)")

s = sns.heatmap(
    phi_final_small,
    square=True,
    cmap=cm.RdBu_r,
    ax=axs[1],
    xticklabels=20,
    yticklabels=20,
    norm=normalize_phis,
)
axs[1].set_title("NMG (domain = 0,1)")


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
plt.suptitle(
    "Final phi after 2000 timesteps \nSD from relaxed IC (n_relax=4), domain=[0,1] vs [0,128]"
)
plt.tight_layout()
plt.show()
# plt.savefig(f"{outdir}/domain_comparison_0_1_vs_0_128.png")

# indir_FD = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function"
# phi_name_FD = "FD_20000_dt_6.10e-07_Nx_128_n_relax_4_gam_1.48e+01_D_16384_final_phi.csv"
# phi_final_FD = np.genfromtxt(
#     f"{indir_FD}/{phi_name_FD}",
#     delimiter=",",
# )
# difference = phi_final_large - phi_final_FD


# vcenter = 0
# vmin, vmax = difference.min(), difference.max()
# normalize_phis = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=-1, vmax=1)
# normalize_diff = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

# fig, axs = plt.subplots(1, 3, figsize=(24, 6))
# s = sns.heatmap(
#     phi_final_large,
#     square=True,
#     cmap=cm.RdBu_r,
#     ax=axs[0],
#     xticklabels=20,
#     yticklabels=20,
#     norm=normalize_phis,
# )
# axs[0].set_title("NMG")

# s = sns.heatmap(
#     phi_final_FD,
#     square=True,
#     cmap=cm.RdBu_r,
#     ax=axs[1],
#     xticklabels=20,
#     yticklabels=20,
#     norm=normalize_phis,
# )
# axs[1].set_title("FD (dt=6.10e-07)")


# s = sns.heatmap(
#     difference,
#     square=True,
#     norm=normalize_diff,
#     cmap=cm.PiYG,
#     ax=axs[2],
#     xticklabels=20,
#     yticklabels=20,
# )
# axs[2].set_title("Difference")

# plt.suptitle(
#     "Final phi after 2000 timesteps \nSD from relaxed IC (n_relax=4), FD vs NMG (domain=[0,128])"
# )
# plt.tight_layout()
# plt.show()
# # plt.savefig(f"{outdir}/method_comparison_NMG_vs_FD_0_128.png")
