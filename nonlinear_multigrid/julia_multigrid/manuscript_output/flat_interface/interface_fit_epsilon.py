import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

nx = 128
max_rows = 128
xx = np.linspace(0, 1 / 128 * (nx - 1), nx)

# indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/flat_interface/output"
# phi_name = "MG_1600_dt_6.25e-6_Nx_128_eps_0.015009369912862116_phi.txt"

# # want t = 0
# skip_header = 0
# phi_initial = np.genfromtxt(
#     f"{indir}/{phi_name}", skip_header=skip_header, max_rows=max_rows
# )

# # want t = end = 1600
# # 1600 timesteps = 160th matrix (every ten steps recorded); 0 included
# skip_header = 128 * 160
# phi_final_NMG = np.genfromtxt(
#     f"{indir}/{phi_name}", skip_header=skip_header, max_rows=max_rows
# )

# skip_header = 128 * 1600
# indir_FD = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/flat_interface"
# phi_FD_name = "FD_1600_dt_6.25e-06_Nx_128_gam_3.69e+00_D_16384_phi.csv"
# phi_final_FD = np.genfromtxt(
#     f"{indir_FD}/{phi_FD_name}",
#     skip_header=skip_header,
#     max_rows=max_rows,
#     delimiter=",",
# )
# skip_header = 0
# phi_initial_FD = np.genfromtxt(
#     f"{indir_FD}/{phi_FD_name}",
#     skip_header=0,
#     max_rows=max_rows,
#     delimiter=",",
# )
# print(phi_initial_FD)

# skip_header = 128 * 16000
# phi_FD_name = "FD_16000_dt_6.25e-07_Nx_128_gam_3.69e+00_D_16384_phi.csv"
# phi_final_FD_small = np.genfromtxt(
#     f"{indir_FD}/{phi_FD_name}",
#     skip_header=skip_header,
#     max_rows=max_rows,
#     delimiter=",",
# )

# skip_header = 128 * 160  # set ns = 1000 instead of 1, so 1600 = 160,000 timesteps
# phi_FD_name = "FD_160000_dt_6.25e-08_Nx_128_gam_3.69e+00_D_16384_phi.csv"
# phi_final_FD_smallest = np.genfromtxt(
#     f"{indir_FD}/{phi_FD_name}",
#     skip_header=skip_header,
#     max_rows=max_rows,
#     delimiter=",",
# )


# plt.plot(
#     xx,
#     phi_initial[:, int(nx / 2)],
#     label="Initial (MG and FD)",
#     linewidth=3,
# )

# plt.plot(
#     xx,
#     phi_final_NMG[:, int(nx / 2)],
#     label="Final (NMG), dt = 6.25e-6",
#     linestyle="--",
#     linewidth=3,
# )
# plt.plot(
#     xx,
#     phi_final_FD[:, int(nx / 2)],
#     label="Final (FD), dt = 6.25e-6",
#     linewidth=3,
# )
# plt.plot(
#     xx,
#     phi_final_FD_small[:, int(nx / 2)],
#     label="Final (FD), dt = 6.25e-7",
#     linestyle="--",
#     linewidth=3,
# )
# plt.plot(
#     xx,
#     phi_final_FD_smallest[:, int(nx / 2)],
#     label="Final (FD), dt = 6.25e-8",
#     linestyle="dotted",
#     linewidth=3,
# )
# M = 8
# gam = M * (1 / nx) / (2 * np.sqrt(2) * np.arctanh(0.9))
# phi_theory = np.tanh((xx - 0.5) / (np.sqrt(2) * gam))

# plt.plot(
#     xx,
#     phi_theory,
#     label="Theory for m = 8",
#     linestyle="solid",
#     color="k",
#     linewidth=1,
#     alpha=0.5,
#     marker="o",
# )

# M = 4
# gam = M * (1 / nx) / (2 * np.sqrt(2) * np.arctanh(0.9))
# phi_theory_4 = np.tanh((xx - 0.5) / (np.sqrt(2) * gam))

# plt.plot(
#     xx,
#     phi_theory_4,
#     label="Theory for m = 4",
#     linestyle="solid",
#     color="k",
#     linewidth=1,
#     alpha=0.5,
#     marker="^",
# )

# plt.axvline(x=0.5, linestyle="--", color="lightgrey")
# plt.legend()
# plt.title(
#     "Initial and Final Phi (1D profile) for Flat Interface \n total time = 0.01, m = 8"
# )


# plt.xlim(0.4, 0.6)
# plt.show()
# # plt.savefig(f"{indir}/interface_profile_multiple_dt_zoomed_+theory.pdf")
# plt.close()

########################################################################################
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/flat_interface"
phi_name = "FD_16000_dt_6.25e-07_Nx_128_gam_3.69e+00_D_16384_phi.csv"
# want t = 0
skip_header = 0
phi_initial = np.genfromtxt(
    f"{indir}/{phi_name}",
    skip_header=skip_header,
    max_rows=max_rows,
    delimiter=",",
)

# want t = end = 1600
# 1600 timesteps = 160th matrix (every ten steps recorded); 0 included
skip_header = 128 * 1600
phi_final_FD = np.genfromtxt(
    f"{indir}/{phi_name}",
    skip_header=skip_header,
    max_rows=max_rows,
    delimiter=",",
)
plt.plot(
    xx,
    phi_initial[:, int(nx / 2)],
    label="Initial (FD)",
    linewidth=3,
)

plt.plot(
    xx,
    phi_final_FD[:, int(nx / 2)],
    label="Final (FD), dt = 6.25e-7",
    linestyle="--",
    linewidth=3,
)

M = 8
gam = M * (1 / nx) / (2 * np.sqrt(2) * np.arctanh(0.9))
phi_theory = np.tanh((xx - 0.5) / (np.sqrt(2) * gam))

plt.plot(
    xx,
    phi_theory,
    label="Theory for m = 8",
    linestyle="solid",
    color="k",
    linewidth=1,
    alpha=0.5,
    marker="o",
)

M = 4
gam = M * (1 / nx) / (2 * np.sqrt(2) * np.arctanh(0.9))
phi_theory_4 = np.tanh((xx - 0.5) / (np.sqrt(2) * gam))

plt.plot(
    xx,
    phi_theory_4,
    label="Theory for m = 4",
    linestyle="solid",
    color="k",
    linewidth=1,
    alpha=0.5,
    marker="^",
)

plt.axvline(x=0.5, linestyle="--", color="lightgrey")
plt.xlim(0.4, 0.6)
plt.legend()
plt.title(
    "Initial and Final Phi (1D profile) for Flat Interface \n total time = 0.01, m = 8"
)
# plt.show()
plt.savefig(f"{indir}/interface_profile_m=8_+theory_fixed_laplace.pdf")
