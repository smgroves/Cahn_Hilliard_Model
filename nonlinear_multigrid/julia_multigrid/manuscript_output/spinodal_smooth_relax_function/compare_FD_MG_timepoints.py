# %% FIGURE 1 COMPARISONS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function"
indir_MG = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_Julia"
# %%
phi_name_MG = "NMG_MATLAB_2000_dt_5.50e-06_Nx_128_n_relax_4_phi.csv"
# indir_MG = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"
# phi_name_MG = "MG_2000_dt_5.5e-6_Nx_128_n_relax_4_eps_0.015009369912862116_phi.txt"

phi_MG = np.genfromtxt(
    f"{indir_MG}/{phi_name_MG}",
)

phi_MG = phi_MG.reshape(-1, 128, 128).transpose(1, 2, 0)

# %%
indir_FD = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output"
phi_name_FD = "FD_MATLAB_2000_dt_5.50e-06_Nx_128_n_relax_4_phi.csv"
# indir_FD = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function"
# phi_name_FD = "FD_2000_dt_5.50e-08_Nx_128_n_relax_4_gam_3.69e+00_D_16384_phi.csv"

phi_FD = np.genfromtxt(f"{indir_FD}/{phi_name_FD}", delimiter=",")

phi_FD = phi_FD.reshape(-1, 128, 128).transpose(1, 2, 0)

# %%
indir_SAV = "/Users/smgroves/Documents/GitHub/CHsolvers_package/output/output_MATLAB"

phi_name_SAV = "SAV_MATLAB_2000_dt_5.50e-06_Nx_128_n_relax_4_phi.csv"
# indir_SAV = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/SAV/output/spinodal_smooth_relax_function"
# phi_name_SAV = "SAV_MATLAB_2000_dt_5.50e-06_Nx_128_n_relax_4_phi.csv"

phi_SAV = np.genfromtxt(f"{indir_SAV}/{phi_name_SAV}", delimiter=",")

phi_SAV = phi_SAV.reshape(-1, 128, 128).transpose(1, 2, 0)

# %% save individual plots
indir_MG = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"

# timepoints = [0, 10, 20, 100, 1000, 2000]
timepoints = [200, 400]
dt_out = 1
dt = 5.5e-6
for timepoint in timepoints:
    normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
    s = sns.heatmap(
        phi_MG[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        norm=normalize_phis,
        cbar=False,
        linewidths=0.0,
    )
    plt.xticks(ticks=[], labels=[])
    plt.yticks(ticks=[], labels=[])
    # plt.title(f"Time= {timepoint*dt}")
    plt.tight_layout()
    plt.savefig(
        f"{indir_MG}/MG_2000_dt_5.5e-6_t_{timepoint*dt*dt_out:.2e}.png",
        bbox_inches="tight",
        pad_inches=0,
        dpi=300,
    )
    plt.close()
    # plt.show()

# %% zoom in to timepoint 10:
timepoint = 1
dt_out = 1
dt = 5.5e-6
normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
s = sns.heatmap(
    phi_SAV[65:80, 50:65, timepoint - 1],
    square=True,
    cmap=cm.RdBu_r,
    norm=normalize_phis,
    cbar=False,
    linewidths=0.0,
)
plt.xticks(ticks=[], labels=[])
plt.yticks(ticks=[], labels=[])
# plt.title(f"Time= {timepoint*dt}")
plt.tight_layout()
plt.savefig(
    f"{indir_SAV}/zoomed_65-80v50-65_SAV_2000_dt_5.5e-6_t_{timepoint*dt*dt_out:.2e}.png",
    bbox_inches="tight",
    pad_inches=0,
    dpi=300,
)
plt.close()
# %%
indir_FD = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function"
phi_name_FD_big = "FD_200_dt_5.50e-07_Nx_128_n_relax_4_gam_3.69e+00_D_16384_phi.csv"

phi_FD_big = np.genfromtxt(f"{indir_FD}/{phi_name_FD_big}", delimiter=",")

phi_FD_big = phi_FD_big.reshape(-1, 128, 128).transpose(1, 2, 0)

# %% save individual plots

timepoints = [0, 500, 750, 1000, 2000]
dt_out = 1
dt = 5.5e-8
for timepoint in timepoints:
    normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
    s = sns.heatmap(
        phi_FD[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        norm=normalize_phis,
        cbar=False,
        linewidths=0.0,
    )
    s.collections[0].cmap.set_bad("grey")
    plt.xticks(ticks=[], labels=[])
    plt.yticks(ticks=[], labels=[])
    # plt.title(f"Time= {timepoint*dt}")
    plt.tight_layout()
    plt.savefig(
        f"{indir_FD}/FD_2000_dt_{dt}_t_{timepoint*dt*dt_out:.2e}.png",
        bbox_inches="tight",
        pad_inches=0,
        dpi=300,
    )
    plt.close()
# %% zoom in to timepoint 75:
timepoint = 75
dt_out = 1
dt = 5.5e-7
normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
s = sns.heatmap(
    phi_FD_big[65:80, 50:65, timepoint - 1],
    square=True,
    cmap=cm.RdBu_r,
    norm=normalize_phis,
    cbar=False,
    linewidths=0.0,
)
plt.xticks(ticks=[], labels=[])
plt.yticks(ticks=[], labels=[])
# plt.title(f"Time= {timepoint*dt}")
plt.tight_layout()
plt.savefig(
    f"{indir_FD}/zoomed_65-80v50-65_FD_200_dt_5.5e-7_t_{timepoint*dt*dt_out:.2e}.png",
    bbox_inches="tight",
    pad_inches=0,
    dpi=300,
)
plt.close()
# %% save individual plots

timepoints = [1, 2, 40, 100, 200]
# timepoints = [40]
dt_out = 10
dt = 5.5e-6
for timepoint in timepoints:
    normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
    s = sns.heatmap(
        # the first one saved by SAV is timestep = 0, so I updated this to no longer subtract 1
        phi_SAV[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        norm=normalize_phis,
        cbar=False,
        linewidths=0.0,
    )
    plt.xticks(ticks=[], labels=[])
    plt.yticks(ticks=[], labels=[])
    # plt.title(f"Time= {timepoint*dt}")
    plt.tight_layout()
    plt.savefig(
        f"{indir_SAV}/SAV_2000_dt_5.5e-6_t_{timepoint*dt*dt_out:.2e}_neumann.png",
        bbox_inches="tight",
        pad_inches=0,
        dpi=300,
    )
    plt.close()
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
dt = 5.5e-8

timepoints = [500, 750, 1000, 2000]
fig, axs = plt.subplots(1, len(timepoints), figsize=(9 * len(timepoints), 8))
for i, timepoint in enumerate(timepoints):
    # for timepoint in [100]:
    normalize_phis = mcolors.TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
    s = sns.heatmap(
        phi_FD[:, :, timepoint],
        square=True,
        cmap=cm.RdBu_r,
        ax=axs[i],
        xticklabels=20,
        yticklabels=20,
        norm=normalize_phis,
    )
    axs[i].set_title("Time={:.2e}".format(timepoint * dt))
plt.tight_layout()
plt.savefig(f"{outdir}/FD_dt5.5e-8_v2.png")
# plt.show()


# %% comparing energy and mass
# FIGURE S1
def get_energy(indir, phi_name, dt, dt_out, variable="energy", suffix="csv", title=""):
    e_name = "_".join(phi_name.split("_")[:-1])
    e_name = e_name + f"_{variable}.{suffix}"
    print(e_name)
    e = pd.read_csv(f"{indir}/{e_name}", header=None, index_col=None)
    l = [x * dt * dt_out for x in range(len(e.iloc[:, 0]))]
    e.columns = [title]
    e["time"] = l
    return e


# %%

e_NMG = get_energy(
    indir_MG, phi_name_MG, dt=5.5e-6, dt_out=1, suffix="txt", title="NMG"
)
e_FD = get_energy(indir_FD, phi_name_FD, dt=5.5e-8, dt_out=1, title="FD")
e_FD_big = get_energy(indir_FD, phi_name_FD_big,
                      dt=5.5e-7, dt_out=1, title="FD_big")

e_SAV = get_energy(indir_SAV, phi_name_SAV, dt=5.5e-6, dt_out=10, title="SAV")

# %%
for i, e in enumerate([e_NMG, e_FD, e_SAV, e_FD_big]):
    plt.figure(figsize=(6, 2))
    ax = sns.lineplot(x=e["time"], y=e.iloc[:, 0], c=sns.color_palette()[i])
    plt.ylabel("Normalized Energy")
    plt.xlabel(r"Time ($t_{char}$)")
    plt.ylim(0.23, 1.05)
    # ax.set(xscale="log")
    # ax.set(xticks = np.arange(0, 0.012, 0.002))
    # plt.title("Normalized (Modified) Energy for Spinodal Decomposition")
    if e.iloc[:, 0].name == "FD":
        title = "FD, dt = 5.5e-8"
        plt.title(title)
        plt.xlim(-0.001, 0.012)

    elif e.iloc[:, 0].name == "FD_big":
        title = "FD, dt = 5.5e-7"
        plt.title(title)
    else:
        title = f"{e.iloc[:, 0].name}, dt=5.5e-6"
        plt.title(title)
        plt.xlim(-0.001, 0.012)
    plt.ticklabel_format(style="plain", axis="x")

    # plt.show()

    plt.savefig(f"{outdir}/{title}_energy.pdf")

# %%
for e in [e_NMG, e_FD, e_SAV]:
    sns.lineplot(x=e["time"], y=e.iloc[:, 0], label=e.columns[0], alpha=0.6)
plt.xlim(0.01, 0.0111)
plt.ylim(0.255, 0.27)
plt.ylabel("")
plt.xlabel("")
# plt.title("Normalized (Modified) Energy for Spinodal Decomposition")
plt.savefig(f"{outdir}/FD_NMG_SAV_energy_zoom_end.pdf")
# %%

m_NMG = get_energy(
    indir_MG,
    phi_name_MG,
    dt=5.5e-6,
    dt_out=1,
    variable="mass",
    suffix="txt",
    title="NMG",
)
m_NMG["NMG"] = m_NMG["NMG"] - m_NMG["NMG"].iloc[0]

m_FD = get_energy(
    indir_FD, phi_name_FD, dt=5.5e-8, variable="mass", dt_out=1, title="FD"
)
m_FD["FD"] = m_FD["FD"] - m_FD["FD"].iloc[0]


m_FD_big = get_energy(
    indir_FD, phi_name_FD_big, dt=5.5e-7, variable="mass", dt_out=1, title="FD_big"
)
m_FD_big["FD_big"] = m_FD_big["FD_big"] - m_FD_big["FD_big"].iloc[0]

m_SAV = get_energy(
    indir_SAV, phi_name_SAV, dt=5.5e-6, variable="mass", dt_out=10, title="SAV"
)
# %%
for m in [m_NMG, m_FD, m_SAV]:
    sns.lineplot(x=m["time"], y=m.iloc[:, 0], label=m.columns[0], alpha=0.6)
plt.ylabel("Normalized Average Mass")
plt.xlabel(r"Time ($t_{char}$)")
plt.title("Normalized Average Mass for Spinodal Decomposition")
plt.savefig(f"{outdir}/FD_NMG_SAV_mass.pdf")
# %%
for i, e in enumerate([m_NMG, m_FD, m_SAV, m_FD_big]):
    plt.figure(figsize=(6, 2))
    ax = sns.lineplot(x=e["time"], y=e.iloc[:, 0], c=sns.color_palette()[i])
    plt.ylabel("Normalized Mass")
    plt.xlabel(r"Time ($t_{char}$)")
    # ax.set(xscale="log")
    # ax.set(xticks = np.arange(0, 0.012, 0.002))
    # plt.title("Normalized (Modified) Energy for Spinodal Decomposition")
    if e.iloc[:, 0].name == "FD":
        title = "FD, dt = 5.5e-8"
        plt.title(title)
        # plt.xlim(-0.001, 0.012)
        plt.ylim(-1e-6, 1e-6)

    elif e.iloc[:, 0].name == "FD_big":
        title = "FD, dt = 5.5e-7"
        plt.title(title)
    else:
        title = f"{e.iloc[:, 0].name}, dt=5.5e-6"
        plt.title(title)
        plt.xlim(-0.001, 0.012)
        plt.ylim(-1e-6, 1e-6)

    # plt.title(f"{e.iloc[:, 0].name}")
    plt.ticklabel_format(style="plain", axis="x")

    plt.savefig(f"{outdir}/{title}_mass_rescaled_Axis.pdf")

# %%
