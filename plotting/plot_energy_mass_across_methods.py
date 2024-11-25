# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting/comparison_sim_results_across_methods"


# %% read_in_files function
def read_in_files(
    path, suffix, timesteps, dt, solver_name, keepevery=10, filetype="csv", plus_one=1
):
    e_df = pd.read_csv(
        f"{path}/{solver_name}_{timesteps}_dt_{dt}_{suffix}_energy.{filetype}",
        header=None,
        index_col=None,
    )
    e_df.columns = ["Energy"]

    m_df = pd.read_csv(
        f"{path}/{solver_name}_{timesteps}_dt_{dt}_{suffix}_mass.{filetype}",
        header=None,
        index_col=None,
    )
    m_df.columns = ["Mass"]

    df = pd.concat([e_df, m_df], axis=1)
    df["Time"] = np.linspace(
        0, float(dt) * (float(timesteps) + plus_one), int(timesteps) + plus_one
    )
    df["Solver"] = solver_name + ", dt=" + dt
    df = df[df.index % keepevery == 0]
    return df


#####################################################################
# %% Large and small droplets
#####################################################################
FD_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/large_and_small_droplets"
FD_suffix = "Nx_128_gam_3.69e+00_D_16384_r1_20_r2_30_space_10"
name = "FD"
# FD_timesteps = "1500000"
# FD_dt = "1.00e-06"
# find all files with this suffix, pull out timesteps and dt
os.chdir(FD_path)
z = []
for file in glob.glob(f"*{FD_suffix}_energy.csv"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])
all_FD = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    if (it == "5000000") and (dt == "1.00e-06"):
        continue
    print(it, dt)

    FD_df = read_in_files(
        FD_path,
        FD_suffix,
        it,
        dt,
        name,
        # keepevery=(round(int(it) / 1000) if int(it) > 10000 else 10),
        keepevery=1,
        filetype="csv",
        plus_one=1,
    )
    FD_df["Mass"] = FD_df["Mass"] - FD_df.iloc[0]["Mass"]

    all_FD = pd.concat([all_FD, FD_df], ignore_index=True)

# %%
MG_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/output"
MG_suffix = "Nx_128_eps_0.015009369912862116_r1_20_r2_30_space_10"
name = "MG"
filetype = "txt"
# find all files with this suffix, pull out timesteps and dt
os.chdir(MG_path)
z = []
for file in glob.glob(f"*{MG_suffix}_energy.{filetype}"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])

all_MG = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    print(it, dt)
    MG_df = read_in_files(
        MG_path,
        MG_suffix,
        it,
        dt,
        name,
        # keepevery=(round(int(it) / 1000) if int(it) > 10000 else 10),
        keepevery=1,
        filetype="txt",
        plus_one=0,
    )
    MG_df["Mass"] = MG_df["Mass"] - MG_df.iloc[0]["Mass"]

    all_MG = pd.concat([all_MG, MG_df], ignore_index=True)

# %%
for i, r in all_MG.iterrows():
    if r["Solver"] == "MG, dt=0.001":
        all_MG.loc[i, "Solver"] = "MG, dt=1.00e-03"
    elif r["Solver"] == "MG, dt=0.0001":
        all_MG.loc[i, "Solver"] = "MG, dt=1.00e-04"

# %%
all_data = pd.concat([all_FD, all_MG], ignore_index=True)


# %%
all_data = all_data.sort_values(["Solver", "Time"])

sns.lineplot(
    all_data,
    x="Time",
    y="Energy",
    hue="Solver",
)
plt.title("Ostwald Ripening of Two Droplets")
plt.ylabel("Normalized Discrete Energy")
plt.savefig(f"{outdir}/two_droplets_totaltime_5_FD_MG_energy.pdf")
# plt.show()
plt.close()
sns.lineplot(
    all_data,
    x="Time",
    y="Mass",
    hue="Solver",
)
plt.title("Ostwald Ripening of Two Droplets")
plt.ylabel("Deviation from Initial Average Mass (M(t)-M(0))")
plt.savefig(f"{outdir}/two_droplets_totaltime_5_FD_MG_mass.pdf")
# plt.show()

# %%
# run with keepevery = 1
short = all_data.loc[all_data["Time"] <= 0.01]

# %%
short = short.sort_values(["Solver", "Time"])
dash_df = {
    "FD, dt=1.00e-05": (3, 2),
    "FD, dt=1.00e-04": (4, 1),
    "FD, dt=1.00e-03": (2, 3),
    "FD, dt=1.00e-06": "",
    "MG, dt=6.103515625e-5": "",
    "MG, dt=1.00e-04": "",
    "MG, dt=1.00e-03": "",
}
sns.lineplot(short, x="Time", y="Energy", hue="Solver", style="Solver", dashes=dash_df)
plt.ylim(0, 1.2)
plt.title("Ostwald Ripening of Two Droplets")
plt.ylabel("Normalized Discrete Energy")
plt.savefig(f"{outdir}/two_droplets_totaltime_.05_FD_MG.pdf")

# plt.show()
#####################################################################
# %% spinodal smoothed
#####################################################################
FD_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_smooth_relax_function"
FD_suffix = "Nx_128_n_relax_2_gam_3.69e+00_D_16384"
FD_timesteps = "200000"
FD_dt = "6.10e-08"
FD_df = read_in_files(
    FD_path,
    FD_suffix,
    FD_timesteps,
    FD_dt,
    "FD",
    keepevery=1000,
    filetype="csv",
    plus_one=1,
)

# %%
MG_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"
MG_suffix = "Nx_128_tol_1.0e-6_n_relax_2"
MG_timesteps = "9600"
MG_dt = "6.25e-6"
MG_df = read_in_files(
    MG_path,
    MG_suffix,
    MG_timesteps,
    MG_dt,
    "MG",
    keepevery=10,
    filetype="txt",
    plus_one=0,
)
MG_df["Mass"] = MG_df["Mass"] / MG_df.iloc[0]["Mass"]
all_data = pd.concat([FD_df, MG_df], ignore_index=True)
sns.lineplot(
    all_data,
    x="Time",
    y="Energy",
    hue="Solver",
)
plt.title("Spinodal pre-smoothened IC (n_relax = 2)")
plt.ylabel("Normalized Discrete Energy")
plt.show()

# %%
sns.lineplot(
    all_data,
    x="Time",
    y="Mass",
    hue="Solver",
)
plt.title("Spinodal pre-smoothened IC (n_relax = 2)")
plt.ylabel("Normalized Average Mass")
plt.show()
# %%
# %% Checkerboard
#####################################################################
FD_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/checkerboard_IC"
FD_suffix = "Nx_128_gam_3.69e+00_D_16384_grid_8"
name = "FD"
# FD_timesteps = "1500000"
# FD_dt = "1.00e-06"
# find all files with this suffix, pull out timesteps and dt
os.chdir(FD_path)
z = []
for file in glob.glob(f"*{FD_suffix}_energy.csv"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])
all_FD = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    print(it, dt)
    FD_df = read_in_files(
        FD_path,
        FD_suffix,
        it,
        dt,
        name,
        keepevery=(round(int(it) / 100) if int(it) > 10000 else 1),
        # keepevery = 1,
        filetype="csv",
        plus_one=1,
    )
    FD_df["Mass"] = FD_df["Mass"] - FD_df.iloc[0]["Mass"]

    all_FD = pd.concat([all_FD, FD_df], ignore_index=True)

# %%
MG_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/checkerboard_IC/output"
MG_suffix = "Nx_128_eps_0.015009369912862116_grid_8"
name = "MG"
filetype = "txt"
# find all files with this suffix, pull out timesteps and dt
os.chdir(MG_path)
z = []
for file in glob.glob(f"*{MG_suffix}_energy.{filetype}"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])

all_MG = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    print(it, dt)
    MG_df = read_in_files(
        MG_path,
        MG_suffix,
        it,
        dt,
        name,
        keepevery=(round(int(it) / 100) if int(it) > 10000 else 1),
        # keepevery = 1,
        filetype="txt",
        plus_one=0,
    )
    MG_df["Mass"] = MG_df["Mass"] - MG_df.iloc[0]["Mass"]

    all_MG = pd.concat([all_MG, MG_df], ignore_index=True)
# %%
all_data = pd.concat([all_FD, all_MG], ignore_index=True)


# %%
for i, r in all_data.iterrows():
    if r["Solver"] == "MG, dt=0.0006103515625":
        all_data.loc[i, "Solver"] = "MG, dt=6.10e-04"
    elif r["Solver"] == "MG, dt=6.103515625e-5":
        all_data.loc[i, "Solver"] = "MG, dt=6.10e-05"
    elif r["Solver"] == "MG, dt=6.103515625e-6":
        all_data.loc[i, "Solver"] = "MG, dt=6.10e-06"

# %%
all_data = all_data.sort_values(["Solver", "Time"])

plt.rcParams["patch.edgecolor"] = "none"

sns.lineplot(
    all_data,
    x="Time",
    y="Energy",
    hue="Solver",
    # marker="o",
    # markersize=6,
    # markeredgecolor="None",
)
plt.ylim(0, 1.2)
plt.title("Checkerboard IC - Energy")
plt.ylabel("Normalized Discrete Energy")
plt.savefig(f"{outdir}/checkerboard_FD_MG_energy.pdf")
# plt.show()
plt.close()

# %%
small = all_data.loc[all_data["Solver"] != "FD, dt=6.10e-05"]
sns.lineplot(small, x="Time", y="Mass", hue="Solver", palette=sns.color_palette()[1:6])
plt.title("Checkerboard IC - Mass")
plt.ylabel("Deviation from Initial Average Mass (M(t)-M(0))")
plt.savefig(f"{outdir}/checkerboard_FD_MG_mass_noFD6.1e-5.pdf")
# plt.show()

#####################################################################
# %% spinodal +1/-1
#####################################################################
FD_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_+1_-1"
FD_suffix = "Nx_128_gam_3.69e+00_D_16384"
name = "FD"
# FD_timesteps = "1500000"
# FD_dt = "1.00e-06"
# find all files with this suffix, pull out timesteps and dt
os.chdir(FD_path)
z = []
for file in glob.glob(f"*{FD_suffix}_energy.csv"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])
all_FD = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    if it == "983040":
        continue
    print(it, dt)
    FD_df = read_in_files(
        FD_path,
        FD_suffix,
        it,
        dt,
        name,
        keepevery=(round(int(it) / 100) if int(it) > 10000 else 10),
        # keepevery = 1,
        filetype="csv",
        plus_one=1,
    )
    FD_df["Mass"] = FD_df["Mass"] - FD_df.iloc[0]["Mass"]

    all_FD = pd.concat([all_FD, FD_df], ignore_index=True)

# %%
MG_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_+1_-1_IC/output"
MG_suffix = "Nx_128_eps_0.015_tol_1.0e-6"
name = "MG"
filetype = "txt"
# find all files with this suffix, pull out timesteps and dt
os.chdir(MG_path)
z = []
for file in glob.glob(f"*{MG_suffix}_energy.{filetype}"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])

all_MG = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    print(it, dt)
    MG_df = read_in_files(
        MG_path,
        MG_suffix,
        it,
        dt,
        name,
        keepevery=(round(int(it) / 100) if int(it) > 10000 else 10),
        # keepevery = 1,
        filetype="txt",
        plus_one=0,
    )
    MG_df["Mass"] = MG_df["Mass"] - MG_df.iloc[0]["Mass"]

    all_MG = pd.concat([all_MG, MG_df], ignore_index=True)
# %%
all_data = pd.concat([all_FD, all_MG], ignore_index=True)

# %%
all_data = all_data.sort_values(["Solver", "Time"])

plt.rcParams["patch.edgecolor"] = "none"

sns.lineplot(
    all_data,
    x="Time",
    y="Energy",
    hue="Solver",
    # marker="o",
    # markersize=6,
    # markeredgecolor="None",
    palette="tab20",
)
plt.ylim(0, 0.1)
plt.title("Spinodal Decomposition (+1/-1 IC) - Energy")
plt.ylabel("Normalized Discrete Energy")
# plt.savefig(f"{outdir}/SD_+1_-1_FD_MG_energy_zoomed.png")
plt.show()
plt.close()
# %%
sns.lineplot(
    all_data, x="Time", y="Mass", hue="Solver"
)  # , palette=sns.color_palette()[1:6]
# )
plt.title("Spinodal Decomposition (+1/-1 IC) - Mass")
plt.ylabel("Deviation from Initial Average Mass (M(t)-M(0))")
plt.savefig(f"{outdir}/SD_+1_-1_FD_MG_mass.pdf")
# plt.show()

# %%
#####################################################################
# %% spinodal t = 0.01 from MG
#####################################################################
FD_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/finite_difference_method/output/spinodal_MG_timepoint_IC/t=0.01"
FD_suffix = "Nx_128_gam_3.69e+00_D_16384"
name = "FD"
# FD_timesteps = "1500000"
# FD_dt = "1.00e-06"
# find all files with this suffix, pull out timesteps and dt
os.chdir(FD_path)
z = []
for file in glob.glob(f"*{FD_suffix}_energy.csv"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])
all_FD = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    print(it, dt)
    FD_df = read_in_files(
        FD_path,
        FD_suffix,
        it,
        dt,
        name,
        # keepevery=(round(int(it) / 100) if int(it) > 10000 else 10),
        keepevery=1,
        filetype="csv",
        plus_one=1,
    )
    FD_df["Mass"] = FD_df["Mass"] - FD_df.iloc[0]["Mass"]

    all_FD = pd.concat([all_FD, FD_df], ignore_index=True)

# %%
MG_path = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_MG_timepoint_IC/output"
MG_suffix = "Nx_128_eps_0.015009369912862116"
name = "MG"
filetype = "txt"
# find all files with this suffix, pull out timesteps and dt
os.chdir(MG_path)
z = []
for file in glob.glob(f"*{MG_suffix}_energy.{filetype}"):
    it_tmp = file.split("_")[1]
    dt_tmp = file.split("_")[3]
    z.append([it_tmp, dt_tmp])

all_MG = pd.DataFrame(columns=["Time", "Solver", "Mass", "Energy"])
# loop through them and add to all_FD
for it, dt in z:
    print(it, dt)
    MG_df = read_in_files(
        MG_path,
        MG_suffix,
        it,
        dt,
        name,
        # keepevery=(round(int(it) / 100) if int(it) > 10000 else 10),
        keepevery=1,
        filetype="txt",
        plus_one=0,
    )
    MG_df["Mass"] = MG_df["Mass"] - MG_df.iloc[0]["Mass"]

    all_MG = pd.concat([all_MG, MG_df], ignore_index=True)
# %%
all_data = pd.concat([all_FD, all_MG], ignore_index=True)

# %%
all_data = all_data.sort_values(["Solver", "Time"])

plt.rcParams["patch.edgecolor"] = "none"

sns.lineplot(
    all_data,
    x="Time",
    y="Energy",
    hue="Solver",
    # marker="o",
    # markersize=6,
    # markeredgecolor="None",
    palette="tab10",
)
plt.ylim(0.5, 1.4)
plt.title("Spinodal Decomposition (t = 0.01 MG IC) - Energy")
plt.ylabel("Normalized Discrete Energy")
plt.savefig(f"{outdir}/SD_t=0.01_IC_FD_MG_energy.png")
# plt.show()
plt.close()
# %%
small = all_data.loc[
    ~(
        all_data["Solver"].isin(
            ["FD, dt=1.00e-04", "FD, dt=1.25e-05", "FD, dt=2.50e-05", "FD, dt=5.00e-05"]
        )
    )
]
sns.lineplot(
    small, x="Time", y="Mass", hue="Solver"
)  # , palette=sns.color_palette()[1:6]
# )
plt.title("Spinodal Decomposition (t = 0.01 MG IC) - Mass")
plt.ylabel("Deviation from Initial Average Mass (M(t)-M(0))")
plt.savefig(f"{outdir}/SD_t=0.01_IC_FD_MG_mass.pdf")
# plt.show()


# %%
