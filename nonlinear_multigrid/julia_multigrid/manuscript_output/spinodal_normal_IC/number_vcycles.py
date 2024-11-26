# %%
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from matplotlib.backends.backend_pdf import PdfPages


def calculate_vcycles(nx, tol, dt, indir, one_timestep=True):
    num_cycles = []
    count = 0
    highest_r = 100000
    num_timesteps = 0
    with open(f"{indir}/residual_{nx}_{tol}_dt_{dt}_mean_0_sd_0.2.txt") as file:
        if one_timestep:
            for count, line in enumerate(file):
                pass
            ave_num_cycles = count + 1
        else:
            for l, line in enumerate(file):
                r = float(line.rstrip())
                if r < highest_r:
                    highest_r = r
                    count += 1
                else:
                    num_cycles.append(count)
                    highest_r = 100000
                    count = 1
                    num_timesteps += 1
                # if num_timesteps == 50:
                # break
            ave_num_cycles = math.ceil(np.mean(num_cycles))
    return ave_num_cycles


indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output/dt_vs_tol_num_vcycles"

table = pd.DataFrame(columns=["nx", "tol", "dt", "num_vcycles"])

nx = 64
tols = np.array(
    ["0.001", "0.0001", "1.0e-5", "1.0e-6", "1.0e-7", "1.0e-8", "1.0e-9", "1.0e-10"]
)
h = 1 / nx
for tol in tols:
    for dt_multiplier in [
        "0",
        "-1",
        "-2",
        "-3",
        "-4",
        "-5",
        "-6",
        "-7",
        "-8",
        "-9",
        "-10",
        "-11",
        "-12",
    ]:

        try:
            if dt_multiplier == "-3":
                dt = "0.006103515625"
            elif dt_multiplier == "-4":
                dt = "0.0006103515625"
            elif dt_multiplier == "-2":
                dt = "0.06103515625"
            elif dt_multiplier == "-1":
                dt = "0.6103515625"
            elif dt_multiplier == "0":
                dt = "6.103515625"
            else:
                if nx == 64:
                    dt = "2.44140625e" + dt_multiplier
                elif nx == 256:
                    dt = "1.52587890625e" + dt_multiplier
                else:
                    dt = "6.103515625e" + dt_multiplier

            ans = calculate_vcycles(nx, tol, dt, indir, one_timestep=True)
            table = pd.concat(
                [
                    table,
                    pd.DataFrame(
                        {
                            "nx": [nx],
                            "tol": [float(tol)],
                            "dt": [float(dt)],
                            "num_vcycles": [ans],
                        }
                    ),
                ],
                ignore_index=True,
            )
        except FileNotFoundError:
            pass

# print(table)

# print(pd.pivot_table(table, columns=["dt"], index="tol", aggfunc="mean"))

print(table)
table = table.drop("nx", axis=1)
final_table = pd.pivot(table, columns="dt", index="tol")
final_table = final_table.sort_index(ascending=False)
final_table = final_table.sort_index(ascending=True, axis=1)
final_table = final_table.droplevel(0, axis=1)

# %% rename columns

dts = final_table.columns
new_names = []
for d in dts:
    factor = d / (1 / nx**2)
    new_names.append(r"{:.0e}h^2".format((factor)))
print(new_names)
final_table.columns = new_names
# %%


def _draw_as_table(df, pagesize):
    alternating_colors = [
        ["white"] * len(df.columns),
        ["lightgray"] * len(df.columns),
    ] * len(df)
    alternating_colors = alternating_colors[: len(df)]
    fig, ax = plt.subplots(figsize=pagesize)
    ax.axis("tight")
    ax.axis("off")
    the_table = ax.table(
        cellText=df.values,
        rowLabels=df.index,
        colLabels=df.columns,
        rowColours=["lightblue"] * len(df),
        colColours=["lightblue"] * len(df.columns),
        cellColours=alternating_colors,
        loc="center",
    )
    return fig


def dataframe_to_pdf(df, filename, numpages=(1, 1), pagesize=(11, 8.5)):
    with PdfPages(filename) as pdf:
        nh, nv = numpages
        rows_per_page = len(df) // nh
        cols_per_page = len(df.columns) // nv
        for i in range(0, nh):
            for j in range(0, nv):
                page = df.iloc[
                    (i * rows_per_page) : min((i + 1) * rows_per_page, len(df)),
                    (j * cols_per_page) : min((j + 1) * cols_per_page, len(df.columns)),
                ]
                fig = _draw_as_table(page, pagesize)
                if nh > 1 or nv > 1:
                    # Add a part/page number at bottom-center of page
                    fig.text(
                        0.5,
                        0.5 / pagesize[0],
                        "Part-{}x{}: Page-{}".format(i + 1, j + 1, i * nv + j + 1),
                        ha="center",
                        fontsize=8,
                    )
                pdf.savefig(fig, bbox_inches="tight")

                plt.close()


final_table.to_csv(f"{indir}/dt_vs_tol_n_vcycles_spinodal_normal_nx_{nx}.csv")

dataframe_to_pdf(final_table, f"{indir}/table_julia_spinodal_normal_nx_{nx}.pdf")

# %%
