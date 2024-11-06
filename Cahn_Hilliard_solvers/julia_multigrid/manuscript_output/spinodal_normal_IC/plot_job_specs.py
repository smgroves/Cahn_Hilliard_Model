import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers"
specs = pd.read_csv(f"{indir}/Job_specs.csv")
specs = specs.loc[
    (
        specs["name"].isin(
            [
                "spinodal_normal_IC",
                # "spinodal_normal_IC_no_printing",
                "spinodal_normal_IC_printing",
            ]
        )
    )
    & (specs["timesteps"].isin([4800, 2400, 1200, 600, 300]))
]

for i, r in specs.iterrows():
    if r["timesteps"] == int(4800):
        if r["dt"] != float(6.25e-06):
            specs.drop(i)

sns.lineplot(
    data=specs,
    x="nx",
    y="time (secs)",
    # hue=specs[["language", "name"]].apply(tuple, axis=1),
    hue="dt",
    markers=True,
    style="language",
    palette=sns.color_palette("muted"),
    # palette={
    #     "MATLAB_SAV": sns.color_palette()[0],
    #     "Julia": sns.color_palette()[1],
    # },
    # units="dt",
    # estimator=None,
)
plt.title("Constant total time with varying timesteps, dt, tol, Nx (eps = 0.015)")
plt.yscale("log")
plt.tight_layout()
plt.legend(fontsize="x-small")
plt.savefig(
    "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_normal_IC/plot_MG_SAV/MG_SAV_log_Nx_vs_time_constant_T_printing.pdf"
)
