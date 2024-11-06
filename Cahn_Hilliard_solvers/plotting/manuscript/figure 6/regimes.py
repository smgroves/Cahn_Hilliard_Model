import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/manuscript/figure 6"
df = pd.read_csv(f"{indir}/alpha_0.0_CPC_droplet_regimes.csv", header=0, index_col=None)

# df['CPC(unitless)'] = df['CPC(unitless)']*3.2
# df['cohesin(unitless)'] = df['cohesin(unitless)']*3.2

plt.figure(figsize=(6, 4), dpi=300)
ax = sns.scatterplot(data=df, x="CPC(unitless)", y="cohesin(unitless)", hue="dynamics")
plt.xlabel(f"Initial radius for IC condensate ($\mu$m)")
plt.ylabel(r"Cohesin width ($\mu$m)")
plt.title("Regimes of Droplet Dynamics for Initial Condition Parameters")
legend_handles, legend_labels = ax.get_legend_handles_labels()
plt.legend(
    legend_handles,
    legend_labels,
    ncol=4,
    title="Dynamics",
    loc="upper center",
    bbox_to_anchor=(0.5, -0.15),
    fontsize="x-small",
)

# sns.move_legend(ax, "lower center", bbox_to_anchor=(1, 1.05))

plt.tight_layout()
# plt.show()
plt.savefig(f"{indir}/regimes.pdf")
