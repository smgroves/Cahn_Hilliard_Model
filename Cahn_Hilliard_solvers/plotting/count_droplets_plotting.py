import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/"
outdir = f"{indir}/radii_over_time_level_set_plots/domain_0_2_e_0.0075_noisy_cohesin/"
df = pd.read_csv(f"{indir}/simulated_droplet_distributions_e_0.0075_noisy_cohesin.csv", header=0)
print(df.head())
long_df = pd.DataFrame(columns=["cpc",'cohesin','Droplet Count', 'Frequency'])

for i, r in df.iterrows():
    if str(r['norm_counts']) == 'zeros(1;0)': continue
    count_list = str(r['norm_counts'])[1:-2].split(";")
    count_list = [float(d) for d in count_list]
    for j, c in enumerate(count_list):
        if j > 3:
            long_df = pd.concat([
                long_df,
                pd.DataFrame({
                    "cpc": [r["CPC(um)"]],
                    "cohesin": [r["cohesin(um)"]],
                    "Droplet Count": [5],
                    "Frequency": [np.sum(count_list[5:])]
                })
            ],
                                ignore_index=True)

            break
        else:
            long_df = pd.concat([
                long_df,
                pd.DataFrame({
                    "cpc": [r["CPC(um)"]],
                    "cohesin": [r["cohesin(um)"]],
                    "Droplet Count": [j + 1],
                    "Frequency": [c]
                })
            ],
                                ignore_index=True)


sns.swarmplot(data=long_df, x='Droplet Count', y='Frequency', hue = 'cohesin', palette = sns.color_palette("muted"), s = 5)
plt.title("Frequency of droplet counts (within 1.6um of IC) in simulations")
plt.xticks([0,1,2,3,4,],['1','2','3','4',"5+"])
plt.show()
# plt.savefig(f"{outdir}/droplet_counts_sim.png")
