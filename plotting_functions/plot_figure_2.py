#%%
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/figure_2_spinodal_decomp"

#%%
name = "phi_128_10000_1.0e-6_"

arr_j = np.genfromtxt(f"{indir}/{name}.txt",
                      delimiter=" ")

arr_j3d = arr_j.reshape(-1, 128, 128)
print(arr_j3d.shape)
t = 1
idx = [i + 1 for i in range(arr_j3d.shape[1])]
ticklabels = [i if np.mod(i,10)==0 else None for i in idx ]
sns.heatmap(pd.DataFrame(arr_j3d[t, :, :], index=idx, columns=idx),
            cmap="RdBu_r",
            vmax=1,
            vmin=-1,
            xticklabels=ticklabels,
            yticklabels=ticklabels)
plt.savefig(f"{indir}/{name}_t={t}.png")

#%%
name = "uf_new_2_v5_"

arr_j = np.genfromtxt(f"{indir}/{name}.csv",
                      delimiter=" ", max_rows = 128)
idx = [i + 1 for i in range(arr_j.shape[1])]
ticklabels = [i if np.mod(i,10)==0 else None for i in idx ]
sns.heatmap(pd.DataFrame(arr_j, index=idx, columns=idx),
            cmap="RdBu_r",
            vmax=1,
            vmin=-1,
            xticklabels=ticklabels,
            yticklabels=ticklabels)
plt.savefig(f"{indir}/{name}_post_relaxation.png")
# %%
