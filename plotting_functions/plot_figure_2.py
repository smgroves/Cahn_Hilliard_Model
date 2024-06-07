#%%
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/figure_2_spinodal_decomp"

#%%
name = "phi_32_10000_1.0e-6_"

arr_j = np.genfromtxt(f"{indir}/{name}.txt",
                      delimiter=" ")

arr_j3d = arr_j.reshape(-1, 32, 32)
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
plt.show()
# np.savetxt(f"{indir}/phi_32_initial.txt",arr_j3d[t,:,:])

#%%
name = "uf_new_2_v5_"

arr_j = np.genfromtxt(f"{indir}/{name}.csv",
                      delimiter=" ", max_rows = 32)
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
#%%
name = "duc"

arr_j = np.genfromtxt(f"{indir}/{name}.csv",
                      delimiter=" ", max_rows = 16)
idx = [i + 1 for i in range(arr_j.shape[1])]
ticklabels = [i if np.mod(i,10)==0 else None for i in idx ]
sns.heatmap(pd.DataFrame(arr_j, index=idx, columns=idx),
            cmap="RdBu_r",
            xticklabels=ticklabels,
            yticklabels=ticklabels)
plt.savefig(f"{indir}/{name}_restricted.png")
# %%
#%%
# coarse grid correction
name = "coarse_grid_correction_u"
name_inter = "coarse_grid_correction_interpolated_u"
skip_rows = 0
for level in range(4):
    print(level)
    ex = level + 1
    nx = 2**ex
    arr_cgc = np.genfromtxt(f"{indir}/{name}.csv",
                            delimiter=" ",
                            max_rows=nx,
                            skip_header=skip_rows)
    arr_cgc_inter = np.genfromtxt(f"{indir}/{name_inter}.csv",
                                  delimiter=" ",
                                  max_rows=2 * nx,
                                  skip_header=2 * skip_rows)

    f, (ax1, ax2) = plt.subplots(1, 2, dpi=200, figsize=(12, 4.5))
    plt.suptitle(f"Coarse grid correction")
    sns.heatmap(pd.DataFrame(arr_cgc),
                cmap="RdBu_r",
                ax=ax1,
                linewidth=.2,
                linecolor='k')
    ax1.set_title("Before interpolation")
    sns.heatmap(pd.DataFrame(arr_cgc_inter),
                cmap="RdBu_r",
                ax=ax2,
                linewidth=.2,
                linecolor='k')
    ax2.set_title("After interpolation")
    plt.title("Coarse grid correction")
    plt.savefig(f"{indir}/cgc_{level}.png")
    plt.show()
    plt.close()
    skip_rows += nx

# %%
