#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

c_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/C/Test_256"
julia_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/Test_256"
version = "v4"
arr_j = np.genfromtxt(f"{julia_dir}/c_new_{version}.csv", delimiter= " ")
arr_j3d = arr_j.reshape(-1,256,256)

arr_c = np.genfromtxt(f"{c_dir}/c_new.csv", delimiter= " ")
arr_c3d = arr_c.reshape(-1,256,256)

# %%
for i in range(8):
    f, (ax1, ax2, ax3) = plt.subplots(1, 3,  dpi = 200, figsize = (18,5))
    diff = arr_j3d[i,:,:] - arr_c3d[i,:,:]
    plt.suptitle(f"VCycle Run #{i+1}")

    sns.heatmap(arr_c3d[i,:,:].transpose(), ax = ax1, xticklabels = 16, yticklabels = 16)
    ax1.set_title("$c-new_{c}$")
    sns.heatmap(diff.transpose(), ax = ax2,xticklabels = 16, yticklabels = 16)
    ax2.set_title("$c-new_{julia} - c-new_{c}$")
    ax3.hist(diff)
    ax3.set_title("Distribution of errors (Julia - C)")
    # plt.tight_layout()
    plt.savefig(f"./julia_c_256_error/{i}_c_new_{version}.png")
    plt.show()
    plt.close()

    
# %%
