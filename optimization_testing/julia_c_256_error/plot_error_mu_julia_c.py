#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

c_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/C/Test_256"
julia_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/Test_256"

# %%
#c_new in cahn
suffix = "tan_IC"
version = "v4"

def plot_c_new(julia_dir, c_dir, version, suffix, out_folder):
    arr_j = np.genfromtxt(f"{julia_dir}/c_new_{version}__{suffix}.csv", delimiter= " ")
    arr_j3d = arr_j.reshape(-1,256,256)

    arr_c = np.genfromtxt(f"{c_dir}/c_new_{suffix}.csv", delimiter= " ")
    arr_c3d = arr_c.reshape(-1,256,256)

    for i in range(8):
        f, (ax1, ax2, ax3) = plt.subplots(1, 3,  dpi = 200, figsize = (18,5))
        diff = arr_j3d[i,:,:] - arr_c3d[i,:,:]
        plt.suptitle(f"VCycle Run #{i+1}")

        sns.heatmap(arr_c3d[i,:,:].transpose(), ax = ax1, xticklabels = 16, yticklabels = 16,
                    center=(None if np.max(diff) == 0 else 0),
                    cmap='RdBu')
        ax1.set_title("$c-new_{c}$")
        sns.heatmap(diff.transpose(), ax = ax2,xticklabels = 16, yticklabels = 16,
                    center=(None if np.max(diff) == 0 else 0),
                    cmap='RdBu')
        ax2.set_title("$c-new_{julia} - c-new_{c}$")
        ax3.hist(diff)
        ax3.set_title("Distribution of errors (Julia - C)")
        # plt.tight_layout()
        outdir = f"/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia_c_256_error/{out_folder}"
        if os.path.exists(outdir):
            pass
        else: os.mkdir(outdir)
        plt.savefig(f"{outdir}/{i}_c_new_{version}_{suffix}.png")
        plt.show()
        plt.close()

j = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia_c_256_error/output_Julia_on_mac_new_MinJhe"
c = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia_c_256_error/output_C_on_mac_MinJhe"
plot_c_new(j, c, version, suffix, out_folder = "MinJhe_Mac")

# %%
#mu in cahn
version = "v4"
arr_j = np.genfromtxt(f"{julia_dir}/mu_{version}_{suffix}.csv", delimiter= " ")
arr_j3d = arr_j.reshape(-1,256,256)

arr_c = np.genfromtxt(f"{c_dir}/mu_{suffix}.csv", delimiter= " ")
arr_c3d = arr_c.reshape(-1,256,256)

for i in range(8):
    f, (ax1, ax2, ax3) = plt.subplots(1, 3,  dpi = 200, figsize = (18,5))
    diff = arr_j3d[i,:,:] - arr_c3d[i,:,:]
    plt.suptitle(f"VCycle Run #{i+1}")

    sns.heatmap(arr_c3d[i,:,:].transpose(), ax = ax1, xticklabels = 16, yticklabels = 16,
                center=(None if np.max(diff) == 0 else 0),
                cmap='RdBu')
    ax1.set_title("$mu_{c}$")
    sns.heatmap(diff.transpose(), ax = ax2,xticklabels = 16, yticklabels = 16,
                center=(None if np.max(diff) == 0 else 0),
                cmap='RdBu')
    ax2.set_title("$mu_{julia} - mu_{c}$")
    ax3.hist(diff)
    ax3.set_title("Distribution of errors (Julia - C)")
    # plt.tight_layout()
    plt.savefig(f"/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia_c_256_error/tan_IC/{i}_mu_{version}_{suffix}.png")
    plt.show()
    plt.close()
    
# %%
#mu and c_new in vcycle (uf new and wf new)
version = 'v4'
arr_j = np.genfromtxt(f"{julia_dir}/uf_new_5_{version}.csv",
                      delimiter=" ",
                      max_rows=256,
                      skip_header=(4 + 8 + 16 + 32 + 64 + 128))
# arr_j = np.genfromtxt(f"{julia_dir}/uf_new_6_{version}.csv", delimiter= " ", max_rows = 256, skip_header=(2+4+8+16+32+64+128))
# arr_j = np.genfromtxt(f"{julia_dir}/uf_new_1_{version}.csv",delimiter=" ",max_rows=256) for 1 2 and 3
# arr_j = np.genfromtxt(f"{julia_dir}/uf_new_4_{version}.csv", delimiter= " ", max_rows = 256, skip_header=(4+8+16+32+64+128)) for 4 and 5

# arr_j3d = arr_j.reshape(-1,256,256)

arr_c = np.genfromtxt(f"{c_dir}/uf_new_5.csv",
                      delimiter=" ",
                      max_rows=256,
                      skip_header=(4 + 8 + 16 + 32 + 64 + 128))

# arr_c3d = arr_c.reshape(-1,256,256)
# for i in range(8):
f, (ax1, ax2, ax3) = plt.subplots(1, 3,  dpi = 200, figsize = (18,5))
diff = arr_j[:,:] - arr_c[:,:]
plt.suptitle(f"VCycle Run #1 after ilevel loop")

sns.heatmap(arr_c[:,:].transpose(), ax = ax1, xticklabels = 16, yticklabels = 16,center=(None if np.max(arr_c) == 0 else 0),
                cmap='RdBu')
ax1.set_title("uf_new_5 in C")
sns.heatmap(diff.transpose(), ax = ax2,xticklabels = 16, yticklabels = 16,center=(None if np.max(diff) == 0 else 0),
                cmap='RdBu')
ax2.set_title("uf_new_5$_{julia}$ - uf_new_5$_{c}$")
ax3.hist(diff)
ax3.set_title("Distribution of errors (Julia - C)")
# plt.tight_layout()
plt.savefig(f"./julia_c_256_error/uf_new_5.png")
plt.show()
plt.close()

# %%
arr_j = np.genfromtxt(f"{julia_dir}/d2f.txt", delimiter= " ", max_rows = 256)

arr_c = np.genfromtxt(f"{c_dir}/d2f.csv", delimiter= " ",max_rows = 256)

f, (ax1, ax2, ax3) = plt.subplots(1, 3,  dpi = 200, figsize = (18,5))
diff = arr_j - arr_c
print(np.max(diff))
print(np.min(diff))
plt.suptitle(f"Relax cycle")

sns.heatmap(arr_c.transpose(), ax = ax1, center=(None if np.max(arr_c) == 0 else 0),
                cmap='RdBu')
ax1.set_title("d2f in C")
sns.heatmap(diff.transpose(), ax = ax2, center=(None if np.max(diff) == 0 else 0),
                cmap='RdBu')
ax2.set_title("d2f$_{julia}$ - d2f$_{c}$")
ax3.hist(diff)
plt.xticks(rotation=90)
ax3.set_title("Distribution of errors (Julia - C)")
# plt.tight_layout()
plt.savefig(f"./julia_c_256_error/d2f.png")
plt.show()
plt.close()

#%%%
arr_c3d = arr_c[:,1].reshape(2, 256, 256)
arr_j3d = arr_j[:,1].reshape(2, 256, 256)
f, axes = plt.subplots(2, 3,  dpi = 200, figsize = (18,12))
for i in range(2):
    diff = arr_j3d[i,:,:] - arr_c3d[i,:,:]
    plt.suptitle(f"Relaxation cycle #{i}")

    sns.heatmap(arr_c3d[i,:,:].transpose(), ax = axes[i,0], xticklabels = 16, yticklabels = 16)
    axes[i,0].set_title("$f[1]_{c}$")
    sns.heatmap(diff.transpose(), ax = axes[i,1],xticklabels = 16, yticklabels = 16)
    axes[i,1].set_title("$f[1]_{julia} - f[1]_{c}$")
    axes[i,2].hist(diff)
    axes[i,2].set_title("Distribution of errors (Julia - C)")
    # plt.tight_layout()
plt.savefig(f"./julia_c_256_error/{i}_f[1]_v5.png")
plt.show()
plt.close()


# %%
#better way to plot than code above
# import matplotlib

# print(matplotlib.__version__)
arr_c = np.genfromtxt(f"{c_dir}/f_full_initial.csv", delimiter= " ", max_rows = 131072)
arr_j = np.genfromtxt(f"{julia_dir}/f_full_initial.txt",
                      delimiter=" ",
                      max_rows=131072)

arr_c3d = arr_c[:,1].reshape(2, 256, 256)
arr_j3d = arr_j[:,1].reshape(2, 256, 256)
fig = plt.figure(constrained_layout=True,dpi = 300, figsize = (18,12))
fig.suptitle('f[1] in relax function', fontsize = 22,fontweight = 'bold')
# create 2x1 subfigs
subfigs = fig.subfigures(nrows=2, ncols=1, hspace=0.08)
for row, subfig in enumerate(subfigs):
    diff = arr_j3d[row, :, :] - arr_c3d[row, :, :]

    subfig.suptitle(f'Relaxation cycle {row}', fontsize=15, fontweight='bold')

    # create 1x3 subplots per subfig
    axs = subfig.subplots(nrows=1, ncols=3)
    for col, ax in enumerate(axs):
        print(col)
        if col == 0:
            sns.heatmap(
                arr_c3d[row, :, :].transpose(),
                ax=ax,
                xticklabels=16,
                yticklabels=16,
                center=(None if np.max(arr_c3d[row, :, :]) == 0 else 0),
                cmap='RdBu')
            ax.set_title("$f[1]_{c}$")
        elif col == 1:
            sns.heatmap(
                diff.transpose(),
                ax=ax,
                xticklabels=16,
                yticklabels=16,
                center = (None if np.max(diff) == 0 else 0), #necessary to fix a bug where all 0 matrix is red not white
                cmap='RdBu')
            ax.set_title("$f[1]_{julia} - f[1]_{c}$")
        elif col == 2:
            ax.hist(diff)

            ax.set_title("Distribution of errors (Julia - C)")
plt.savefig(f"./julia_c_256_error/relax/f[1]_initial.png")
plt.show()

plt.close()

# %%
name = "c_new_before_update"
arr_c = np.genfromtxt(f"{c_dir}/{name}.txt", delimiter= " ", max_rows = 512)

arr_j = np.genfromtxt(f"{julia_dir}/{name}.txt",
                      delimiter=" ",
                      max_rows=512)

arr_c3d = arr_c.reshape(2, 256, 256)
arr_j3d = arr_j.reshape(2, 256, 256)
# print([f"{i:16.15f}" for i in arr_j3d[1, 125:130,0]])
# print([f"{i:16.15f}" for i in arr_c3d[1, 125:130,0]])

fig = plt.figure(constrained_layout=True,dpi = 300, figsize = (18,12))
fig.suptitle(f'{name}', fontsize = 22,fontweight = 'bold')
# create 2x1 subfigs
subfigs = fig.subfigures(nrows=2, ncols=1, hspace=0.08)
for row, subfig in enumerate(subfigs):
    diff = arr_j3d[row, :, :] - arr_c3d[row, :, :]

    subfig.suptitle(f'Relaxation cycle {row}', fontsize=15, fontweight='bold')

    # create 1x3 subplots per subfig
    axs = subfig.subplots(nrows=1, ncols=3)
    for col, ax in enumerate(axs):
        print(col)
        if col == 0:
            sns.heatmap(
                arr_c3d[row, :, :].transpose(),
                ax=ax,
                xticklabels=16,
                yticklabels=16,
                center=(None if np.max(arr_c3d[row, :, :]) == 0 else 0),
                cmap='RdBu')
            ax.set_title(f"{name} in C")
        elif col == 1:
            sns.heatmap(
                diff.transpose(),
                ax=ax,
                xticklabels=16,
                yticklabels=16,
                center = (None if np.max(diff) == 0 else 0), #necessary to fix a bug where all 0 matrix is red not white
                cmap='RdBu')
            ax.set_title(f"{name} (julia) - {name} (c)")
        elif col == 2:
            sns.heatmap(
                arr_j3d[row, :, :].transpose(),
                ax=ax,
                xticklabels=16,
                yticklabels=16,
                center=(None if np.max(arr_j3d[row, :, :]) == 0 else 0),
                cmap='RdBu')
            ax.set_title(f"{name} in Julia")

# plt.savefig(f"./julia_c_256_error/relax/{name}.png")
plt.show()

plt.close()

# %%
# troubleshooting d2f function


def print_problem(name, c = 'txt'):
    arr_c = np.genfromtxt(f"{c_dir}/{name}.{c}", delimiter=" ", max_rows=512)

    arr_j = np.genfromtxt(f"{julia_dir}/{name}.txt",
                          delimiter=" ",
                          max_rows=512)

    arr_c3d = arr_c.reshape(2, 256, 256)
    arr_j3d = arr_j.reshape(2, 256, 256)
    print([f"{i:19.18f}" for i in arr_j3d[1, 124:144, 2]])
    print([f"{i:16.15f}" for i in arr_c3d[1, 124:144, 2]])


print_problem("c_new_before_update")
print_problem("d2f", c='csv')
# %%
#plot residuals
version = "v4"
jul_err = pd.read_csv(f"{julia_dir}/error2/res2_{version}_{suffix}.csv", header = None)
c_err = pd.read_csv(f"{c_dir}/res2_{suffix}.csv", header = None, sep = ' ')
errors = jul_err.iloc[0:8]
errors.columns = ['err']
errors['Code'] = "Julia"
c_err = c_err.iloc[0:8].drop([1,2], axis = 1)
print(c_err)
c_err.columns = ['err']
c_err['Code'] = "C"
errors = pd.concat([errors, c_err])
print(errors)
sns.lineplot(x=errors.index, y=errors['err'], hue = errors['Code'])
plt.yscale('log')
plt.title("Residual for each vcycle iteration (C vs Julia)")
plt.axhline(y=1e-6, linestyle = '--', c='lightgray')
plt.xlabel("VCycle Iteration")
plt.ylabel("Residual")
plt.savefig(f"/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia_c_256_error/tan_IC/residual_plot_{suffix}.png")
# %%
