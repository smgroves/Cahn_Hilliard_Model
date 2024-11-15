import numpy as np
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as m


def plot(indir, grid_size, num_it, suffix, fileprefix = "ave_mass", name = "Average Mass", subset = True):
    file_match = f"{fileprefix}_{grid_size}_{num_it}_*.txt"
    path = (os.path.join(indir, file_match))
    files = glob.glob(os.path.join(indir, file_match))
    print(files)
    lines = pd.DataFrame(columns=[fileprefix,'tol'])
    for f in files:
        if subset:
            if '0.01' in f: continue
        print(f.split("/")[-1])
        tmp = pd.read_csv(f, header=None, index_col=None)
        tmp.columns = [fileprefix]
        tmp["tol"] = float(f.split(f"{suffix}.txt")[0].split("_")[-1])
        print(tmp.shape)
        lines = pd.concat([lines, tmp])
        lines['Timestep'] = lines.index
        lines = lines.reindex()
    sns.lineplot(data=lines,
                x="Timestep",
                y=fileprefix,
                hue="tol",
                palette='viridis',
                hue_norm=m.colors.LogNorm())
    plt.title(f"{name} for Grid Size {grid_size}")
    plt.tight_layout()
    plt.ylabel(name)
    if subset:
        plt.savefig(f"{indir}/{fileprefix}_{grid_size}_{num_it}_subset.png")
    else:
        plt.savefig(f"{indir}/{fileprefix}_{grid_size}_{num_it}.png")

    plt.close()

# indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/tan_IC_ave_mass_128/"
# grid_size = 128
# num_it = 15000
# suffix = "_tan_IC"
# plot(indir, grid_size, num_it, suffix, name = "Average Mass ($M(phi)/(h^2*nx*ny)$)")
# plot(indir, grid_size, num_it, suffix,"discrete_norm_e", "Discrete Normalized Energy")


indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/tan_IC_ave_mass_256/"
grid_size = 256
num_it = 15000
suffix = "_tan_IC"
plot(indir, grid_size, num_it, suffix, name = "Average Mass ($M(phi)/(h^2*nx*ny)$)", subset = False)
plot(indir, grid_size, num_it, suffix,"discrete_norm_e", "Discrete Normalized Energy",subset= False)
