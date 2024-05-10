import numpy as np
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as m


def plot(indir, grid_size, num_it, suffix, fileprefix = "ave_mass", name = "Average Mass"):
    file_match = f"{fileprefix}_{grid_size}_{num_it}_*.txt"
    path = (os.path.join(indir, file_match))
    files = glob.glob(os.path.join(indir, file_match))
    print(files)
    lines = pd.DataFrame(columns=[fileprefix,'tol'])
    for f in files:
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
    plt.title(f"{name} ($M(phi)/(h^2*nx*ny)$) for Grid Size {grid_size}")
    plt.ylabel(name)
    plt.savefig(f"{indir}/{fileprefix}_{grid_size}_{num_it}_subset.png")
    plt.close()

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/tan_IC_ave_mass_128/"
grid_size = 128
num_it = 15000
suffix = "_tan_IC"
plot(indir, grid_size, num_it, suffix)
plot(indir, grid_size, num_it, suffix,"discrete_norm_e", "Discrete Normalized Energy")
