import numpy as np
import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as m
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/ave_mass/"

grid_size = 128
num_it = 300
file_match = f"ave_mass_psi_{grid_size}_{num_it}_*.txt"
path = (os.path.join(indir, file_match))
files = glob.glob(os.path.join(indir, file_match))
print(files)


lines = pd.DataFrame(columns=['ave_mass','tol'])

for f in files:
    # if f.endswith('0.1.txt'): continue
    print(f.split("/")[-1])
    tmp = pd.read_csv(f, header=None, index_col=None)
    tmp.columns = ["ave_mass"]
    tmp["tol"] = float(f.split("_")[-1].split(".txt")[0])
    print(tmp.shape)
    lines = pd.concat([lines, tmp])
    lines['Timestep'] = lines.index
    lines = lines.reindex()


sns.lineplot(data=lines,
             x="Timestep",
             y='ave_mass',
             hue="tol",
             palette='viridis',
             hue_norm=m.colors.LogNorm())
plt.title(f"Average Mass ($M(psi)/(h^2*nx*ny)$) for Grid Size {grid_size}")
plt.ylabel("Average Mass")
plt.savefig(f"{indir}/mass_psi_{grid_size}_{num_it}.png")
