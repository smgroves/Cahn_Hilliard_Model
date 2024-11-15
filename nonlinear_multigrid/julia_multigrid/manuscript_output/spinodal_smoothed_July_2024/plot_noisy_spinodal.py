import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colors import ListedColormap, BoundaryNorm

nx = 128
mean = -0.58
sd = "0.8"
color_decider = "Final"
choose_colormap = "1"
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/noisy_IC"

# indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/smoothed_with_multigrid"
data = np.loadtxt(
    f"{indir}/phi_128_3200_1.0e-6_dt_6.25e-6_eps_0.015009369912862116_mean_{mean}_sd_{sd}.txt"
)
# data = np.loadtxt(
#     f"{indir}/phi_128_3200_1.0e-6_dt_6.25e-6_initial_5.txt"
# )
reshaped_data = data.reshape(-1, nx, nx)
final_values = reshaped_data[-1, :, :]  # Shape: (X, Y)
initial_values =reshaped_data[0, :, :]
if color_decider == "Final":
    color_values = final_values
elif color_decider == "Initial":
    color_values = initial_values
###############
# Color map
################
if choose_colormap == "1":
    norm = Normalize(vmin=np.min(color_values), vmax=np.max(color_values))
    colormap = cm.viridis

################
elif choose_colormap == "2":
    colors = [
        '#0000ff',  # Blue (for values < -0.58)
        '#008000',  # Green (for values between -0.58 and 0)
        '#ffff00',  # Yellow (for values between 0 and 0.58)
        '#ff0000'   # Red (for values > 0.58)
    ]
    colormap = ListedColormap(colors)
    bounds = [np.min(color_values), -0.58, 0, 0.58, np.max(color_values)]
    norm = BoundaryNorm(bounds, colormap.N)


plt.figure(figsize = [32,20], dpi = 300)
for i in range(nx):  # Iterate over X
    for j in range(nx):  # Iterate over Y
        color = colormap(norm(color_values[i, j]))

        plt.plot(reshaped_data[:, i, j],
                 linewidth=.5,
                 alpha=0.5,
                 color=color)
        plt.xlabel('Timestep')
        plt.ylabel('Phi')
plt.axhline(y=mean, linestyle="--", color="grey")
plt.axhline(y=-.58, linestyle = "--", color = "black")
plt.axhline(y=.58, linestyle = "--", color = "black")
plt.axhline(y=0, linestyle = "--", color = "black")

plt.title(
    f"Evolution of spinodal decomposition with IC mean = {mean}, sd = {sd}")
sm = cm.ScalarMappable(cmap=colormap, norm=norm)

sm.set_array([])
plt.colorbar(sm, label=f'{color_decider} Value')

plt.savefig(f"{indir}/mean_{mean}_sd_{sd}_lineplot_{color_decider}_range.png")
# plt.savefig(f"{indir}/smoothed_with_multigrid_t_5_lineplot_{color_decider}_range.png")
