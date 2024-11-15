import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/critical_radius/flat_interface"
nx = 128
alpha = -0.1
offset = -2
tmp = pd.read_csv(f"{indir}/nx_{nx}_alpha_{alpha}_eps_0.015009_equilibrium.csv",
            header=None,
            index_col=None)

#beta, phi_1, phi_2
param_dict = {-0.5:[0.679355124356581,-1.2074996663997333,0.8741663330663997],
              -0.1:[0.7024393586862704,-1.0349986134211147,0.968331946754448],
              0.1:[0.7059312080025176,-0.9683319467544478, 1.0349986134211147],
              0.5:[0.679366220486757,-0.8741663330663997, 1.2074996663997333]}

def equilibrium_interface(x,beta, phi_1, phi_2):
    return (phi_2 + phi_1*np.exp(-x/beta))/(1+np.exp(-x/beta))

h = 1 / nx
M = 8
gam = M * (1 / nx) / (2 * np.sqrt(2) * np.arctanh(0.9))
interface_location=nx/2
xx = np.linspace(0,h*(nx-1),nx)
yy = [equilibrium_interface((x-(interface_location+offset)*h)/gam, param_dict[alpha][0],param_dict[alpha][1],param_dict[alpha][2]) for x in xx]
plt.plot(yy, label = "theory")
plt.plot(tmp.iloc[int(nx / 2)], label="simulation")
residual = np.sum([t - s for t,s in zip(yy, tmp.iloc[int(nx / 2)])])
plt.axhline(y=0, linestyle="--", color='lightgrey')
plt.axvline(x=int(interface_location), linestyle="--", color='lightgrey')
plt.title(
    f"Theory versus simulation of equilibrium interface (last timepoint) \nalpha = {alpha}, offset = {offset}"
)
plt.annotate(f"Residual = {round(residual,3)}",xy=(1, 0.01), xycoords='axes fraction',ha = 'right')
plt.xlabel("Gridpoint")
plt.ylabel(r"$\phi$")
plt.legend()
plt.savefig(f"{indir}/{alpha}_equilibrium_interface.png")
