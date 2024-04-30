#%%

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#%%
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/SAV"
mat = pd.read_csv(f"{indir}/MATLAB/D.csv", header = None, index_col= None)
mat = mat.T

jul = pd.read_csv(f"{indir}/Julia/D.csv", header = None, index_col= None)

all = pd.merge(mat, jul, left_index=True, right_index=True, suffixes=['_MATLAB', '_Julia'])
all["timestep"] = range(15001)
all.columns = [ 'MATLAB','Julia','timestep']
print(all.head())

# %%

all_long = pd.melt(all, id_vars= 'timestep', value_vars=["MATLAB",'Julia'], var_name="Code")
sns.lineplot(data=all_long, x = "timestep", y ='value', hue = 'Code')
plt.ylabel('(r-sqrt(E1))/sqrt(E1)')
plt.title("Relative Error of r from E1")
plt.savefig("error.png")

# %%
all['diff'] =all['MATLAB'] - all['Julia']
print(all['diff'].max())
# %%
f, ax = plt.subplots()
sns.lineplot(data = all, x = 'timestep', y = 'diff')
max = all['diff'].max()
min = all['diff'].min()
plt.title("Difference in Relative Error of r from E1")
plt.ylabel("$E_{MATLAB} - E_{Julia}$")
plt.text(.99, .99,  f'Max: {float(f"{max:.4g}"):g} \n Min: {float(f"{min:.4g}"):g}', ha='right', va='top', transform=ax.transAxes)
plt.savefig("matlab_vs_julia_error.png")
# %%

mat = pd.read_csv(f"{indir}/MATLAB/E.csv", header = None, index_col= None)
mat = mat.T
jul = pd.read_csv(f"{indir}/Julia/E.csv", header = None, index_col= None)
all = pd.merge(mat, jul, left_index=True, right_index=True, suffixes=['_MATLAB', '_Julia'])
all["timestep"] = range(15001)
all.columns = [ 'MATLAB','Julia','timestep']
print(all.head())
all_long = pd.melt(all, id_vars= 'timestep', value_vars=["MATLAB",'Julia'], var_name="Code")
sns.lineplot(data=all_long, x = "timestep", y ='value', hue = 'Code')
plt.ylabel('E')
plt.title("Modified Energy")
plt.savefig("energy.png")
all['diff'] =all['MATLAB'] - all['Julia']
print(all['diff'].max())
f, ax = plt.subplots()
sns.lineplot(data = all, x = 'timestep', y = 'diff')
max = all['diff'].max()
min = all['diff'].min()
plt.title("Difference in Modified Energy")
plt.ylabel("$E_{MATLAB} - E_{Julia}$")
plt.text(.99, .99,  f'Max: {float(f"{max:.4g}"):g} \n Min: {float(f"{min:.4g}"):g}', ha='right', va='top', transform=ax.transAxes)
plt.savefig("matlab_vs_julia_energy.png")

#%%

m_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/SAV/MATLAB"
julia_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/SAV/Julia"

arr_m = np.genfromtxt(f"{m_dir}/Phi_MATLAB.csv", delimiter= ",")
arr_j = np.genfromtxt(f"{julia_dir}/Phi_julia.csv",
                      delimiter=",")

arr_m3d = arr_m.reshape(-1, 128, 128)
arr_j3d = arr_j.reshape(-1, 128, 128)
print(arr_j3d.shape)
# %%
def calc_resid(arr1, arr2):
    resids = []
    for i in range(301):
        r = np.sum((arr1[i,:,:]-arr2[i,:,:])**2)
        resids.append(r)
    return resids

print(calc_resid(arr_m3d, arr_j3d))

# %%
rs = calc_resid(arr_m3d, arr_j3d)
f, ax = plt.subplots()
sns.lineplot( x = [q*50 for q in range(301)], y = rs)
max = np.max(rs)
min = np.min(rs)
plt.title("Residual (MATLAB vs Julia)")
plt.xlabel("Timestep")
plt.ylabel("Sum of squared error: $\\sum(\\phi_{MATLAB}-\\phi_{Julia})^2$")
plt.text(.01, .99,  f'Max: {float(f"{max:.4g}"):g} \n Min: {float(f"{min:.4g}"):g}', ha='left', va='top', transform=ax.transAxes)
plt.savefig("matlab_vs_julia_residuals.png")
plt.show()
# %%
