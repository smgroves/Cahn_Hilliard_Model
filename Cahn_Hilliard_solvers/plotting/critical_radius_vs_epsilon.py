#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
indir ="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output"
#%%
tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_0.015009.txt",header = 0, index_col=None)
print(tmp.shape)
# %%
tmp['R0'].unique()

#%%
sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Epsilon = 0.015009")
plt.savefig(f"{indir}/critical_radius_vs_epsilon_0.015009.pdf")

plt.show()
# %%
tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_0.030019.txt",header = 0, index_col=None)
print(tmp.shape)
sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Epsilon = 0.030019")
plt.savefig(f"{indir}/critical_radius_vs_epsilon_0.030019.pdf")
# plt.show()
# %%
tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_0.060037.txt",header = 0, index_col=None)
print(tmp.shape)
#%%
sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Epsilon = 0.060037")
plt.savefig(f"{indir}/critical_radius_vs_epsilon_0.060037.pdf")
# plt.show()
#%%
## Reuse this one
alpha = "-0.5"

for epsilon in ["0.045028"]:#,"0.060037"]:#,"0.04","0.075047","0.090056"
    folder=f"critical_radius/"
    tmp = pd.read_csv(f"{indir}/{folder}/radius_0.5_level_set_epsilon_{epsilon}_alpha_{alpha}.txt",header = 0, index_col=None, sep =",",
                    on_bad_lines='skip')
    print(tmp.shape)
    sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 'small')
    plt.title(f"Epsilon = {epsilon}, alpha = {alpha}")
    plt.tight_layout()
    plt.savefig(f"{indir}/{folder}/critical_radius_vs_epsilon_{epsilon}_alpha_{alpha}.pdf")
    plt.close()
#%%
epsilon = "0.030019"
radius = 0.0148
tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_{epsilon}.txt",header = 0, index_col=None)
tmp_old = pd.read_csv(f"{indir}/from_Rivanna_old/radius_0.5_level_set_epsilon_{epsilon}.txt",header = 0, index_col=None)
tmp = tmp.loc[tmp['R0']==0.148].dropna()
tmp_old = tmp_old.loc[tmp_old['R0']==0.148].dropna()
plt.plot(tmp['time'].values,tmp['radius'].values, label = 'new IC')
plt.plot(tmp_old['time'].values,tmp_old['radius'].values, label = 'old IC')
plt.legend()
plt.title("R0 = 0.148 for epsilon 0.033019")
plt.xlabel("Time")
plt.ylabel("Radius")
plt.savefig(f"{indir}/r0_{radius}_new_old_IC.pdf")

# %%
#from scipy.ndimage import gaussian_filter1d

# data = tmp.loc[tmp['R0']==0.133].dropna()['radius'].values
# data = data/np.max(data)

# smooth = gaussian_filter1d(data, 100)
# smooth_d1 = (np.gradient(data))

# # compute second derivative
# smooth_d2 = np.diff(np.diff(data))
# print(smooth_d2)
# # find switching points
# infls = np.where(np.diff(np.sign(smooth_d2)))[0]

# # plot results
# plt.plot(data, label='Data')
# plt.show()
# # plot results
# plt.plot(smooth_d1, label='Data')
# plt.ylim(-0.005,0.01)

# plt.show()
# plt.plot(smooth_d2, label='Second Derivative (scaled)')
# for i, infl in enumerate(infls, 1):
#     plt.axvline(x=infl, color='k', label=f'Inflection Point {i}')
# plt.ylim(-0.005,0.01)
# plt.legend(bbox_to_anchor=(1.55, 1.0))
# %%
epsilon = np.array([0.015009, 0.030019, 0.060037])
critical_radius = np.array([0.07,0.095,0.14])
plt.plot(epsilon, critical_radius, marker = "o",label = "Data")
from sklearn.linear_model import LinearRegression
x = epsilon.reshape((-1, 1))
y = critical_radius
model = LinearRegression().fit(x,y)
y_pred = model.intercept_ + model.coef_ * x
plt.plot(epsilon, y_pred, linestyle = "--", label = "Regression")
plt.legend(bbox_to_anchor=(1.55, 1.0))
eps =((0.108-model.intercept_)/model.coef_)[0]
plt.annotate(f"Eps($R_c$= .108)={np.round(eps,6)}", xy = (eps, 0.108), 
             xytext=(1.1*eps, 0.108),
             va="center", ha="left",arrowprops=dict(arrowstyle="-|>"))
plt.title("Epsilon versus critical radius for single droplet")
plt.savefig(f"{indir}/epsilon_vs.critical_radius.pdf")

# %%
