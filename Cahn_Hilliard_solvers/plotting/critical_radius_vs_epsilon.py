### Plotting for each individual alpha and epsilon. This will plot the radius over time to figure out what the critical radius for each epsilon is. 
# After this, use the plot_inflection_pt.m code to find the inflection points, and record these in critical_radii_epsilon.csv.
#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
indir ="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output"
# #%%
# tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_0.015009.txt",header = 0, index_col=None)
# print(tmp.shape)
# # %%
# tmp['R0'].unique()

# #%%
# sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Epsilon = 0.015009")
# plt.savefig(f"{indir}/critical_radius_vs_epsilon_0.015009.pdf")

# plt.show()
# # %%
# tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_0.030019.txt",header = 0, index_col=None)
# print(tmp.shape)
# sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Epsilon = 0.030019")
# plt.savefig(f"{indir}/critical_radius_vs_epsilon_0.030019.pdf")
# # plt.show()
# # %%
# tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_0.060037.txt",header = 0, index_col=None)
# print(tmp.shape)
# #%%
# sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title("Epsilon = 0.060037")
# plt.savefig(f"{indir}/critical_radius_vs_epsilon_0.060037.pdf")
# # plt.show()
#%%
## Reuse this one
alpha = "0.0"

for epsilon in ["0.011257"]:#"0.0037523","0.0056285", "0.0075047"]:#,"0.060037"]:#,"0.04","0.075047","0.090056"
    folder=f"critical_radius"
    tmp = pd.read_csv(f"{indir}/{folder}/alpha_{alpha}/radius_0.5_level_set_epsilon_{epsilon}_alpha_{alpha}_nx_256.txt",header = 0, index_col=None, sep =",",
                    on_bad_lines='skip')
    # tmp = tmp.drop(tmp[tmp["R0"].isin([0.1,0.11, 0.1111, 0.1112, 0.1113, 0.1114, 0.1115])].index)
    print(tmp.shape)
    tmp = tmp.sort_values("R0")
    sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 'small')
    plt.title(f"Epsilon = {epsilon}, alpha = {alpha}")
    plt.tight_layout()
    plt.savefig(f"{indir}/{folder}/alpha_{alpha}/critical_radius_vs_epsilon_{epsilon}_alpha_{alpha}_nx_256.pdf")
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
#############################################
# Troubleshooting the problems with e= 0.045
#############################################
alpha = "-0.5"

epsilon = "0.045028"
folder=f"critical_radius/"
tmp = pd.read_csv(f"{indir}/{folder}/radius_0.5_level_set_epsilon_{epsilon}_alpha_{alpha}.txt",header = 0, index_col=None, sep =",",
                on_bad_lines='warn')
print(tmp.shape)
## USE WARN TO DELETE THOSE ROWS

#%%
g = sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
g.set_xticks([tmp["time"].min(), tmp["time"].max()])
g.set_yticks([tmp["radius"].min(), tmp["radius"].max()])

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 'small')
plt.title(f"Epsilon = {epsilon}, alpha = {alpha}")
plt.tight_layout()
plt.savefig(f"{indir}/{folder}/critical_radius_vs_epsilon_{epsilon}_alpha_{alpha}.pdf")
plt.close()
# %%
for i,r in tmp.iterrows():
    try:
        f = float(r["time"])
    except:
        print(i)
        print(r)

##### OUTPUT FOR TIME
# 34654 X
# radius       NaN
# time      5..061
# R0           NaN
# Name: 34654, dtype: object
# 35695 X
# radius          NaN
# time      6.717.061
# R0              NaN
# Name: 35695, dtype: object

###### OUTPUT FOR RADIUS
# 34134 X
# radius     NaaN
# time       0.62
# R0        0.061
# Name: 34134, dtype: object
# 34907 X
# radius    Na3725
# time        0.06
# R0           NaN
# Name: 34907, dtype: object
# 38034 X
# radius    N5.67
# time      0.061
# R0          NaN
# Name: 38034, dtype: object

###### OUTPUT FOR R0
# 36737 X
# radius        NaN
# time       8.0625
# R0        0.0.061
# Name: 36737, dtype: object

# %%
