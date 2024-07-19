#%%
# Plot relationship between critical radius and epsilon values for a given alpha. 
# This gives a line plot that allows us to approximate the epsilon value we should use to get a critical radius equal to the CPC critical radius from images.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%%
########################################
# Trend of R0 vs equilibrium radii
########################################

# Plot the final equilibrium radius for a given alpha and epsilon for various R0. We want to see if the relationship
# between critical eq radius and R0 is linear or not, to see if the cases that grow from R0 to critical eq radius
# are converging to some critical eq radius as the R0 approaches the critical initial radius.
indir ="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output"
alpha = "-0.2"
epsilon = "0.060037"
df = pd.read_csv(f"{indir}/critical_radius/alpha_{alpha}/radius_0.5_level_set_epsilon_{epsilon}_alpha_{alpha}.txt", header = 0, index_col=None)
print(df.head())
# wide_df = df.pivot(index='R0', columns='time', values='radius')
wide_df = df.pivot_table(values='radius',index='R0',columns='time')

x = []
y = []
for i,r in wide_df.iterrows():
    if r.hasnans: pass
    else:
        end_idx = r.last_valid_index()
        y.append(r[end_idx])
        x.append(i)

xs = np.array(x)
ys = np.array(y)
coef = np.polyfit(xs,ys,1)
poly1d_fn = np.poly1d(coef) 
# poly1d_fn is now a function which takes in x and returns an estimate for y
f, ax = plt.subplots()
plt.plot(xs,ys, 'o', xs, poly1d_fn(x), '--')
plt.text(0.95,0.1, f"y = {round(coef[0],3)}x+ {round(coef[1],3)}",
         horizontalalignment='right',
      verticalalignment='center',
      transform = ax.transAxes)
plt.xlabel("R0")
plt.ylabel("Final (Equilibrium) radius at T=10")
plt.title(f"Epsilon: {epsilon}, alpha: {alpha}")
plt.savefig(f"{indir}/critical_radius/alpha_{alpha}/final_radius_vs_R0_eps_{epsilon}.pdf")
plt.show()
print(y[0])
# %%
# import scipy.optimize
# xs = np.array(x)
# ys = np.array(y)
# def monoExp(x, m, t, b):
#     return m * np.exp(t * x) + b
# # perform the fit
# p0 = (1, 0.5, -1) # start with values near those we expect
# params, cv = scipy.optimize.curve_fit(monoExp, xs, ys, p0)
# m, t, b = params
# sampleRate = 20_000 # Hz
# tauSec = (1 / t) / sampleRate

# # determine quality of the fit
# squaredDiffs = np.square(ys - monoExp(xs, m, t, b))
# squaredDiffsFromMean = np.square(ys - np.mean(ys))
# rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
# print(f"R² = {rSquared}")

# # plot the results
# plt.plot(xs, ys, '.', label="data")
# plt.plot(xs, monoExp(xs, m, t, b), '--', label="fitted")
# plt.title("Fitted Exponential Curve")


# # inspect the parameters
# print(f"Y = {m} * e^(-{t} * x) + {b}")
# print(f"Tau = {tauSec * 1e6} µs")


#%%
########################################
# all alpha
########################################
df = pd.read_csv("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon.csv", header = 0, index_col=None)

print(df.head())

sns.lineplot(data = df, y = "critical equilibrium radius", x = "epsilon", hue = 'alpha', linestyle = "-",markers='o')
plt.show()

#%%
########################################
# only alpha = -0.5, with max and min
########################################
df = pd.read_csv("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon_-0.5.csv", header = 0, index_col=None)

print(df.head())
y = "critical equilibrium radius (min)"
sns.lineplot(data = df, y = y, x = "epsilon",  marker='o')
plt.title(f"Critical radius vs. epsilon \n {y}, alpha = -0.5")
plt.savefig(f"{y}_vs_epsilon.png")
plt.show()

#%%

plt.figure()
plt.plot(df['epsilon'], df['critical equilibrium radius (min)'], label = "Min (inflection point)", marker = "o")
plt.plot(df['epsilon'], df['critical equilibrium radius (max)'], label = "Max (equilibrium radius)", marker = "o")
plt.legend(loc = "upper right")
plt.title("Critical radius vs. epsilon, alpha = -0.5")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.show()
xs = np.array(df['epsilon'].values[1:4])
ys = np.array(df['critical equilibrium radius (min)'].values[1:4])
coef = np.polyfit(xs,ys,1)
poly1d_fn = np.poly1d(coef) 

xs_max = np.array(df['epsilon'].values)
ys_max = np.array(df['critical equilibrium radius (max)'].values)
coef_max = np.polyfit(xs_max,ys_max,1)
poly1d_fn_max = np.poly1d(coef_max) 

# poly1d_fn is now a function which takes in x and returns an estimate for y
f, ax = plt.subplots()
plt.plot(xs,ys, 'o',label = "Minimum (inflection points)")
plt.plot( [0.030019, 0.09, 0.15], poly1d_fn([0.030019, 0.09, 0.15]), '-', c = '#1f77b4')
plt.plot(xs_max,ys_max, 'o', label = "Maximum ($R_{eq}$ minimum)")
plt.plot( [0.015009, 0.09, 0.15], poly1d_fn_max([0.015009, 0.09, 0.15]), '-', c = '#ff7f0e')
plt.text(0.99,0.25, f"y = {round(coef[0],3)}x+ {round(coef[1],3)}",
         horizontalalignment='right',
      verticalalignment='center',
      transform = ax.transAxes, c = '#1f77b4')
plt.text(.99,0.3, f"y = {round(coef_max[0],3)}x+ {round(coef_max[1],3)}",
         horizontalalignment='right',
      verticalalignment='center',
      transform = ax.transAxes, c = '#ff7f0e')
plt.axhline(y = 0.108, label = "Experimental CPC $R_{critical}$", c = 'k', linestyle = "--")
plt.legend(loc ="lower right")

plt.title("Critical radius vs. epsilon, alpha = -0.5")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.savefig("Critical equilibrium radius (min and max)_vs_epsilon.png")

#%%
########################################
# all alpha, with max and min
########################################
df = pd.read_csv("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon copy.csv", header = 0, index_col=None)
print(df.head())
y = "critical equilibrium radius (max)"
sns.lineplot(data = df, y = y, x = "epsilon",  marker='o', hue = 'alpha')
plt.title(f"Critical radius vs. epsilon \n {y}, alpha = -0.5")
# plt.savefig(f"{y}_vs_epsilon.png")
plt.show()

#%%
df = pd.read_csv("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon copy.csv", header = 0, index_col=None)

alpha = 0
tmp = df.loc[df['alpha']==alpha]
# plt.figure()
# plt.plot(tmp['epsilon'], tmp['critical equilibrium radius (min)'], label = "Min (inflection point)", marker = "o")
# plt.plot(tmp['epsilon'], tmp['critical equilibrium radius (max)'], label = "Max (equilibrium radius)", marker = "o")
# plt.legend(loc = "upper right")
# plt.title(f"Critical radius vs. epsilon, alpha = {alpha}")
# plt.xlabel("Epsilon")
# plt.ylabel("Critical Radius")
# plt.show()
xs = np.array(tmp['epsilon'])#[0:2])
ys = np.array(tmp['critical equilibrium radius (min)'])#[0:2])
coef = np.polyfit(xs,ys,1)
poly1d_fn = np.poly1d(coef) 

# xs_max = np.array(tmp['epsilon'].values)
# ys_max = np.array(tmp['critical equilibrium radius (max)'].values)
# coef_max = np.polyfit(xs_max,ys_max,1)
# poly1d_fn_max = np.poly1d(coef_max) 

# poly1d_fn is now a function which takes in x and returns an estimate for y
f, ax = plt.subplots()
plt.plot(xs,ys, 'o')#,label = "Minimum (inflection points)")
plt.plot( [0, 0.09], poly1d_fn([0, 0.09]), '-', c = '#1f77b4')
# plt.plot(xs_max,ys_max, 'o', label = "Maximum ($R_{eq}$ minimum)")
# plt.plot( [0, 0.09, 0.15], poly1d_fn_max([0, 0.09, 0.15]), '-', c = '#ff7f0e')
plt.text(0.99,0.25, f"y = {round(coef[0],3)}x+ {round(coef[1],3)}",
         horizontalalignment='right',
      verticalalignment='center',
      transform = ax.transAxes, c = '#1f77b4')
# plt.text(.99,0.3, f"y = {round(coef_max[0],3)}x+ {round(coef_max[1],3)}",
#          horizontalalignment='right',
#       verticalalignment='center',
#       transform = ax.transAxes, c = '#ff7f0e')
plt.axhline(y = 0.054, label = "Experimental CPC $R_{critical}$", c = 'k', linestyle = "--")
plt.legend(loc ="lower right")

plt.title(f"Critical radius vs. epsilon, alpha = {alpha}")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.show()
# plt.savefig(f"Critical equilibrium radius (min)_vs_epsilon_alpha_{alpha}.png")
#%%
import scipy.optimize

def monoExp(x, m, t, b):
#     return m * np.exp(t * x) + b
      #b(1).*(b(2).*xdata./(b(3) + xdata) + xdata);
      return m*(t*x/(b+x)+x)## HYPERBOLIC FIT
      # return 
# perform the fit
p0 = (0,0,0) # start with values near those we expect
params, cv = scipy.optimize.curve_fit(monoExp, xs, ys, p0)
print(params)
m, t, b = params
sampleRate = 20_000 # Hz
tauSec = (1 / t) / sampleRate

# determine quality of the fit
squaredDiffs = np.square(ys - monoExp(xs, m, t, b))
squaredDiffsFromMean = np.square(ys - np.mean(ys))
rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
print(f"R² = {rSquared}")

# plot the results
fig = plt.figure(figsize = (5,4))
plt.text(0.75,0.024, f"R² = {round(rSquared,3)}",
         horizontalalignment='right',
      verticalalignment='center',
      transform = ax.transAxes)
plt.plot(xs, ys, '.', label="Data")
plt.plot(np.linspace(0, 0.09), monoExp(np.linspace(0, 0.09), m, t, b), '--', label="Fit")
print(monoExp(0.0125, m, t, b))
plt.title("Hyperbolic-to-Linear Fit of \n Critical Radius vs Epsilon")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.legend()
plt.savefig(f"Critical equilibrium radius_vs_epsilon_alpha_{alpha}_hyperlin.png")
plt.close()
plt.show()

# %%
