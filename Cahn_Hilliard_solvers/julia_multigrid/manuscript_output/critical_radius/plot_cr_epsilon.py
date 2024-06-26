#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%%
df = pd.read_csv("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon.csv", header = 0, index_col=None)

print(df.head())

sns.lineplot(data = df, y = "critical equilibrium radius", x = "epsilon", hue = 'alpha', linestyle = "-",markers='o')
plt.show()

#%%
indir ="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output"
alpha = "-0.5"
epsilon = "0.015009"
df = pd.read_csv(f"{indir}/critical_radius/alpha_{alpha}/radius_0.5_level_set_epsilon_{epsilon}_alpha_{alpha}.txt", header = 0, index_col=None)
print(df.head())
wide_df = df.pivot(index='R0', columns='time', values='radius')


# %%
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

