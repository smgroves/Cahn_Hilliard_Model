# %%
# Plot relationship between critical radius and epsilon values for a given alpha.
# This gives a line plot that allows us to approximate the epsilon value we should use to get a critical radius equal to the CPC critical radius from images.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams["font.family"] = "Arial"

# %%
########################################
# Trend of R0 vs equilibrium radii
########################################

# Plot the final equilibrium radius for a given alpha and epsilon for various R0. We want to see if the relationship
# between critical eq radius and R0 is linear or not, to see if the cases that grow from R0 to critical eq radius
# are converging to some critical eq radius as the R0 approaches the critical initial radius.
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output"
alpha = "-0.2"
epsilon = "0.060037"
df = pd.read_csv(
    f"{indir}/critical_radius/alpha_{alpha}/radius_0.5_level_set_epsilon_{epsilon}_alpha_{alpha}.txt",
    header=0,
    index_col=None,
)
print(df.head())
# wide_df = df.pivot(index='R0', columns='time', values='radius')
wide_df = df.pivot_table(values="radius", index="R0", columns="time")

x = []
y = []
for i, r in wide_df.iterrows():
    if r.hasnans:
        pass
    else:
        end_idx = r.last_valid_index()
        y.append(r[end_idx])
        x.append(i)

xs = np.array(x)
ys = np.array(y)
coef = np.polyfit(xs, ys, 1)
poly1d_fn = np.poly1d(coef)
# poly1d_fn is now a function which takes in x and returns an estimate for y
f, ax = plt.subplots()
plt.plot(xs, ys, "o", xs, poly1d_fn(x), "--")
plt.text(
    0.95,
    0.1,
    f"y = {round(coef[0],3)}x+ {round(coef[1],3)}",
    horizontalalignment="right",
    verticalalignment="center",
    transform=ax.transAxes,
)
plt.xlabel("R0")
plt.ylabel("Final (Equilibrium) radius at T=10")
plt.title(f"Epsilon: {epsilon}, alpha: {alpha}")
plt.savefig(
    f"{indir}/critical_radius/alpha_{alpha}/final_radius_vs_R0_eps_{epsilon}.pdf"
)
plt.show()
print(y[0])
# %%
# import scipy.optimize
# xs = np.array(x)
# ys = np.array(y)
# def h2l(x, m, t, b):
#     return m * np.exp(t * x) + b
# # perform the fit
# p0 = (1, 0.5, -1) # start with values near those we expect
# params, cv = scipy.optimize.curve_fit(h2l, xs, ys, p0)
# m, t, b = params
# sampleRate = 20_000 # Hz
# tauSec = (1 / t) / sampleRate

# # determine quality of the fit
# squaredDiffs = np.square(ys - h2l(xs, m, t, b))
# squaredDiffsFromMean = np.square(ys - np.mean(ys))
# rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
# print(f"R² = {rSquared}")

# # plot the results
# plt.plot(xs, ys, '.', label="data")
# plt.plot(xs, h2l(xs, m, t, b), '--', label="fitted")
# plt.title("Fitted Exponential Curve")


# # inspect the parameters
# print(f"Y = {m} * e^(-{t} * x) + {b}")
# print(f"Tau = {tauSec * 1e6} µs")


# %%
########################################
# all alpha
########################################
df = pd.read_csv(
    "/nonlinear_multigrid/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon.csv",
    header=0,
    index_col=None,
)

print(df.head())

sns.lineplot(
    data=df,
    y="critical equilibrium radius",
    x="epsilon",
    hue="alpha",
    linestyle="-",
    markers="o",
)
plt.show()

# %%
########################################
# only alpha = -0.5, with max and min
########################################
df = pd.read_csv(
    "/nonlinear_multigrid/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon_-0.5.csv",
    header=0,
    index_col=None,
)

print(df.head())
y = "critical equilibrium radius (min)"
sns.lineplot(data=df, y=y, x="epsilon", marker="o")
plt.title(f"Critical radius vs. epsilon \n {y}, alpha = -0.5")
plt.savefig(f"{y}_vs_epsilon.png")
plt.show()

# %%

plt.figure()
plt.plot(
    df["epsilon"],
    df["critical equilibrium radius (min)"],
    label="Min (inflection point)",
    marker="o",
)
plt.plot(
    df["epsilon"],
    df["critical equilibrium radius (max)"],
    label="Max (equilibrium radius)",
    marker="o",
)
plt.legend(loc="upper right")
plt.title("Critical radius vs. epsilon, alpha = -0.5")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.show()
xs = np.array(df["epsilon"].values[1:4])
ys = np.array(df["critical equilibrium radius (min)"].values[1:4])
coef = np.polyfit(xs, ys, 1)
poly1d_fn = np.poly1d(coef)

xs_max = np.array(df["epsilon"].values)
ys_max = np.array(df["critical equilibrium radius (max)"].values)
coef_max = np.polyfit(xs_max, ys_max, 1)
poly1d_fn_max = np.poly1d(coef_max)

# poly1d_fn is now a function which takes in x and returns an estimate for y
f, ax = plt.subplots()
plt.plot(xs, ys, "o", label="Minimum (inflection points)")
plt.plot([0.030019, 0.09, 0.15], poly1d_fn([0.030019, 0.09, 0.15]), "-", c="#1f77b4")
plt.plot(xs_max, ys_max, "o", label="Maximum ($R_{eq}$ minimum)")
plt.plot(
    [0.015009, 0.09, 0.15], poly1d_fn_max([0.015009, 0.09, 0.15]), "-", c="#ff7f0e"
)
plt.text(
    0.99,
    0.25,
    f"y = {round(coef[0],3)}x+ {round(coef[1],3)}",
    horizontalalignment="right",
    verticalalignment="center",
    transform=ax.transAxes,
    c="#1f77b4",
)
plt.text(
    0.99,
    0.3,
    f"y = {round(coef_max[0],3)}x+ {round(coef_max[1],3)}",
    horizontalalignment="right",
    verticalalignment="center",
    transform=ax.transAxes,
    c="#ff7f0e",
)
plt.axhline(y=0.108, label="Experimental CPC $R_{critical}$", c="k", linestyle="--")
plt.legend(loc="lower right")

plt.title("Critical radius vs. epsilon, alpha = -0.5")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.savefig("Critical equilibrium radius (min and max)_vs_epsilon.png")

# %%
########################################
# all alpha, REUSE THIS ONE
########################################
df = pd.read_csv(
    "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/critical_radius/critical_radii_epsilon copy.csv",
    header=0,
    index_col=None,
)
print(df.head())
df_0 = df.loc[df["alpha"] == 0]
y = "critical equilibrium radius (min)"
sns.lineplot(
    data=df_0,
    y=y,
    x="epsilon",
    markers=True,
    style="Nx",
    hue="Nx",
    palette="muted",
    alpha=0.6,
)
plt.title(f"Critical radius vs. epsilon")
plt.ylabel("Critical equilibrium radius")
plt.savefig(
    f"/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/critical_radius/cr_vs_epsilon_alpha_0.pdf"
)
plt.close()
plt.show()

# %%
df = pd.read_csv(
    "critical_radii_epsilon copy.csv",
    header=0,
    index_col=None,
)

alpha = 0
tmp = df.loc[df["alpha"] == alpha]
# plt.figure()
# plt.plot(tmp['epsilon'], tmp['critical equilibrium radius (min)'], label = "Min (inflection point)", marker = "o")
# plt.plot(tmp['epsilon'], tmp['critical equilibrium radius (max)'], label = "Max (equilibrium radius)", marker = "o")
# plt.legend(loc = "upper right")
# plt.title(f"Critical radius vs. epsilon, alpha = {alpha}")
# plt.xlabel("Epsilon")
# plt.ylabel("Critical Radius")
# plt.show()
xs = np.array(tmp["epsilon"])  # [0:2])
ys = np.array(tmp["critical equilibrium radius (min)"])  # [0:2])
coef = np.polyfit(xs, ys, 1)
poly1d_fn = np.poly1d(coef)

# xs_max = np.array(tmp['epsilon'].values)
# ys_max = np.array(tmp['critical equilibrium radius (max)'].values)
# coef_max = np.polyfit(xs_max,ys_max,1)
# poly1d_fn_max = np.poly1d(coef_max)

# poly1d_fn is now a function which takes in x and returns an estimate for y
f, ax = plt.subplots()
plt.plot(xs, ys, "o")  # ,label = "Minimum (inflection points)")
plt.plot([0, 0.09], poly1d_fn([0, 0.09]), "-", c="#1f77b4")
# plt.plot(xs_max,ys_max, 'o', label = "Maximum ($R_{eq}$ minimum)")
# plt.plot( [0, 0.09, 0.15], poly1d_fn_max([0, 0.09, 0.15]), '-', c = '#ff7f0e')
plt.text(
    0.99,
    0.25,
    f"y = {round(coef[0],3)}x+ {round(coef[1],3)}",
    horizontalalignment="right",
    verticalalignment="center",
    transform=ax.transAxes,
    c="#1f77b4",
)
# plt.text(.99,0.3, f"y = {round(coef_max[0],3)}x+ {round(coef_max[1],3)}",
#          horizontalalignment='right',
#       verticalalignment='center',
#       transform = ax.transAxes, c = '#ff7f0e')
plt.axhline(y=0.054, label="Experimental CPC $R_{critical}$", c="k", linestyle="--")
plt.legend(loc="lower right")

plt.title(f"Critical radius vs. epsilon, alpha = {alpha}")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.show()
# plt.savefig(f"Critical equilibrium radius (min)_vs_epsilon_alpha_{alpha}.png")
# %% with hyperbolic to linear fit
########################################
# FINAL FIGURE CR vs E
########################################
import scipy.optimize

df = pd.read_csv(
    "critical_radii_epsilon copy.csv",
    header=0,
    index_col=None,
)

alpha = 0
tmp = df.loc[df["alpha"] == alpha]
xs = np.array(tmp["epsilon"])  # [0:2])
ys = np.array(tmp["critical equilibrium radius (min)"])  # [0:2])
# ys = np.array(tmp["critical initial radius"])  # [0:2])

print(df.head())


def h2l(x, m, t, b):
    #     return m * np.exp(t * x) + b
    # b(1).*(b(2).*xdata./(b(3) + xdata) + xdata);
    return m * (t * x / (b + x) + x)  ## HYPERBOLIC FIT
    # return


# hyperbolic to linear
# perform the fit
p0 = (0, 0, 0)  # start with values near those we expect
params, cv = scipy.optimize.curve_fit(h2l, xs, ys, p0)
print(params)
m, t, b = params
# determine quality of the fit
squaredDiffs = np.square(ys - h2l(xs, m, t, b))
squaredDiffsFromMean = np.square(ys - np.mean(ys))
rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
print(f"H2L R² = {rSquared}")

# %%
# plot the results
f, ax = plt.subplots(figsize=(5, 4))

# fig = plt.figure(figsize = (5,4))
plt.text(
    0.95,
    0.024,
    f"R² = {round(rSquared,3)}",
    horizontalalignment="right",
    verticalalignment="center",
    transform=ax.transAxes,
)
# plt.plot(xs, ys, '.', label="Simulation results")
markers = {128: "o", 256: "^"}
plt.plot(
    np.linspace(0, 0.09),
    h2l(np.linspace(0, 0.09), m, t, b),
    "--",
    label="Fit",
    c="gray",
)
# sns.scatterplot(
#     data=tmp,
#     x="epsilon",
#     y="critical equilibrium radius (min)",
#     # hue="Nx",
#     palette=sns.color_palette("bright"),
#     edgecolor="k",
#     markers=markers,
#     style="Nx",
#     c="k",
# )
for cat, group in tmp.groupby("Nx"):
    ax.scatter(
        group["epsilon"],
        group["critical equilibrium radius (min)"],
        edgecolors="k",
        facecolors="none",
        marker=markers[cat],
        label=cat,
    )


print(h2l(0.0125, m, t, b))
plt.title("Hyperbolic-to-Linear Fit of\n Critical Radius vs Epsilon")
plt.xlabel(r"Epsilon ($ \epsilon $)")
plt.ylabel(r"Critical Equilibrium Radius ($R_c$)")
plt.legend()
plt.tight_layout()
plt.savefig(
    f"Critical equilibrium radius_vs_epsilon_alpha_{alpha}_hyperlin_128_256_black.pdf"
)
plt.close()
# plt.show()

# %%


def linFit(x, m, b):
    return m * x + b


def logFit(x, a, b):
    return a * np.log(x) + b


def hyperbolic(x, m, t, b):
    return m * (t * x / (b + x))  ## HYPERBOLIC FIT
    # return


def calculate_bic(xs, ys, func, p0):
    """
    Calculate the best fit parameters and BIC for a given function.

    Parameters:
    xs (array): Independent variable data.
    ys (array): Dependent variable data.
    func (callable): The model function.
    p0 (list): Initial guess for parameters.

    Returns:
    dict: Dictionary containing best fit parameters, BIC, and residuals.
    """
    # Fit the function to the data
    params, cv = scipy.optimize.curve_fit(func, xs, ys, p0=p0)

    # Calculate residuals
    residuals = ys - func(xs, *params)

    # Sum of squared residuals
    ssr = np.sum(residuals**2)

    # Number of observations
    n = len(ys)

    # Number of parameters
    k = len(params)

    # Bayesian Information Criterion
    bic = n * np.log(ssr / n) + k * np.log(n)

    return {"best_fit_parameters": params, "bic": bic, "residuals": residuals}


initial_guess = [1.0, 0.0]
results = calculate_bic(xs, ys, linFit, initial_guess)
print("Linear", results["bic"])
plt.plot(
    np.linspace(0, 0.09),
    linFit(np.linspace(0, 0.09), *results["best_fit_parameters"]),
    "--",
    label=f"Linear Fit, BIC = {round(results['bic'],2)}",
    c=sns.color_palette()[3],
    alpha=0.6,
)
initial_guess = [0.0, 0.0, 0.0]
results = calculate_bic(xs, ys, h2l, initial_guess)
print("H2L", results["bic"])
plt.plot(
    np.linspace(0, 0.09),
    h2l(np.linspace(0, 0.09), *results["best_fit_parameters"]),
    "--",
    label=f"H2L Fit, BIC = {round(results['bic'],2)}",
    c=sns.color_palette()[4],
    alpha=0.6,
)

initial_guess = [1, 0.0, 0.0]
results = calculate_bic(xs, ys, hyperbolic, initial_guess)
print("Hyperbolic", results["bic"])
plt.plot(
    np.linspace(0, 0.09),
    hyperbolic(np.linspace(0, 0.09), *results["best_fit_parameters"]),
    "--",
    label=f"Hyperbolic Fit, BIC = {round(results['bic'],2)}",
    c=sns.color_palette()[5],
    alpha=0.6,
)

initial_guess = [1.0, 0.0]
results = calculate_bic(xs, ys, logFit, initial_guess)
print("Log", results["bic"])
plt.plot(
    np.linspace(0, 0.09),
    logFit(np.linspace(0, 0.09), *results["best_fit_parameters"]),
    "--",
    label=f"Log Fit, BIC = {round(results['bic'],2)}",
    c=sns.color_palette()[6],
    alpha=0.6,
)

markers = {128: "o", 256: "X"}

sns.scatterplot(
    data=tmp,
    x="epsilon",
    y="critical equilibrium radius (min)",
    hue="Nx",
    palette=sns.color_palette("bright"),
    edgecolor="k",
    markers=markers,
    style="Nx",
)
print(h2l(0.0125, m, t, b))
plt.title("Different Fits of\n Critical Radius vs Epsilon")
plt.xlabel("Epsilon")
plt.ylabel("Critical Radius")
plt.legend()
plt.savefig(
    f"Critical equilibrium radius_vs_epsilon_alpha_{alpha}_multi_fits_128_256.png"
)
plt.close()
plt.show()


# %%
tmp = df.loc[df["alpha"] == 0]


def crit_init_r_theory(e, V=1):
    return np.cbrt((np.sqrt(6) / (8 * np.pi)) * V * e)


f, ax = plt.subplots(figsize=(5, 4))

plt.plot(
    np.linspace(0, 0.09),
    crit_init_r_theory(np.linspace(0, 0.09)),
    "--",
    label=f"Theory",
    c="grey",
    alpha=0.6,
)

markers = {128: "o", 256: "^"}

# sns.scatterplot(
#     data=tmp,
#     x="epsilon",
#     y="critical initial radius",
#     # hue="Nx",
#     palette=sns.color_palette("bright"),
#     edgecolor="k",
#     markers=markers,
#     style="Nx",
# )
for cat, group in tmp.groupby("Nx"):
    ax.scatter(
        group["epsilon"],
        group["critical initial radius"],
        edgecolors="k",
        facecolors="none",
        marker=markers[cat],
        label=cat,
    )
plt.legend()
plt.title("Critical Initial Radius vs Epsilon")
plt.xlabel(r"Epsilon ($ \epsilon $)")
plt.ylabel(r"Critical Initial Radius ($R_0$)")
plt.tight_layout()
plt.savefig(f"Critical initial radius_vs_epsilon_w_theory_128_256_black.pdf")
plt.close()
# plt.show()
# %%
