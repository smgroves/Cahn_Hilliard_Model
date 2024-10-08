import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def free_energy(x, alpha):
    return (1 / 12) * x * (-4 * alpha * (-3 + x**2) + 3 * x * (-2 + x**2))


def modified_free_energy(x, alpha):
    param_dict = {
        -0.5: [0.679355124356581, -1.2074996663997333, 0.8741663330663997],
        -0.1: [0.7024393586862704, -1.0349986134211147, 0.968331946754448],
        0: [0, -1, 1],
        0.1: [0.7059312080025176, -0.9683319467544478, 1.0349986134211147],
        0.5: [0.679366220486757, -0.8741663330663997, 1.2074996663997333]
    }
    return (1 / 4) * ((param_dict[alpha][1] - x)**2) * (
        (param_dict[alpha][2] - x)**2)
xs = np.linspace(-2, 2, num=400)
full_df = pd.DataFrame(columns = ["x","y","alpha"])
for alpha in [-1, -0.5, -0.1, 0, 0.1, 0.5, 1]:
    ys = [free_energy(x, alpha) for x in xs]
    full_df = pd.concat([full_df, pd.DataFrame({'x': xs, 'y': ys, 'alpha': [alpha]*len(xs)})], ignore_index = True)

ax = sns.lineplot(data=full_df, x="x", y="y", hue='alpha', palette="tab20")
plt.axvline(0, linestyle = "--", color = 'gray')
plt.axhline(0, linestyle = "--", color = 'gray')
plt.title("Free Energy Functional")
plt.xlabel(r"$\phi$")
plt.ylabel(r"Free Energy F($\phi$)")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/manuscript/figure_4/Free_energy.pdf")
plt.close()

xs = np.linspace(-1.5, 1.5, num=400)
full_df = pd.DataFrame(columns = ["x","y","alpha"])
for alpha in [ -0.5, -0.1, 0, 0.1, 0.5]:
    ys = [modified_free_energy(x, alpha) for x in xs]
    full_df = pd.concat([full_df, pd.DataFrame({'x': xs, 'y': ys, 'alpha': [alpha]*len(xs)})], ignore_index = True)

ax = sns.lineplot(data=full_df, x="x", y="y", hue='alpha', palette="coolwarm")
plt.axvline(0, linestyle = "--", color = 'gray')
plt.axhline(0, linestyle = "--", color = 'gray')
plt.title("Modified Free Energy Functional")
plt.xlabel(r"$\phi$")
plt.ylabel(r"Modified Free Energy $\tilde{F}(\phi$)")
sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/manuscript/figure_4/modified_free_energy.pdf")
plt.close()
