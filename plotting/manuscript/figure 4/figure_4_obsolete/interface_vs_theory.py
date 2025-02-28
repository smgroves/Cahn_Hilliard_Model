# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

nx = 128
alpha = -0.1
offset = -2

# beta, phi_1, phi_2
param_dict = {
    -0.5: [0.679355124356581, -1.2074996663997333, 0.8741663330663997],
    -0.1: [0.7024393586862704, -1.0349986134211147, 0.968331946754448],
    0: [1 / np.sqrt(2), -1, 1],
    0.1: [0.7059312080025176, -0.9683319467544478, 1.0349986134211147],
    0.5: [0.679366220486757, -0.8741663330663997, 1.2074996663997333],
}


def equilibrium_interface(x, beta, phi_1, phi_2):
    return (phi_2 + phi_1 * np.exp(-x / beta)) / (1 + np.exp(-x / beta))


h = 1 / nx
M = 8
gam = M * (1 / nx) / (2 * np.sqrt(2) * np.arctanh(0.9))
interface_location = nx / 2
xx = np.linspace(0, h * (nx - 1), nx)
yy = [
    equilibrium_interface(
        (x - (interface_location) * h) / gam,
        param_dict[alpha][0],
        param_dict[alpha][1],
        param_dict[alpha][2],
    )
    for x in xx
]
plt.plot(yy, label="theory")
plt.axhline(y=0, linestyle="--", color="lightgrey")
plt.axvline(x=int(interface_location), linestyle="--", color="lightgrey")
plt.title(f"Interface profile")
plt.xlabel("Gridpoint")
plt.ylabel(r"$\phi$")
plt.legend()
plt.show()
# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# beta, phi_1, phi_2
param_dict = {
    -0.5: [0.679355124356581, -1.2074996663997333, 0.8741663330663997],
    -0.1: [0.7024393586862704, -1.0349986134211147, 0.968331946754448],
    0: [1 / np.sqrt(2), -1, 1],
    0.1: [0.7059312080025176, -0.9683319467544478, 1.0349986134211147],
    0.5: [0.679366220486757, -0.8741663330663997, 1.2074996663997333],
}


def initialization_from_function_v2(nx, ny, beta, phi_1, phi_2, R0=0.1, gam=0.01):
    # Step size
    h = 1.0 / nx

    # Create x and y arrays
    x = h * np.arange(nx)
    y = h * np.arange(ny)

    # Create meshgrid
    xx, yy = np.meshgrid(x, y, indexing="ij")

    # Compute the distance from the center (0.5, 0.5)
    R = np.sqrt((xx - 0.5) ** 2 + (yy - 0.5) ** 2)

    # Compute phi based on the given formula
    phi = (phi_2 + phi_1 * np.exp((R - R0) / (beta * gam))) / (
        1 + np.exp((R - R0) / (beta * gam))
    )

    return phi


# Parameters
nx, ny = 256, 256

R0 = 0.1

M = 8
gam = M * (1 / nx) / (2 * np.sqrt(2) * np.arctanh(0.9))

full_df = pd.DataFrame(columns=["x", "y", "alpha"])
# for alpha in [-0.5, -0.1, 0, 0.1, 0.5]:
alpha = 0
beta, phi_1, phi_2 = param_dict[alpha]
# Generate phi
phi = initialization_from_function_v2(nx, ny, beta, phi_1, phi_2, R0, gam)

# Plot 1D slice across the middle (halfway through the matrix)
middle_index = ny // 2
phi_slice = phi[:, middle_index]
xs = np.linspace(0, 1, nx)
full_df = pd.concat(
    [full_df, pd.DataFrame({"x": xs, "y": phi_slice, "alpha": [alpha] * len(xs)})]
)

# Plot
# sns.lineplot(data=full_df, x="x", y="y", hue="alpha", palette="coolwarm")
sns.lineplot(data=full_df, x="x", y="y")
plt.xlabel("Domain")
plt.ylabel(r"$\phi$")
plt.title(r"Interface Profile for R0 = 0.1, $\epsilon$ = 0.015009")
plt.grid(True)
plt.savefig(
    "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting/manuscript/figure_4/interface_alpha_0.pdf"
)

# %%
