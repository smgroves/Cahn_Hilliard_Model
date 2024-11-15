import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools


def read_specific_lines(file_path, line_numbers):
    result = []
    with open(file_path, "r") as file:
        for current_line_number, line in enumerate(file):
            if current_line_number in line_numbers:
                digits = [float(i) for i in line.strip("\n").split(" ")]
                result.append(digits)
            if current_line_number > max(line_numbers):
                break
    result = np.array(result)
    return result


def redblue(m=None):
    # If m is not specified, use the current figure's colormap size
    if m is None:
        m = plt.get_cmap().N

    if m % 2 == 0:
        # Even case: From [0, 0, 1] to [1, 1, 1], then [1, 1, 1] to [1, 0, 0]
        m1 = m // 2
        r = np.linspace(0, 1, m1)
        g = r
        r = np.concatenate((r, np.ones(m1)))
        g = np.concatenate((g, np.flipud(g)))
        b = np.flipud(r)
    else:
        # Odd case: From [0, 0, 1] to [1, 1, 1], then [1, 1, 1] to [1, 0, 0]
        m1 = m // 2
        r = np.linspace(0, 1, m1)
        g = r
        r = np.concatenate((r, np.ones(m1 + 1)))
        g = np.concatenate((g, [1], np.flipud(g)))
        b = np.flipud(r)

    # Combine r, g, b to create the colormap array
    c = np.vstack((r, g, b)).T
    return c


def plot_snapshot(
    indir, name, Nx, dt, dt_out=10, time=None, timepoint=None, save=False, outdir=""
):
    if time != None:
        timepoint = int(time / (dt * dt_out))
    elif timepoint != None:
        time = timepoint * dt
    else:
        raise Exception("Either time or timepoint must be set. time = timepoint * dt.")
    print(timepoint)
    first_line = timepoint * Nx
    last_line = first_line + Nx

    line_list = range(first_line, last_line)

    fig = plt.figure(figsize=(3, 4), dpi=300)
    snapshot = read_specific_lines(f"{indir}/{name}", line_list)
    sns.heatmap(
        snapshot[128:384].T,
        square=True,
        cmap=plt.cm.colors.ListedColormap(redblue(100)),
        xticklabels=100,
        yticklabels=100,
        linewidths=0,
    )
    plt.title(f"t = {time}")
    plt.tight_layout()
    if save:
        plt.savefig(f"{outdir}/t={time}_{name}.png")
        plt.close()
    else:
        plt.show()


indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
name = "phi_512_19661_1.0e-5__CPC_0.173_cohesin_0.1_eps_0.0075_alpha_0_domain_0_2.txt"
Nx = 512

# for time in [0.005, 0.01, 0.015, 0.02]:
plot_snapshot(
    indir,
    name,
    Nx,
    time=0.015106,
    dt=0.000001525878906,
    save=True,
    outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/plotting/manuscript/figure 6",
)
