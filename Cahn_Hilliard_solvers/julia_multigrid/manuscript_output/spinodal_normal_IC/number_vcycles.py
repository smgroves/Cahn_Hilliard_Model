import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def calculate_vcycles(nx, tol, dt, indir, one_timestep=True):
    num_cycles = []
    count = 0
    highest_r = 100000
    num_timesteps = 0
    with open(f"{indir}/residual_{nx}_{tol}_dt_{dt}_mean_0_sd_0.2.txt") as file:
        if one_timestep:
            for count, line in enumerate(file):
                pass
            ave_num_cycles = count
        else:
            for l, line in enumerate(file):
                r = float(line.rstrip())
                if r < highest_r:
                    highest_r = r
                    count += 1
                else:
                    num_cycles.append(count)
                    highest_r = 100000
                    count = 1
                    num_timesteps += 1
                # if num_timesteps == 50:
                # break
            ave_num_cycles = math.ceil(np.mean(num_cycles))
    return ave_num_cycles


indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_normal_IC/output/dt_vs_tol_num_vcycles"

table = pd.DataFrame(columns=["nx", "tol", "dt", "num_vcycles"])

for nx in [64, 12, 256]:
    # for tol in ["0.0001", "1.0e-5", "1.0e-6", "1.0e-7", "1.0e-8"]:
    # for dt in ["6.25e-6", "1.25e-5", "2.5e-5", "5.0e-5", "0.0001"]:
    h = 1 / nx
    for tol in ["1.0e-5", "1.0e-6", "1.0e-7", "1.0e-8", "1.0e-9", "1.0e-10"]:
        for dt_multiplier in ["-6", "-7", "-8", "-9", "-10", "-11", "-12"]:

            try:
                dt = str(h**2 * 1e5) + "e" + dt_multiplier
                ans = calculate_vcycles(nx, tol, dt, indir)
                table = pd.concat(
                    [
                        table,
                        pd.DataFrame(
                            {"nx": [nx], "tol": [tol], "dt": [dt], "num_vcycles": [ans]}
                        ),
                    ],
                    ignore_index=True,
                )
            except FileNotFoundError:
                pass

# print(table)
print(pd.pivot(table, columns=["nx", "dt"], index="tol"))

# print(pd.pivot_table(table, columns=["dt"], index="tol", aggfunc="mean"))
