import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(
    "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/t_1/discrete_norm_e_128_3200_1.0e-6_dt_6.25e-6_eps_0.015009369912862116.txt",
    header = 0
)

df.columns = ["e"]
plt.plot(df["e"])
plt.title("Normalized Energy for Multigrid Solver")
plt.ylabel("Normalized Energy")
plt.xlabel("Timesteps")
plt.savefig("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/manuscript/figure 2/discrete_norm_energy.pdf")
