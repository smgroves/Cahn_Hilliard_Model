import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
r = pd.read_csv('/Users/smgroves/Documents/GitHub/jlCHSolver/residuals_256.txt', sep = "\t", header = None)
r.columns = ["residuals"]
r = r.fillna(1000)
print(np.max(r.residuals))
print(r.head())
print(len(r.index))
p = sns.lineplot(data = r, y = "residuals", x = range(len(r.index)))
p.set(yscale='log')
plt.axhline(y = 1e-6, linestyle="--", color = 'grey')
plt.ylabel("Residual for each iteration")
plt.xlabel("Iteration #")
# plt.ylim([1e-7, 1e5])
plt.title("Residuals for Solver on 256x256 Grid")

plt.savefig("../plots/residuals_256.png")
# plt.show()