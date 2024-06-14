#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
indir ="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius/"
tmp = pd.read_csv(f"{indir}/from_Rivanna/radius_0.5_level_set_epsilon_0.015009.txt",header = 0, index_col=None)
print(tmp.shape)
# %%
tmp['R0'].unique()

#%%
sns.lineplot(data = tmp, x = 'time', y = 'radius', hue = 'R0', palette = 'tab20')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title("Epsilon = 0.015009")
plt.show()
# %%
