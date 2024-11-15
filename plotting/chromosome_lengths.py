import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv
lengths = []
indir = '/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/haspin_stripe_linescans'
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/plotting/image_analysis"
for image in [0, 1, 2, 5, 7, 8, 9]:
    tmp = pd.read_csv(f"{indir}/analysis/count_peaks_image{image}_.csv", header=0, index_col=0)
    for i,r in tmp.iterrows():
        lengths.append(r['l_arm_length'])
        lengths.append(r['r_arm_length'])
print(lengths)
# sns.histplot(lengths, bins = 30)
# plt.title(f"Histogram of chromosome arm lengths in images \n n = 7 images, {len(lengths)} chromosomes; mean = {round(np.mean(lengths),3)}")
# plt.savefig(f"{indir}/chr_arm_lengths_hist.png")

np.savetxt(f"{outdir}/chromosome_lengths.csv", lengths, fmt = "%f", delimiter=",")

