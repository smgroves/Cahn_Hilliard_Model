import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv
lengths = []
indir = '/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/MCF10A/MCF10A_CPC_analysis/analysis/linescans'
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting/image_analysis"
for image in range(1, 51):
    try:
        tmp = pd.read_csv(
            f"{indir}/count_peaks_image{image}_.csv", header=0, index_col=0)
        for i, r in tmp.iterrows():
            lengths.append(r['l_arm_length'])
            lengths.append(r['r_arm_length'])
    except FileNotFoundError:
        print(f"Image {image} not found")
sns.histplot(lengths, bins=30)
plt.title(
    f"Histogram of chromosome arm lengths in images \n n = 7 images, {len(lengths)} chromosomes; mean = {round(np.mean(lengths),3)}")
plt.savefig(f"{outdir}/chr_arm_lengths_hist_MCF10A.png")

np.savetxt(f"{outdir}/chromosome_lengths_MCF10A.csv",
           lengths, fmt="%f", delimiter=",")
