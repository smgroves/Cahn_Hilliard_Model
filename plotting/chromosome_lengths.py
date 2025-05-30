# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import csv
lengths = []
total_lengths = []
indir = '/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/MCF10A/MCF10A_CPC_analysis/analysis/linescans'
# indir = "/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/haspin_stripe_linescans/analysis/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting/image_analysis"
cell_line = "MCF10A"
image_cnt = 0
for image in range(50):
    try:
        tmp = pd.read_csv(
            f"{indir}/count_peaks_image{image}_.csv", header=0, index_col=0)
        for i, r in tmp.iterrows():
            lengths.append(r['l_arm_length'])
            lengths.append(r['r_arm_length'])
            total_lengths.append(r['l_arm_length'] + r['r_arm_length'])
        image_cnt += 1
    except FileNotFoundError:
        print(f"Image {image} not found")
sns.histplot(lengths, bins=30)
plt.axvline(np.percentile(lengths, 95), color='red',
            linestyle='--', label='95th percentile')
plt.axvline(np.mean(lengths), color='grey', linestyle='--', label='mean')
plt.legend()
plt.title(
    f"Histogram of chromosome arm lengths in images \n n = {image_cnt} images, {len(lengths)//2} chromosomes; mean = {round(np.mean(lengths),3)} (um, one-arm length)")
plt.savefig(f"{outdir}/chr_arm_lengths_hist_{cell_line}.pdf")
plt.close()
sns.histplot(total_lengths, bins=30)
plt.axvline(np.percentile(total_lengths, 95), color='red',
            linestyle='--', label='95th percentile')
plt.axvline(np.mean(total_lengths), color='grey', linestyle='--', label='mean')
plt.legend()
plt.title(
    f"Histogram of chromosome sizes in images \n n = {image_cnt} images, {len(total_lengths)} chromosomes; mean = {round(np.mean(total_lengths),3)} (um, full length)")
# plt.show()
plt.savefig(f"{outdir}/chr_full_lengths_hist_{cell_line}.pdf")
plt.close()
# np.savetxt(f"{outdir}/chromosome_lengths_{cell_line}_v2.csv",
#            lengths, fmt="%f", delimiter=",")

# %%
