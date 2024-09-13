import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/"
outdir = f"{indir}/radii_over_time_level_set_plots/domain_0_2_e_0.0075_noisy_cohesin/"

# dist = pd.read_csv(f"{indir}/distance_between_droplets.csv",converters={'distances': pd.eval})
dist = pd.read_csv(f"{indir}/simulated_droplet_distances_e_0.0075_noisy_cohesin.csv",
                   header=0)
long_dist_df = pd.DataFrame(columns=["cpc",'cohesin','distance'])

for i, r in dist.iterrows():
    dist_list = str(r['distances'])[1:-2].split(" ")
    try:
        dist_list = [float(d) for d in dist_list]
        for d in dist_list:
            long_dist_df = pd.concat([
                long_dist_df,
                pd.DataFrame({
                    "cpc": [r["cpc"]],
                    "cohesin": [r["cohesin"]],
                    "distance": [d]
                })
            ],
                                    ignore_index=True)
    except ValueError: pass

# long_dist_df = pd.DataFrame(columns=["cpc",'cohesin','distance'])
# for i,r in dist.iterrows():
#     for d in r['distances']:
#         long_dist_df = pd.concat([long_dist_df,pd.DataFrame({"cpc":[r["cpc"]], "cohesin":[r["cohesin"]], "distance":[d]})], ignore_index=True)

sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
plt.title("Distances between droplets by CPC radius and Cohesin width")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
plt.savefig(f"{outdir}distances_between_droplets_swarmplot_0.0075_noisy_cohesin.png")


# the below code is pulled from distances_between_droplets.py in the CPC_condensate_images folder on Box.
def inter_droplet_distance(indir, image):
    distance_dict = {}
    tmp = pd.read_csv(f"{indir}/count_peaks_image{image}_.csv", header = 0, index_col = 0,
        converters={'IC_peaks': pd.eval, "left_peaks":pd.eval, "right_peaks":pd.eval})
    for i,r in tmp.iterrows():
        ic = list(r['IC_peaks'])
        left = list(r['left_peaks'])
        right = list(r['right_peaks'])
        [ic.extend(l) for l in (left,right)]
        all_peaks = sorted(ic)

        distances = []
        for j in range(len(all_peaks)-1):
            d = (all_peaks[j+1] - all_peaks[j])*.06013
            distances.append(d)
        distance_dict[i] = distances
    return distance_dict

indir2 = "/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/haspin_stripe_linescans/"
all_images = []
for image in [0,1,2,5,7,8,9]:
    distance_dict = inter_droplet_distance(indir2, image=image)

    all_ = []
    for k in distance_dict.keys():
        all_.extend(distance_dict[k])
    all_images.extend(all_)

print(all_images)

long_dist_df['Category'] = 'sim'
for d in all_images:
    long_dist_df = pd.concat([long_dist_df,pd.DataFrame({"Category":["exp"],"cpc":[0], "cohesin":[0], "distance":[d]})], ignore_index=True)

# dist_df = pd.DataFrame(columns=['distance', 'chr'])
sns.histplot(data= long_dist_df, x = 'distance', hue = 'Category', bins = 30)
plt.title(f"Distances between droplets\n 7 Images versus Simulations (eps = 0.0096)")
plt.xlabel("Distance (um)")
plt.savefig(f"{outdir}distances_between_droplets_histplot_image_vs_sim_.png")
plt.show()

sns.histplot(data= long_dist_df.loc[long_dist_df['Category']== 'exp'], x = 'distance', bins = 30)
plt.title(f"Distances between droplets\n 7 Images, median = {round(long_dist_df.loc[long_dist_df['Category']== 'exp']['distance'].median(),3)}")
plt.xlabel("Distance (um)")
plt.savefig(f"{outdir}distances_between_droplets_histplot_image_.png")
plt.show()

def inter_droplet_distance_with_meta(indir, image):
    tmp = pd.read_csv(f"{indir}/count_peaks_image{image}_.csv",
                      header=0,
                      index_col=0,
                      converters={
                          'IC_peaks': pd.eval,
                          "left_peaks": pd.eval,
                          "right_peaks": pd.eval
                      })
    distances = pd.DataFrame(
            columns=["distance", "i"])
    for i, r in tmp.iterrows():
        ic = list(r['IC_peaks'])
        left = list(r['left_peaks'])
        right = list(r['right_peaks'])
        [ic.extend(l) for l in (left, right)]
        all_peaks = sorted(ic)

        for j in range(len(all_peaks) - 1):
            d = (all_peaks[j + 1] - all_peaks[j]) * .06013
            if (all_peaks[j + 1] in r['IC_peaks']) and (all_peaks[j]
                                                        in r['IC_peaks']):
                IC_distance = "inner-IC distance"
            elif (all_peaks[j + 1] in r['IC_peaks']) or (all_peaks[j]
                                                         in r['IC_peaks']):
                IC_distance = "IC distance"
            else:
                IC_distance = "Outside IC"
            distances = pd.concat([
                distances,
                pd.DataFrame({
                    "distance": [d],
                    "IC_distance": [IC_distance],
                    "i": [i]
                })
            ],
                                  ignore_index=True)
    return distances

# indir2 = "/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/haspin_stripe_linescans/"
# all_distances = pd.DataFrame(columns=["distance", "IC_distance", "i"])
# for image in [0,1,2,5,7,8,9]:
#     all_distances = pd.concat([all_distances,inter_droplet_distance_with_meta(indir2, image=image)], ignore_index = True)

# for i, d in long_dist_df.iterrows():
#     all_distances = pd.concat([all_distances,pd.DataFrame({"distance":d['distance'], "IC_distance":["sim"], "i":0})], ignore_index=True)


# means = all_distances.groupby("IC_distance").distance.mean()

# print(means)
# sns.kdeplot(data=all_distances, x='distance', hue="IC_distance", common_norm=True)
# # plt.axvline(x=means.loc['IC distance'], c=sns.color_palette("deep")[1], linestyle="--")
# # plt.axvline(x = means.loc['Outside IC'], c=sns.color_palette("deep")[0], linestyle = "--")
# # plt.axvline(x = means.loc['inner-IC distance'], c=sns.color_palette("deep")[2], linestyle = "--")

# plt.title(f"Distances between droplets\n All condensates from 7 images, separated by proximity to inner centromere")
# plt.xlabel("Distance (um)")
# plt.ylabel("Frequency")
# # plt.show()
# plt.savefig(
#     f"{outdir}distances_between_droplets_kdeplot_image_vs_sim_grouped.png")
