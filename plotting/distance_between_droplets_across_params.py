# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting"


def make_long_dist_df(indir, file):
    # dist = pd.read_csv(f"{indir}/distance_between_droplets.csv",converters={'distances': pd.eval})
    dist = pd.read_csv(f"{indir}/{file}", header=0)
    long_dist_df = pd.DataFrame(
        # columns=["seed", "cpc", 'cohesin', 'time', 'distance'])
        columns=["cpc", "cohesin", "distance"]
    )

    for i, r in dist.iterrows():
        dist_list = str(r["distances"])[1:-2].split(" ")
        try:
            dist_list = [float(d) for d in dist_list]
            for d in dist_list:
                long_dist_df = pd.concat(
                    [
                        long_dist_df,
                        pd.DataFrame(
                            {
                                # "seed": [r["seed"]],
                                # "time": [r["time"]],
                                "cpc": [r["cpc"]],
                                "cohesin": [r["cohesin"]],
                                "distance": [d],
                            }
                        ),
                    ],
                    ignore_index=True,
                )
        except ValueError:
            pass
    return long_dist_df


# %%

long_dist_df = make_long_dist_df(
    indir,
    "simulated_droplet_distributions/simulated_droplet_distances_e_0.01_noisy_cohesin_chr_lengths.csv",
)
sns.histplot(
    data=long_dist_df,
    x="distance",
    palette=sns.color_palette("muted"),
    binwidth=0.1,
    stat="probability",
    common_norm=False,
    kde=True,
)

# sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
plt.title("Distances between droplets by CPC radius and Cohesin width")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.tight_layout()
plt.show()
# plt.savefig(f"{outdir}distances_between_droplets_swarmplot_0.0075.png")
# %%
# outdir = (
#     f"{indir}/radii_lineplots_kymographs/domain_0_2_from_rivanna_kymographs_e_0.0075/"
# )

# long_dist_df = make_long_dist_df(
#     indir,
#     file="simulated_droplet_distributions/simulated_droplet_distances_e_0.0075_domain_0_2_chr_lengths.csv",
# )
# sns.histplot(
#     data=long_dist_df, x="distance", palette=sns.color_palette("muted"), bins=30
# )

# # sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
# plt.title("Distances between droplets by CPC radius and Cohesin width")
# plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
# plt.tight_layout()
# plt.show()
# %%
# outdir = f"{indir}/radii_lineplots_kymographs/domain_0_1_eps_0.0075/"

# long_dist_df2 = make_long_dist_df(
#     indir,
#     file="simulated_droplet_distributions/simulated_droplet_distances_e_0.0075_domain_0_1_chr_lengths.csv",
# )
# sns.histplot(
#     data=long_dist_df2, x="distance", palette=sns.color_palette("muted"), bins=30
# )

# # sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
# plt.title("Distances between droplets by CPC radius and Cohesin width")
# plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
# plt.tight_layout()
# plt.show()

# %%

long_dist_df2 = make_long_dist_df(
    indir,
    file="simulated_droplet_distributions/simulated_droplet_distances_e_0.0075_noisy_cohesin_sd_0.11.csv",
)
# %%
outdir = f"{indir}/radii_lineplots_kymographs/domain_0_2_noisy_cohesin_sd_0.11"

sns.histplot(
    data=long_dist_df2,
    x="distance",
    palette=sns.color_palette("muted"),
    binwidth=0.05,
    stat="probability",
    common_norm=False,
    kde=True,
)
plt.xlabel("Distance (um)")
plt.ylim(0, 0.16)
plt.xlim(0, 3.2)
# plt.tight_layout()
# plt.savefig(f"{outdir}/sim_histplot_eps_0.0075_bin_0.05.pdf")
plt.show()

# %%
sns.swarmplot(
    data=long_dist_df,
    x="cohesin",
    y="distance",
    hue="cpc",
    palette=sns.color_palette("muted"),
    size=4,
)
plt.title("Distances between droplets by CPC radius and Cohesin width")
plt.ylim(0, 0.16)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.tight_layout()
plt.show()
# %%

long_dist_df3 = make_long_dist_df(
    indir,
    file="simulated_droplet_distributions/simulated_droplet_distances_e_0.005_noisy_cohesin_chr_lengths.csv",
)
sns.histplot(
    data=long_dist_df3,
    x="distance",
    palette=sns.color_palette("muted"),
    binwidth=0.1,
    stat="probability",
    common_norm=False,
    kde=True,
)

# sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
plt.title("Distances between droplets by CPC radius and Cohesin width")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.tight_layout()
plt.show()

# %%

long_dist_df4 = make_long_dist_df(
    indir,
    file="simulated_droplet_distributions/simulated_droplet_distances_e_0.008_noisy_cohesin_chr_lengths.csv",
)
sns.histplot(
    data=long_dist_df4,
    x="distance",
    palette=sns.color_palette("muted"),
    binwidth=0.1,
    stat="probability",
    common_norm=False,
    kde=True,
)

# sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
plt.title("Distances between droplets by CPC radius and Cohesin width")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.tight_layout()
plt.show()

# %%

long_dist_df5 = make_long_dist_df(
    indir,
    file="simulated_droplet_distributions/simulated_droplet_distances_e_0.0067_noisy_cohesin_chr_lengths.csv",
)
sns.histplot(
    data=long_dist_df5,
    x="distance",
    palette=sns.color_palette("muted"),
    binwidth=0.1,
    stat="probability",
    common_norm=False,
    kde=True,
)

# sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
plt.title("Distances between droplets by CPC radius and Cohesin width")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
plt.tight_layout()
plt.show()

# %% compare two simulation types
outdir = f"{indir}/radii_lineplots_kymographs/domain_0_2_noisy_cohesin_sd_0.11"

binwidth = 0.05
long_dist_df["Category"] = "eps = 0.01"
long_dist_df2["Category"] = "eps = 0.0075"
# long_dist_df3["Category"] = "eps = 0.005"
long_dist_df4["Category"] = "eps = 0.008"
long_dist_df5["Category"] = "eps = 0.0067"

long_dist_df = pd.concat(
    [long_dist_df, long_dist_df2, long_dist_df4, long_dist_df5], ignore_index=True
)
sns.histplot(
    data=long_dist_df,
    x="distance",
    hue="Category",
    binwidth=binwidth,
    stat="probability",
    common_norm=False,
    kde=True,
)
# plt.savefig(f"{outdir}/distances_between_droplets_histplot_image_vs_sim_{binwidth}_v3.png")
plt.show()
# %%
# outdir = f"{indir}/radii_over_time_level_set_plots/domain_0_2_noisy_cohesin_sd_0.25/"

# long_dist_df = make_long_dist_df(
#     indir,
#     file="simulated_droplet_distances/simulated_droplet_distances_e_0.0075_noisy_cohesin_sd_0.25.csv",
# )
# sns.histplot(
#     data=long_dist_df, x="distance", palette=sns.color_palette("muted"), bins=30
# )

# # sns.swarmplot(data= long_dist_df, x = 'cohesin', y = 'distance', hue = 'cpc', palette=sns.color_palette("muted"), size = 4)
# plt.title("Distances between droplets by CPC radius and Cohesin width")
# plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0)
# plt.tight_layout()
# plt.show()

# %%
# TODO: redo this section to use bootstrapping on all images


# the below code is pulled from distances_between_droplets.py in the CPC_condensate_images folder on Box.
def inter_droplet_distance(indir, image):
    distance_dict = {}
    tmp = pd.read_csv(
        f"{indir}/count_peaks_image{image}_.csv",
        header=0,
        index_col=0,
        converters={"IC_peaks": pd.eval, "left_peaks": pd.eval, "right_peaks": pd.eval},
    )
    for i, r in tmp.iterrows():
        ic = list(r["IC_peaks"])
        left = list(r["left_peaks"])
        right = list(r["right_peaks"])
        [ic.extend(l) for l in (left, right)]
        all_peaks = sorted(ic)

        distances = []
        for j in range(len(all_peaks) - 1):
            d = (all_peaks[j + 1] - all_peaks[j]) * 0.06013
            distances.append(d)
        distance_dict[i] = distances
    return distance_dict


indir2 = "/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/haspin_stripe_linescans/analysis"
all_images = []
for image in [0, 1, 2, 5, 7, 8, 9]:
    distance_dict = inter_droplet_distance(indir2, image=image)

    all_ = []
    for k in distance_dict.keys():
        all_.extend(distance_dict[k])
    all_images.extend(all_)

print(all_images)

# long_dist_df["Category"] = "Simulation"
for d in all_images:
    long_dist_df = pd.concat(
        [
            long_dist_df,
            pd.DataFrame(
                {
                    "Category": ["Experiment_HeLa"],
                    "seed": [0],
                    "time": [0],
                    "cpc": [0],
                    "cohesin": [0],
                    "distance": [d],
                }
            ),
        ],
        ignore_index=True,
    )

# %%
indir3 = "/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/haspin_stripe_linescans/MCF10A"
all_images = []
for image in [10]:
    distance_dict = inter_droplet_distance(indir3, image=image)

    all_ = []
    for k in distance_dict.keys():
        all_.extend(distance_dict[k])
    all_images.extend(all_)

print(all_images)

# long_dist_df["Category"] = "Simulation"
for d in all_images:
    long_dist_df = pd.concat(
        [
            long_dist_df,
            pd.DataFrame(
                {
                    "Category": ["Experiment_MCF10A"],
                    "seed": [0],
                    "time": [0],
                    "cpc": [0],
                    "cohesin": [0],
                    "distance": [d],
                }
            ),
        ],
        ignore_index=True,
    )

# %%
# Earth movers distance
import scipy.stats as ss
from itertools import combinations

distances = pd.DataFrame(
    np.zeros((5, 5)),
    index=long_dist_df["Category"].unique(),
    columns=long_dist_df["Category"].unique(),
)

for i, j in list(combinations(long_dist_df["Category"].unique(), 2)):
    i_dist = long_dist_df.loc[long_dist_df["Category"] == i]["distance"].values
    j_dist = long_dist_df.loc[long_dist_df["Category"] == j]["distance"].values

    distance = ss.wasserstein_distance(i_dist, j_dist)
    distances.loc[j, i] = distance

mask = np.triu(np.ones_like(distances))

# plotting a triangle correlation heatmap
# %%
dataplot = sns.heatmap(distances, cmap="Reds_r", annot=True, mask=mask)
plt.tight_layout()
plt.show()
# plt.savefig(f"{outdir}/wasserstein_distance_image_vs_sim_wMCF10A.png")
# plt.show()
# %%
# dist_df = pd.DataFrame(columns=['distance', 'chr'])
binwidth = 0.05
sns.histplot(
    data=long_dist_df,
    x="distance",
    hue="Category",
    binwidth=binwidth,
    stat="probability",
    common_norm=False,
    kde=True,
)

# num_sim_timepoints = len(long_dist_df[['seed', 'cpc', 'cohesin', 'time']].value_counts().reset_index().index) -1
num_sim_timepoints = (
    (len(long_dist_df["seed"].unique()) - 1)
    * (len(long_dist_df["cpc"].unique()) - 1)
    * (len(long_dist_df["cohesin"].unique()) - 1)
    * 254
)
# num_sim_timepoints = 33
plt.title(
    f"Distances between droplets: 7 Images (n = {len(long_dist_df.loc[long_dist_df['Category']=='Experiment'].index)}) \nvs {num_sim_timepoints} Simulation Timepoints ($n_{{0.01}}$ = {len(long_dist_df.loc[long_dist_df['Category']=='eps = 0.01'].index)},$n_{{0.0075}}$ = {len(long_dist_df.loc[long_dist_df['Category']!='eps = 0.0075'].index)})"
)

plt.xlabel("Distance (um)")
plt.xlim(0, 3.2)
plt.ylim(0, 0.2)
# plt.savefig(
#     f"{outdir}/distances_between_droplets_histplot_image_vs_sim_{binwidth}_v3.png"
# )
plt.show()
# plt.close()

# sns.histplot(
#     data=long_dist_df.loc[long_dist_df["Category"] == "Experiment"],
#     x="distance",
#     bins=30,
#     stat="probability",
#     binwidth=binwidth,
#     common_norm=False,
#     color=sns.color_palette()[1],
#     kde=True,
# )
# plt.xlim(0, 3.2)
# # plt.ylim(0, 0.16)
# plt.title(
#     f"Distances between droplets\n 7 Images, median = {round(long_dist_df.loc[long_dist_df['Category']== 'exp']['distance'].median(),3)}"
# )
# plt.xlabel("Distance (um)")
# # plt.savefig(f"{outdir}distances_between_droplets_histplot_image_{binwidth}.png")
# plt.show()
# %%
sns.kdeplot(
    data=long_dist_df,
    x="distance",
    hue="Category",
    common_norm=False,
)
plt.show()

# %%
outdir = f"{indir}/radii_lineplots_kymographs/domain_0_2_noisy_cohesin_sd_0.11"

sns.kdeplot(
    data=long_dist_df,
    x="distance",
    hue="Category",
    common_norm=False,
)
plt.savefig(
    f"{outdir}/distances_between_droplets_kde_image_vs_sim_{binwidth}_withMCF10A.png"
)


# plt.show()


# %%
## Plot images by location of distance (IC-neighbor distance, outside IC-outside IC distance, or IC-IC distance)
def inter_droplet_distance_with_meta(indir, image):
    tmp = pd.read_csv(
        f"{indir}/count_peaks_image{image}_.csv",
        header=0,
        index_col=0,
        converters={"IC_peaks": pd.eval, "left_peaks": pd.eval, "right_peaks": pd.eval},
    )
    distances = pd.DataFrame(columns=["distance", "i"])
    for i, r in tmp.iterrows():
        ic = list(r["IC_peaks"])
        left = list(r["left_peaks"])
        right = list(r["right_peaks"])
        [ic.extend(l) for l in (left, right)]
        all_peaks = sorted(ic)

        for j in range(len(all_peaks) - 1):
            d = (all_peaks[j + 1] - all_peaks[j]) * 0.06013
            if (all_peaks[j + 1] in r["IC_peaks"]) and (all_peaks[j] in r["IC_peaks"]):
                IC_distance = "inner-IC distance"
            elif (all_peaks[j + 1] in r["IC_peaks"]) or (all_peaks[j] in r["IC_peaks"]):
                IC_distance = "IC distance"
            else:
                IC_distance = "Outside IC"
            distances = pd.concat(
                [
                    distances,
                    pd.DataFrame(
                        {"distance": [d], "IC_distance": [IC_distance], "i": [i]}
                    ),
                ],
                ignore_index=True,
            )
    return distances


indir2 = "/Users/smgroves/Box/CPC_Model_Project/CPC_condensate_images/haspin_stripe_linescans/"
all_distances = pd.DataFrame(columns=["distance", "IC_distance", "i"])
for image in [0, 1, 2, 5, 7, 8, 9]:
    all_distances = pd.concat(
        [all_distances, inter_droplet_distance_with_meta(indir2, image=image)],
        ignore_index=True,
    )

for i, d in long_dist_df.iterrows():
    all_distances = pd.concat(
        [
            all_distances,
            pd.DataFrame({"distance": d["distance"], "IC_distance": ["sim"], "i": 0}),
        ],
        ignore_index=True,
    )


means = all_distances.groupby("IC_distance").distance.mean()

print(means)
sns.kdeplot(data=all_distances, x="distance", hue="IC_distance", common_norm=True)
# plt.axvline(x=means.loc['IC distance'], c=sns.color_palette("deep")[1], linestyle="--")
# plt.axvline(x = means.loc['Outside IC'], c=sns.color_palette("deep")[0], linestyle = "--")
# plt.axvline(x = means.loc['inner-IC distance'], c=sns.color_palette("deep")[2], linestyle = "--")

plt.title(
    f"Distances between droplets\n All condensates from 7 images, separated by proximity to inner centromere"
)
plt.xlabel("Distance (um)")
plt.ylabel("Frequency")
# plt.show()
plt.savefig(f"{outdir}distances_between_droplets_kdeplot_image_vs_sim_grouped.png")

# %%
