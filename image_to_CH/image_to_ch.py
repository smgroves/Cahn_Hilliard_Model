# %%
from PIL import Image

# indir = (
#     "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting/manuscript/figure 5/"
# )
# im = Image.open(f"{indir}/T3_overlay_cropped_2_v2_CPC.jpg")
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/image_to_CH/input"
im = Image.open(f"{indir}/T6_I8_Lng-cropped-v2_CPC.jpg")
im.show()

res = im.resize((256, 128), resample=Image.NEAREST)


# %%
import numpy as np

# Convert the image to grayscale
img_gray = res.convert("L")

# Convert the image to a numpy array
img_array = np.array(img_gray)
# %%

import seaborn as sns
import matplotlib.pyplot as plt

sns.heatmap(img_array, square=True)
plt.show()
# %%
# rescale to 0 to 1
img_array_psi = (img_array - np.min(img_array)) / (
    np.max(img_array) - np.min(img_array)
)
# rescale to phi  = -1 to 1
im_array_phi = 2 * img_array_psi - 1

# %%
sns.heatmap(im_array_phi, square=True)
plt.show()

# %%
im_array_phi = im_array_phi.T


# %%
def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get("padder", 0)
    vector[: pad_width[0]] = pad_value
    if pad_width[1] != 0:  # <-- the only change (0 indicates no padding)
        vector[-pad_width[1] :] = pad_value


width = np.abs(im_array_phi.shape[0] - im_array_phi.shape[1])
im_array_phi_padded = np.pad(
    im_array_phi, ((0, 0), (int(width / 2), int(width / 2))), pad_with, padder=-1
)
sns.heatmap(im_array_phi_padded, square=True, cmap="RdBu")
plt.show()
# %%
np.savetxt(
    "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/image_to_CH/input/MCF10A_T6I8_chromosome_phi_IC_256.csv",
    im_array_phi_padded,
    delimiter=",",
)

# %%
