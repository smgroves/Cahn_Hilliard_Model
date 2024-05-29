import matplotlib.pyplot as plt
import numpy as np

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/figure_2_spinodal_decomp"
name = "phi_16_10_1.0e-5_SD"

arr_j = np.genfromtxt(f"{indir}/{name}.txt",
                      delimiter=" ",
                      max_rows=512)

arr_j3d = arr_j.reshape(2, 256, 256)