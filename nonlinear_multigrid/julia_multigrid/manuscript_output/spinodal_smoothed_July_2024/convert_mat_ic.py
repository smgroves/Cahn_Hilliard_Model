import scipy.io
import numpy as np

# Load .mat file
mat = scipy.io.loadmat("./nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/t=1/SAV_smoothed_IC_time500_dt_6p25e-6.mat")

# Access matrix (assuming it's named 'matrix_data')
matrix_data = mat['phi']  # Replace 'matrix_data' with the actual matrix name

# Save as .txt
np.savetxt('./nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/SAV_smoothed_IC_time500_dt_6p25e-6.txt', matrix_data, delimiter=' ')