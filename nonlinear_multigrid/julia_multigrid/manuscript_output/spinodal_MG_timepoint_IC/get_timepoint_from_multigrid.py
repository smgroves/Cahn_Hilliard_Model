import numpy as np

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_+1_-1_IC/output"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_MG_timepoint_IC/IC"

# phi_name = "spinodal__same_total_time/phi_128_3200_1.0e-6_dt_6.25e-6.txt"

# #t = 500 = 50th matrix
# skip_header = 128*50 #(t=0 is included, skip t=0 to t=49)
# max_rows = 128
# phi = np.genfromtxt(f"{indir}/{phi_name}", skip_header=skip_header, max_rows=max_rows)

# np.savetxt(f"{indir}/spinodal_smoothed/smoothed_with_multigrid/500_dt_6p25e-6.txt",phi)


phi_name = "phi_128_4800_1.0e-6_dt_1.25e-5.txt"
# want t = 0.01 -> 800 timesteps = 80th matrix (every ten steps recorded)
skip_header = 128 * 80  # (t=0 is included, skip t=0 to t=79)
max_rows = 128
phi = np.genfromtxt(f"{indir}/{phi_name}", skip_header=skip_header, max_rows=max_rows)

np.savetxt(f"{outdir}/t=0.01_dt_1.25e-5.txt", phi)
