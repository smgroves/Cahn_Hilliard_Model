import numpy as np

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/"

# phi_name = "spinodal__same_total_time/phi_128_3200_1.0e-6_dt_6.25e-6.txt"

# #t = 500 = 50th matrix
# skip_header = 128*50 #(t=0 is included, skip t=0 to t=49)
# max_rows = 128
# phi = np.genfromtxt(f"{indir}/{phi_name}", skip_header=skip_header, max_rows=max_rows)

# np.savetxt(f"{indir}/spinodal_smoothed/smoothed_with_multigrid/500_dt_6p25e-6.txt",phi)


phi_name = "spinodal__same_total_time/phi_128_320_1.0e-6_dt_6.25e-6_everytimestep.txt"
#t = 5 = 5th matrix
skip_header = 128*1 #(t=0 is included, skip t=0 to t=49)
max_rows = 128
phi = np.genfromtxt(f"{indir}/{phi_name}", skip_header=skip_header, max_rows=max_rows)

np.savetxt(f"{indir}/spinodal_smoothed/smoothed_with_multigrid/1_dt_6p25e-6.txt",phi)