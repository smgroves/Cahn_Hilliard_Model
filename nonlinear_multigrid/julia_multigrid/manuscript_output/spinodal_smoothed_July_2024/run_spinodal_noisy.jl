#%%
using Dates
include("../../CH_multigrid_solver.jl")
################################
# Spinodal decomposition noise around 0
################################
nx = 128
tol = 1e-6
dt = 6.25e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/noisy_IC"
total_time = 0.02
mean = -0.58
std_dev = 0.8
phi = mean .+ std_dev .* randn(nx, nx)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_eps_$(epsilon)_mean_$(mean)_sd_$(std_dev)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_smoothed" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end