using DelimitedFiles
using Dates
using BenchmarkTools

include("../../CH_multigrid_solver_optimized.jl")
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/optimized_testing/"

nx = 64
tol = 1e-6
dt = 5e-5
m = 8
mean = 0
std_dev = 0.2
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
total_time = 0.01
max_it = Int.(total_time / dt)
date_time = now()
phi = initialization_from_file("$(indir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
@profview multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, suffix="dt_$(dt)_mean_$(mean)_sd_$(std_dev)", check_dir=false)
# open("$(outdir)Job_specs_mem_allocation_optimization.csv", "a", lock=false) do f
#     writedlm(f, [date_time "5A - threaded" "Julia" nx epsilon dt tol max_it 10000 mem_allocated / 1e6], ",")
# end
