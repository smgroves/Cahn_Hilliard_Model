
using DelimitedFiles
using Dates
using Plots

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output/"

include("../../CH_multigrid_solver_large_domain.jl")
#domain = 0 to 128
nx = 128
dt = 1e-1 / (nx^2)
n_relax = 4
ny = nx
tol = 1e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
epsilon

gamma = epsilon^2 / h^2;
# total_time = 0.01
total_time = (0.1 / 64^2) * 500
max_it = Int.(total_time / dt)
date_time = now()
name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_n_relax_$(n_relax)_eps_$(epsilon)"
phi = initialization_from_file("$(indir)initial_phi_$(nx)_smooth_n_relax_$(n_relax).csv", nx, nx)
final_phi = multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=true, suffix=name, check_dir=false)
writedlm("$(outdir)/$(name)_final_phi.csv", final_phi, ',')
