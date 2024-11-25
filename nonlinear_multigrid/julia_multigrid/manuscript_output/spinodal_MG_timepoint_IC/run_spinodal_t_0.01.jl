#%%

using DelimitedFiles
using Dates

################################################################################################################################
# Spinodal decomposition +1 and -1 starting from MG t = 0.01 (after FD mistake with increasing energy) to compare to FD
################################################################################################################################

include("../../CH_multigrid_solver.jl")
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_MG_timepoint_IC/"

# tol = 1e-6
for dt in [0.0001, 5e-5, 2.5e-5, 1.25e-5, 6.25e-6]
    tol = 1e-6
    nx = 128
    m = 8
    epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
    total_time = 0.06
    max_it = Int.(total_time / dt)
    date_time = now()
    name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_eps_$(epsilon)"
    phi = initialization_from_file("$(outdir)/IC/t=0.01_dt_1.25e-5.csv", nx, nx)
    time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, print_phi=true, suffix=name, check_dir=false)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
        writedlm(f, [date_time "spinodal_t_0.01_IC" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
    end
end



