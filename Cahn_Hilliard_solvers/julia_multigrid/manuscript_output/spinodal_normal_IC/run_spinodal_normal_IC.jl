#%%
using DelimitedFiles
using Dates

nx = 1024
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_normal_IC/"
mean = 0
std_dev = 0.2
phi = mean .+ std_dev .* randn(nx, nx)
writedlm("$(outdir)initial_phi_$(nx)_mean_$(mean)_sd_$(std_dev).csv", phi, ',')
#%%
################################
# Spinodal decomposition noise around 0
################################
include("../../CH_multigrid_solver.jl")
# tol = 1e-6
dt = 6.25e-6
for nx in [32, 64, 128, 256, 512]
    for tol in [1e-4, 1e-5, 1e-6, 1e-7]
        for dt in [6.25e-6, 1.25e-5, 2.5e-5, 5e-5, 0.0001]
            m = 8
            epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
            total_time = 0.03
            max_it = Int.(total_time / dt)
            date_time = now()
            phi = initialization_from_file("$(outdir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
            time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_mean_$(mean)_sd_$(std_dev)", check_dir=false)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, [date_time "spinodal_normal_IC" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
            end
        end
    end
end