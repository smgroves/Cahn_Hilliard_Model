#%%
using DelimitedFiles
using Dates

# nx = 1024
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output"
# mean_ = 0
# std_dev = 0.2
# phi = mean .+ std_dev .* randn(nx, nx)
# writedlm("$(outdir)initial_phi_$(nx)_mean_$(mean)_sd_$(std_dev).csv", phi, ',')
#%%
################################
# Spinodal decomposition noise around 0
################################
using DelimitedFiles
using Dates
include("../../CH_multigrid_solver.jl")
# tol = 1e-6
for nx in [32, 64, 128, 256]
    for tol in [1e-4, 1e-5, 1e-6]
        for dt in [1.25e-5, 2.5e-5, 5e-5, 0.0001]
            m = 8
            epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
            # total_time = 0.03
            # max_it = Int.(total_time / dt)
            max_it = 4800
            date_time = now()
            phi = initialization_from_file("$(outdir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
            time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=false, suffix="dt_$(dt)_mean_0_sd_0.2", check_dir=false)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
                writedlm(f, [date_time "spinodal_normal_IC_no_printing" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
            end
        end
    end
end

#%%
using DelimitedFiles
using Dates
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output"

include("../../CH_multigrid_solver.jl")
# tol = 1e-6
for nx in [32, 64, 128, 256]
    for tol in [1e-4, 1e-5, 1e-6]
        for dt in [3e-9]
            m = 8
            epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
            # total_time = 0.003
            # max_it = Int.(total_time / dt)
            max_it = 4800
            date_time = now()
            phi = initialization_from_file("$(outdir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
            time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, ns=1000, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=false, suffix="dt_$(dt)_mean_0_sd_0.2", check_dir=false)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
                writedlm(f, [date_time "spinodal_normal_IC_no_printing" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
            end
        end
    end
end

#%%
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
nx = 256
tol = 0.0001
dt = 2.5e-5
max_it = 4800
date_time = now()
phi = initialization_from_file("$(outdir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=false, suffix="dt_$(dt)_mean_0_sd_0.2", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_normal_IC_no_printing" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end

#%%
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
nx = 256
tol = 1e-6
dt = 5e-5
max_it = 4800
date_time = now()
phi = initialization_from_file("$(outdir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=false, suffix="dt_$(dt)_mean_0_sd_0.2", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_normal_IC_no_printing" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end

#%% testing the number of vcycles
using DelimitedFiles
using Dates
include("../../CH_multigrid_solver.jl")
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output"

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output/dt_vs_tol_num_vcycles"
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
for nx in [64, 256]
    h = 1 / nx
    for tol in [1e-4, 1e-3, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
        # for dt in [1e-7 * h^2, 1e-6 * h^2, 1e-5 * h^2, 1e-4 * h^2, 1e-3 * h^2, 1e-2 * h^2, 1e-1 * h^2]
        for dt in [1e3 * h^2, 1e4 * h^2, 1e5 * h^2]
            max_it = 1
            date_time = now()
            phi = initialization_from_file("$(indir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
            time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, max_it_CH=10000, print_mass=false, print_e=false, overwrite=false, print_r=true, print_phi=false, suffix="residual_$(nx)_$(tol)_dt_$(dt)_mean_0_sd_0.2", check_dir=false)
            # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
            #     writedlm(f, [date_time "spinodal_normal_IC_1_iteration" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
            # end
        end
    end
end

#%%  comparing to FDM, not run yet
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
nx = 128
tol = 1e-6
dt = 5e-7
max_it = 60000
date_time = now()
phi = initialization_from_file("$(outdir)/initial_phi_$(nx)_mean_0_sd_0.2.csv", nx, nx)
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=true, suffix="dt_$(dt)_mean_0_sd_0.2", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_normal_IC_" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end
