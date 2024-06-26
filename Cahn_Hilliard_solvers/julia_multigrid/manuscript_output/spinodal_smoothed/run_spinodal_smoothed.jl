#%%
using Dates
include("../../CH_multigrid_solver.jl")
################################
# Spinodal decomposition smoothed t = 500
################################
nx = 128
tol = 1e-6
dt = 6.25e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/"
total_time = 0.02
phi = initialization_from_file("$outdir/SAV_smoothed_IC_time500_dt_6p25e-6.txt", nx, nx, delim=' ', transpose_matrix=true)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_eps_$(epsilon)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_smoothed" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end

################################
# Spinodal decomposition smoothed
################################
nx = 128
tol = 1e-6
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/"
total_time = 0.02
phi = initialization_from_file("$outdir/SAV_smoothed_IC_time500_dt_6p25e-6.txt", nx, nx, delim=' ', transpose_matrix=true)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, suffix="dt_$(dt)_eps_$(epsilon)_transposed", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_smoothed" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end

################################
# Spinodal decomposition smoothed
################################
nx = 128
tol = 1e-6
dt = 2.5e-5
m = 4
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/"
total_time = 0.02
phi = initialization_from_file("$outdir/SAV_smoothed_IC_time500_dt_6p25e-6.txt", nx, nx, delim=' ', transpose_matrix=true)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, suffix="dt_$(dt)_eps_$(epsilon)_transposed", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_smoothed" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end

#%%
using Dates
include("../../CH_multigrid_solver.jl")
################################
# Spinodal decomposition smoothed - t=1
################################
nx = 128
tol = 1e-6
dt = 6.25e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/t_5"
total_time = 0.02
phi = initialization_from_file("$outdir/5_dt_6p25e-6.txt", nx, nx, delim=' ', transpose_matrix=true)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_eps_$(epsilon)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_smoothed_t=1" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end

################################
# Spinodal decomposition smoothed
################################
nx = 128
tol = 1e-6
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/t_5"
total_time = 0.02
phi = initialization_from_file("$outdir/5_dt_6p25e-6.txt", nx, nx, delim=' ', transpose_matrix=true)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, suffix="dt_$(dt)_eps_$(epsilon)_transposed", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_smoothed_t=1" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end


#%%
################################
# Spinodal decomposition example - periodic boundary condition
################################
using Dates

include("../../CH_multigrid_solver_periodic_bc.jl")
nx = 128
tol = 1e-6
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/t_500"
total_time = 0.02
phi = initialization_from_file("$outdir/SAV_smoothed_IC_time500_dt_6p25e-6.txt", nx, nx, delim=' ', transpose_matrix=true)
println(dt)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver_with_periodic(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_periodic", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_periodic" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")

end

#%%
using Dates

include("../../CH_multigrid_solver_periodic_bc.jl")
nx = 128
tol = 1e-6
dt = 6.25e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/t_500"
total_time = 0.02
phi = initialization_from_file("$outdir/SAV_smoothed_IC_time500_dt_6p25e-6.txt", nx, nx, delim=' ', transpose_matrix=true)
println(dt)
max_it = Int.(total_time / dt)
date_time = now()
time_passed = @elapsed multigrid_solver_with_periodic(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_periodic", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_periodic" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")

end
