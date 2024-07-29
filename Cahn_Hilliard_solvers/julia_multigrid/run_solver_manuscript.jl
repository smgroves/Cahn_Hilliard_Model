#%%
# using Pkg
# Pkg.add("DataFrames")
# Pkg.add("BenchmarkTools")
# Pkg.add("StaticArrays")
# Pkg.add("ProfileView")
# Pkg.add("DelimitedFiles")
# Pkg.add("LinearAlgebra")
# Pkg.add("Printf")
using Dates
include("./CH_multigrid_solver.jl")

#%%
#########################
# Single droplet example
#########################
nx = 128
tol = 1e-6
dt_0 = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
dts = [0.25 * dt_0, 0.5 * dt_0, dt_0, 2 * dt_0, 4 * dt_0]
total_time = 1
phi = initialization(nx, nx, method="droplet", h=1 / 128, epsilon=epsilon, R0=0.2)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/droplet_R0_0.2_v3"
for dt in dts
    println(dt)
    max_it = Int.(total_time / dt)
    date_time = now()
    max_it_CH = 10000
    time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, suffix="_dt_$dt", check_dir=false)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, [date_time "droplet_R0_0.2" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end

#%%
################################
# Spinodal decomposition example
################################
nx = 128
tol = 1e-6
dt_0 = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
dts = [0.25 * dt_0, 0.5 * dt_0, dt_0, 2 * dt_0, 4 * dt_0]
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_same_total_time"
total_time = 0.02
phi = initialization_from_file("$outdir/phi_128_initial.txt", nx, nx, delim=' ')
for dt in dts
    println(dt)
    max_it = Int.(total_time / dt)
    date_time = now()
    time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$dt", check_dir=false)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, [date_time "droplet_R0_0.2" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end

#%%
################################
# Spinodal decomposition example - every output
################################
nx = 128
tol = 1e-6
dt_0 = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
dts = [0.25 * dt_0, dt_0]
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal__same_total_time"
total_time = 0.002
phi = initialization_from_file("$outdir/phi_128_initial.txt", nx, nx, delim=' ')
for dt in dts
    println(dt)
    max_it = Int.(total_time / dt)
    date_time = now()
    time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_everytimestep", check_dir=false, ns=1)
    # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    # writedlm(f, [date_time "droplet_R0_0.2" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    # end
end

#%%
################################
# Spinodal decomposition example - from smoothed (multigrid) IC
################################
nx = 128
tol = 1e-6
dt_0 = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
dts = [0.25 * dt_0, dt_0]
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/smoothed_with_multigrid/"
total_time = 0.02
initial_time = 1
for initial_time in [1, 5, 500]
    phi = initialization_from_file("$outdir/$(initial_time)_dt_6p25e-6.txt", nx, nx, delim=' ')
    for dt in dts
        println(dt)
        max_it = Int.(total_time / dt)
        date_time = now()
        time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$(dt)_initial_$(initial_time)", check_dir=false)
        # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        # writedlm(f, [date_time "droplet_R0_0.2" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
        # end
    end
end



#%%
################################
# alpha scan
################################
include("./CH_multigrid_solver_with_alpha.jl")

nx = 128
tol = 1e-6
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
max_it = 10

ns = 10
# alphas = [-0.4, -0.2, 0.0, 0.2, 0.4]
alphas = [0.4]
phi = initialization(nx, nx, method="droplet", h=1 / 128, gam=epsilon)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius_alpha"
for alpha in alphas
    time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, max_it=max_it, dt=dt, gam=epsilon, print_mass=true, print_e=true, overwrite=false, suffix="_alpha_$alpha", check_dir=false, alpha=alpha)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, ["critical radius alpha" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end

#%%
################################
# alpha scan v2
################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

nx = 128
tol = 1e-6
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
total_time = 0.25
max_it = Int.(round(total_time / dt))
ns = 10
alphas = [-0.4, -0.2, 0.0, 0.2, 0.4, -0.002, 0.002]
# alphas = [0.4]
phi = initialization(nx, nx, method="droplet", h=1 / 128, gam=epsilon)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius_alpha_v2"
for alpha in alphas
    time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, max_it=max_it, dt=dt, gam=epsilon, print_mass=true, print_e=true, overwrite=false, suffix="_alpha_$alpha", check_dir=false, alpha=alpha)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, ["critical radius alpha" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end
#%%
################################
# alpha scan v2
################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

nx = 256
tol = 1e-6
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
total_time = 0.25
max_it = Int.(round(total_time / dt))
ns = 10
alphas = [-0.4, -0.2, 0.0, 0.2, 0.4, -0.002, 0.002]
# alphas = [0.4]
phi = initialization(nx, nx, method="droplet", h=1 / 256, gam=epsilon)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius_alpha_v2"
for alpha in alphas
    time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, max_it=max_it, dt=dt, gam=epsilon, print_mass=true, print_e=true, overwrite=false, suffix="_alpha_$(alpha)_eps_$(epsilon)", check_dir=false, alpha=alpha)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, ["critical radius alpha" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end

#%%
################################
# alpha scan v2
################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

nx = 512
tol = 1e-6
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
total_time = 0.25
max_it = Int.(round(total_time / dt))
ns = 10
alphas = [-0.4, -0.2, 0.0] #0.2, 0.4, -0.002, 0.002 run on Rivanna
# alphas = [0.4]
phi = initialization(nx, nx, method="droplet", h=1 / 512, gam=epsilon)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius_alpha_v2"
for alpha in alphas
    time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, max_it=max_it, dt=dt, gam=epsilon, print_mass=true, print_e=true, overwrite=false, suffix="_alpha_$(alpha)_eps_$(epsilon)", check_dir=false, alpha=alpha)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, ["critical radius alpha" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end

#%%
################################
# CPC geometry 
################################
include("./CH_multigrid_solver.jl")

nx = 256
ny = 256
tol = 1e-5
dt = 0.1 * (1 / 256)^2
m = 8
epsilon = m * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
# total_time = 0.25
# max_it = Int.(round(total_time / dt))
max_it = 19660
ns = 10
# CPC_width = 20
# experiments = .173 um = .108 of width 
cohesin_width = 4
for CPC_width in [10]
    phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
    # phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
    outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
    time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, print_r=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)", check_dir=false)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, ["CPC_geometry" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end

#%%
################################
# CPC geometry: compare grid sizes
################################
include("./CH_multigrid_solver.jl")

nx = 128
ny = 128
tol = 1e-5
dt = 0.1 * (1 / 256)^2
m = 8
epsilon = m * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
# total_time = 0.03
# max_it = Int.(round(total_time / dt))
max_it = 19660
ns = 10
CPC_width = 5
# experiments = .173 um = .108 of width 
cohesin_width = 2
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, print_r=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
end

#%%
################################
# CPC geometry: compare grid sizes
################################
include("./CH_multigrid_solver.jl")

nx = 128
ny = 128
tol = 1e-5
dt = 0.1 * (1 / 256)^2
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
# total_time = 0.03
# max_it = Int.(round(total_time / dt))
max_it = 19660
ns = 10
CPC_width = 5
# experiments = .173 um = .108 of width 
cohesin_width = 2
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, print_r=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
end


#%%
################################
# CPC geometry: compare grid sizes
################################
include("./CH_multigrid_solver.jl")

nx = 512
ny = 512
tol = 1e-5
dt = 0.1 * (1 / 256)^2
m = 8
epsilon = m * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
# total_time = 0.03
# max_it = Int.(round(total_time / dt))
max_it = 19660
ns = 10
CPC_width = 20
# experiments = .173 um = .108 of width 
cohesin_width = 8
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry_grid_sizes"
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, print_r=true, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
end

#%%
################################
# CPC geometry: compare grid sizes
################################
include("./CH_multigrid_solver.jl")

nx = 1024
ny = 1024
tol = 1e-5
dt = 0.5 * 0.1 * (1 / 256)^2
m = 8
epsilon = m * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
# total_time = 0.03
# max_it = Int.(round(total_time / dt))
max_it = 19660 * 2 #twice at many iterations for smaller timestep
ns = 10
CPC_width = 40
# experiments = .173 um = .108 of width 
cohesin_width = 16
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, print_r=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
end


#%%
################################
# CPC geometry : compare to SAV
################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
m = 8
epsilon = m * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
total_time = 0.025
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
factor = Int.(nx / 128)
cohesin_width = 2 * factor
CPC_width = 5 * factor
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
alpha = -0.5
time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_alpha_-0.5" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
end

#%%
################################
# CPC geometry : compare to SAV
################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
total_time = 0.025
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
factor = Int.(nx / 128)
cohesin_width = 2 * factor
CPC_width = 5 * factor
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
alpha = -0.5
time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_alpha_-0.5" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
end

#%%
###############################################
# CPC geometry : alpha = -0.5, epsilon scan
###############################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
# epsilon = 0.14
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
factor = Int.(nx / 128)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"

for epsilon in [0.08, 0.14]
    for c in [2, 6, 8, 10]
        for CPC in [5, 10, 12, 14, 16]
            println("c=$(c), CPC=$(CPC), epsilon=$(epsilon)")
            if (c == 8) & (CPC == 14)
                println("Already completed")
            elseif (c == 2) & (CPC == 5) & (epsilon == 0.08)
                println("Already completed")
            else
                cohesin_width = c * factor #changed from 2
                CPC_width = CPC * factor #changed from 5
                phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
                # phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
                alpha = -0.5
                time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
                open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                    writedlm(f, ["CPC_geometry_alpha_-0.5" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                end
            end
        end
    end
end

#%%
###############################################
# CPC geometry : alpha = -0.5, epsilon scan
###############################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
# epsilon = 0.14
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
factor = Int.(nx / 128)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"

for epsilon in [0.14]
    for c in [2, 6, 8] #number of grid points in 128 case
        for CPC in [18, 20, 22, 24, 26]
            println("c=$(c), CPC=$(CPC), epsilon=$(epsilon)")
            if (c == 8) & (CPC == 14)
                println("Already completed")
            elseif (c == 2) & (CPC == 5) & (epsilon == 0.08)
                println("Already completed")
            else
                cohesin_width = c * factor #changed from 2
                CPC_width = CPC * factor #changed from 5
                phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
                # phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
                alpha = -0.5
                time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
                open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                    writedlm(f, ["CPC_geometry_alpha_-0.5" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                end
            end
        end
    end
end

#%%
###############################################
# CPC geometry : alpha = -0.5, epsilon scan
###############################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
# epsilon = 0.14
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
factor = Int.(nx / 128)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"

for epsilon in [0.01, 0.02]
    for c in [2, 6, 10] #number of grid points in 128 case
        for CPC in [5, 10, 14, 16]
            println("c=$(c), CPC=$(CPC), epsilon=$(epsilon)")
            # if (c == 8) & (CPC == 14)
            #     println("Already completed")
            # elseif (c == 2) & (CPC == 5) & (epsilon == 0.08)
            #     println("Already completed")
            # else
            cohesin_width = c * factor #changed from 2
            CPC_width = CPC * factor #changed from 5
            phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
            # phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
            alpha = 0
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end

#%%
###############################################
# CPC geometry : alpha = -0.5, epsilon scan, fixed IC
###############################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"

for epsilon in [0.04]
    for cohesin_width in [0.1, 0.2, 0.3] #in um: radius of CPC droplet, experimental = 0.2
        for CPC_width in [0.12, 0.173, 0.22, 0.25, 0.3] #in um: total width of cohesin stripe, experimental = 0.173
            println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
            phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=3.2)
            alpha = -0.5
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end
#%%
###############################################
# CPC geometry : alpha = -0.5, epsilon scan, fixed IC
###############################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.2"

for epsilon in [0.0125]
    for cohesin_width in [0.1, 0.2, 0.3] #in um: radius of CPC droplet, experimental = 0.2
        for CPC_width in [0.12, 0.173, 0.22, 0.25, 0.3] #in um: total width of cohesin stripe, experimental = 0.173
            println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
            phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=3.2)
            alpha = -0.2
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end

#%%
###############################################
# CPC geometry : alpha = -0.5, epsilon scan, fixed IC
###############################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"

for epsilon in [0.0031]
    for cohesin_width in [0.1, 0.2, 0.3] #in um: radius of CPC droplet, experimental = 0.2
        for CPC_width in [0.12, 0.173, 0.22, 0.25, 0.3] #in um: total width of cohesin stripe, experimental = 0.173
            println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
            phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=3.2)
            alpha = 0
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end

#%%
###############################################
# CPC geometry : change domain to [0,2]x[0,2]
# length scale: 1 grid point = 25 nm (not 12.5)
# domain width = 6.4um
# epsilon is halved from normal length scale
###############################################
include("./CH_multigrid_solver_with_alpha_change_domain.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"

for epsilon in [0.00375, 0.0075, 0.015009]
    for cohesin_width in [0.1, 0.2] #in um: radius of CPC droplet
        for CPC_width in [0.125, 0.173, 0.22] #in um: total width of cohesin stripe, experimental = 0.173
            println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
            phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=6.4)
            alpha = 0
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)_domain_0_2", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end

################################################################################################################################################################
# VCell IC

#%%
###############################################
# CPC geometry : VCell IC
###############################################
using Glob
datadir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/post_transition/256"
ic_files = glob("*8.4max.csv", datadir)
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx

for file in ic_files
    # if file == "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/03_25_24_CPC_relaxed_RefModel_128x64_03_25_24_relaxed_RefModel_100_256x256_100s_8.4max.csv"
    #     println("already run")
    # else
    for epsilon in [0.00375, 0.006, 0.0075, 0.015009, 0.030019]
        for alpha in [-0.5, -0.2, 0, 0.2] #in um: radius of CPC droplet
            phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
            s = split(last(split(file, "/")), ".csv")[1]
            outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
            outdir = "$(outdir_prefix)/$(s)"
            println(outdir)
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end
# end

#%%
###############################################
# CPC geometry : VCell IC
###############################################
using Glob
datadir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/pre_transition_relaxed/"
ic_files = glob("*8.4max.csv", datadir)
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx

for file in ic_files
    # if file == "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/03_25_24_CPC_relaxed_RefModel_128x64_03_25_24_relaxed_RefModel_100_256x256_100s_8.4max.csv"
    #     println("already run")
    # else
    for epsilon in [0.00375, 0.006, 0.0075, 0.015009, 0.030019]
        for alpha in [-0.5, -0.2, 0, 0.2] #in um: radius of CPC droplet
            phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
            s = split(last(split(file, "/")), ".csv")[1]
            outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
            outdir = "$(outdir_prefix)/$(s)"
            println(outdir)
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end
# end

#%%
###############################################
# CPC geometry : VCell IC after 70s relaxed
###############################################
using Glob
datadir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/post_transition_70s/"
ic_files = glob("*8.4max.csv", datadir)
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx

for file in ic_files
    # if file == "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/03_25_24_CPC_relaxed_RefModel_128x64_03_25_24_relaxed_RefModel_100_256x256_100s_8.4max.csv"
    #     println("already run")
    # else
    for epsilon in [0.00375, 0.006, 0.0075, 0.0125, 0.015009, 0.030019, 0.04]
        for alpha in [-0.5, -0.2, 0, 0.2] #in um: radius of CPC droplet
            phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
            s = split(last(split(file, "/")), ".csv")[1]
            outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
            outdir = "$(outdir_prefix)/$(s)"
            println(outdir)
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end
# end

#%%

include("./Cahn_Hilliard_solvers/julia_multigrid/CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx

file = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/pre_transition_relaxed/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_70_256x256_70s_8.4max.csv"
for epsilon in [0.00375, 0.006, 0.0075, 0.0125, 0.015009, 0.030019, 0.04]
    for alpha in [-0.5, -0.2, 0, 0.2] #in um: radius of CPC droplet
        phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
        s = split(last(split(file, "/")), ".csv")[1]
        outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
        outdir = "$(outdir_prefix)/$(s)"
        println(outdir)
        time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
        open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
            writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
            # end
        end
    end
end

#%%

include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx

file = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/post_transition_70s/10_24_23_CPC_tensed_RefModel_128x64_post_transition_07_14_24_500s_post_transition_base_0_256x256_0s_8.4max.csv"
for epsilon in [0.00375, 0.006, 0.0075, 0.0125, 0.015009, 0.030019, 0.04]
    for alpha in [-0.5, -0.2, 0, 0.2] #in um: radius of CPC droplet
        phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
        s = split(last(split(file, "/")), ".csv")[1]
        outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
        outdir = "$(outdir_prefix)/$(s)"
        println(outdir)
        time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
        open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
            writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
            # end
        end
    end
end

#%%

include("./Cahn_Hilliard_solvers/julia_multigrid/CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx

file = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/pre_transition_relaxed/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_transactiv_70_256x256_70s_8.4max.csv"
for epsilon in [0.00375, 0.006, 0.0075, 0.0125, 0.015009, 0.030019, 0.04]
    for alpha in [-0.5, -0.2, 0, 0.2] #in um: radius of CPC droplet
        phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
        s = split(last(split(file, "/")), ".csv")[1]
        outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
        outdir = "$(outdir_prefix)/$(s)"
        println(outdir)
        time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
        open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
            writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
            # end
        end
    end
end

#%%

include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/"
files = [#"$(indir)/pre_transition_relaxed/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_70_256x256_70s_8.4max.csv",
    "$(indir)/post_transition_70s/10_24_23_CPC_tensed_RefModel_128x64_post_transition_07_14_24_500s_post_transition_base_20Pac_0_256x256_0s_8.4max.csv"]
for file in files
    for epsilon in [0.016, 0.018, 0.02, 0.025]
        alpha = -0.2
        phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
        s = split(last(split(file, "/")), ".csv")[1]
        outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
        outdir = "$(outdir_prefix)/$(s)"
        println(outdir)
        time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
        open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
            writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
            # end
        end
    end
end

#%%

include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/"
files = [#"$(indir)/pre_transition_relaxed/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_70_256x256_70s_8.4max.csv",
    "$(indir)/post_transition_70s/10_24_23_CPC_tensed_RefModel_128x64_post_transition_07_14_24_500s_post_transition_base_20Pac_0_256x256_0s_8.4max.csv"]
for file in files
    for epsilon in [0.024, 0.023, 0.022, 0.021]
        alpha = -0.2
        phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
        s = split(last(split(file, "/")), ".csv")[1]
        outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
        outdir = "$(outdir_prefix)/$(s)"
        println(outdir)
        time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
        open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
            writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
            # end
        end
    end
end


#%%
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.5
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 256
ny = nx
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/"
# file = "$(indir)/pre_transition_relaxed/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_100_256x256_100s_8.4max.csv"
epsilon = 0.0096
alpha = 0
# phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
# s = split(last(split(file, "/")), ".csv")[1]
# outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
# outdir = "$(outdir_prefix)/$(s)"
# println(outdir)
# time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
# open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
#     writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
#     # end
# end

file = "$(indir)/post_transition_100s/256/10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_0_256x256_0s_8.4max.csv"
phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
s = split(last(split(file, "/")), ".csv")[1]
outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
outdir = "$(outdir_prefix)/$(s)"
println(outdir)
time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    # end
end

#%%
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-6
max_it_CH = 10000
total_time = 0.02
max_it = Int.(round(total_time / dt))
ns = 10
nx = 256
ny = nx
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/"
# file = "$(indir)/pre_transition_relaxed/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_100_256x256_100s_8.4max.csv"
epsilon = 0.030019
alpha = -0.2
# phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
# s = split(last(split(file, "/")), ".csv")[1]
# outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
# outdir = "$(outdir_prefix)/$(s)"
# println(outdir)
# time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
# open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
#     writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
#     # end
# end

file = "$(indir)/post_transition_100s/256/10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_0_256x256_0s_8.4max.csv"
phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
s = split(last(split(file, "/")), ".csv")[1]
outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
outdir = "$(outdir_prefix)/$(s)"
println(outdir)
time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    # end
end


#%%
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
dt = 2.5e-6
max_it_CH = 10000
total_time = 0.02
max_it = Int.(round(total_time / dt))
ns = 10
nx = 256
ny = nx
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/data/"
file = "$(indir)/pre_transition_relaxed/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_100_256x256_100s_8.4max.csv"
epsilon = 0.030019
alpha = -0.2
phi = initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
s = split(last(split(file, "/")), ".csv")[1]
outdir_prefix = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC"
outdir = "$(outdir_prefix)/$(s)"
println(outdir)
time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_VCell 70s" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    # end
end

################################################################################################################################################################
# Back to CPC geometric IC

#%%
################################
# CPC geometry: compare grid sizes
################################
include("./CH_multigrid_solver_with_alpha_v2.jl")

tol = 1e-5
# dt = 0.0000125
dt = 0.1 * (1 / 256)^2

max_it_CH = 10000
total_time = 0.03
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 512
ny = nx
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"

for epsilon in [0.0096]
    for cohesin_width in [0.05, 0.15, 0.2] #in um: radius of CPC droplet, experimental = 0.2
        for CPC_width in [0.125, 0.1, 0.12, 0.15, 0.173, 0.22] #in um: total width of cohesin stripe, experimental = 0.173
            println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
            phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=3.2)
            alpha = 0
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end

#after this, still need to do .173/.1 and .22/.1 combos

#%%
###############################################
# CPC geometry : change domain to [0,2]x[0,2]
# length scale: 1 grid point = 25 nm (not 12.5)
# domain width = 6.4um
# epsilon should be doubled from normal length scale
# grid size is doubled to keep the same resolution (512)
###############################################
include("./Cahn_Hilliard_solvers/julia_multigrid/CH_multigrid_solver_with_alpha_change_domain.jl")

tol = 1e-5
dt = 2.5e-5
max_it_CH = 10000
total_time = 0.05
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
nx = 512
ny = nx
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"

for epsilon in [0.0192, 0.015]
    for cohesin_width in [0.1, 0.2] #in um: radius of CPC droplet
        for CPC_width in [0.125, 0.173, 0.22] #in um: total width of cohesin stripe, experimental = 0.173
            println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
            phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=6.4)
            alpha = 0
            time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)_domain_0_2", check_dir=false, alpha=alpha)
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
                writedlm(f, ["CPC_geometry_alpha_0" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
                # end
            end
        end
    end
end
