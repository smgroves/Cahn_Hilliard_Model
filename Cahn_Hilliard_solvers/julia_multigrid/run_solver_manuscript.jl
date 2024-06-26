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
include("./CH_multigrid_solver.jl")

tol = 1e-5
dt = 2.5e-5
m = 8
epsilon = m * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
total_time = 0.025
max_it = Int.(round(total_time / dt))
# max_it = 19660
ns = 10
for nx in [128, 256, 512]
    ny = nx
    factor = Int.(nx / 128)
    cohesin_width = 2 * factor
    CPC_width = 5 * factor
    phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
    # phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
    outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
    time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, print_r=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)", check_dir=false)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, ["CPC_geometry" "Julia" nx epsilon dt tol max_it max_it_CH time_passed], ",")
    end
end