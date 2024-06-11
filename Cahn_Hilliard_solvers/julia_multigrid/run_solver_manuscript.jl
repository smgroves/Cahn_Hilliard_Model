#%%
# using Pkg
# Pkg.add("DataFrames")
# Pkg.add("BenchmarkTools")
# Pkg.add("StaticArrays")
# Pkg.add("ProfileView")
# Pkg.add("DelimitedFiles")
# Pkg.add("LinearAlgebra")
# Pkg.add("Printf")
include("./CH_multigrid_solver.jl")


#########################
# Single droplet example
#########################
#%%
nx = 128
tol = 1e-6
dt_0 = 2.5e-5
M = 8
gam = M * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
dts = [0.25 * dt_0, 0.5 * dt_0, dt_0, 2 * dt_0, 4 * dt_0]
total_time = 1
max_it = Int.(total_time / dt)
phi = initialization(nx, nx, method="droplet", h=1 / 128, gam=gam)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/droplet_2"
for dt in dts
    println(dt)
    main(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="_dt_$dt", check_dir=false)
end

#%%
################################
# Spinodal decomposition example
################################
nx = 128
tol = 1e-6
dt_0 = 2.5e-5
M = 8
gam = M * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
dts = [0.25 * dt_0, 0.5 * dt_0, dt_0, 2 * dt_0, 4 * dt_0]
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_2"
total_time = 0.02
max_it = Int.(total_time / dt)
phi = initialization_from_file("$outdir/phi_128_initial.txt", nx, nx, delim=' ')
for dt in dts
    println(dt)
    main(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, suffix="dt_$dt", check_dir=false)
end


###################################
# Critical radius - epsilon scan
###################################

### COPY FROM RIVANNA WHEN DONE

#%%
################################
# alpha scan
################################
include("./CH_multigrid_solver_with_alpha.jl")

nx = 128
tol = 1e-6
dt = 2.5e-5
M = 8
gam = M * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
max_it = 10000

ns = 10
# alphas = [-0.4, -0.2, 0.0, 0.2, 0.4]
alphas = [0.0]
phi = initialization(nx, nx, method="droplet", h=1 / 128)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius"
for alpha in alphas
    time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=true, print_e=true, overwrite=false, suffix="_alpha_$alpha", check_dir=false, alpha=alpha)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
        writedlm(f, ["critical radius" "Julia" nx gam dt tol max_it max_it_CH time_passed], ",")
    end
end