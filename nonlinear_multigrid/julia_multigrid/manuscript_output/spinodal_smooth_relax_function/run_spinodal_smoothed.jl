#%% FIGURE 2 GENERATING DIFFERENT IC FOR SPINODAL DECOMPOSITION

#%%
using DelimitedFiles
using Dates

#use the official final code for smoothing
include("/Users/smgroves/Documents/GitHub/CHsolvers_package/CahnHilliard_Julia_solvers/nmg_solver.jl")
include("/Users/smgroves/Documents/GitHub/CHsolvers_package/CahnHilliard_Julia_solvers/relax.jl")

# indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_+1_-1_IC/output/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC_FIGURE1_FIGURE2/"

#%%
# include("../../CH_multigrid_solver.jl")
#%%
# tol = 1e-6
nx = 128
ny = nx
perc = "50p"

date_time = now()

using Random
Random.seed!(1234)

function biased_spin_matrix(nx, p_plus=0.75)
    total = nx * nx
    n_plus = round(Int, p_plus * total)
    n_minus = total - n_plus

    # Create the array
    values = vcat(fill(1.0, n_plus), fill(-1.0, n_minus))
    shuffle!(values)

    # Reshape into matrix
    return reshape(values, nx, nx)
end

if perc == "25p"
    p_plus = 0.25
elseif perc == "50p"
    p_plus = 0.5
elseif perc == "75p"
    p_plus = 0.75
end
# phi = biased_spin_matrix(nx, p_plus)
# writedlm("$(outdir)$(perc)/initial_phi_$(nx)_$(perc).csv", phi, ',')
#%%
n_relax = 4
boundary = "periodic"

domain_left = 0
domain_right = 1
tol = 1e-6
dt = 6.25e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
Cahn = epsilon^2  # ϵ^2 defined as a global variable

phi = initialization_from_file("$(outdir)$(perc)/initial_phi_$(nx)_$(perc).csv", nx, nx)
sc, smu = source(phi, nx, ny, dt, domain_right, domain_left, domain_right, domain_left, boundary)
mu = zeros(Float64, nx, ny)
relax!(phi, mu, sc, smu, nx, nx, n_relax, domain_right, domain_left, domain_right, domain_left, dt, Cahn, boundary)

# phi_smooth, mu_smooth = relax(phi, mu, sc, smu, nx, nx, n_relax, domain_right, domain_left, domain_right, domain_left, dt, Cahn, boundary)
writedlm("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary)_v2.csv", phi, ',')


#%% visualize the initial condition from file
# import Pkg;
# Pkg.add("Plots");
using Plots
using DelimitedFiles
perc = "50p"
nx = 128
boundary = "periodic"
n_relax = 16
phi_smooth = initialization_from_file("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary)_v2.csv", nx, nx)
gr()
heatmap(1:size(phi_smooth, 1),
    1:size(phi_smooth, 2), phi_smooth,
    c=cgrad([:blue, :white, :red]),
    title="Initial condition for phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary)",
    aspect_ratio=1,
    clim=(-1, 1),
    xlims=(1, nx),
    ylims=(1, nx))

#%% comparing v1 and v2 for different levels of smoothing
using Plots
using DelimitedFiles
perc = "50p"
nx = 128
boundary = "periodic"
n_relax = 4
phi_smooth_v1 = initialization_from_file("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary).csv", nx, nx)
phi_smooth_v2 = initialization_from_file("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary)_v2.csv", nx, nx)

gr()
heatmap(phi_smooth_v1 - phi_smooth_v2,
    c=cgrad([:blue, :white, :red]),
    title="Initial condition for phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary)",
    aspect_ratio=1,
    clim=(-1, 1),
    xlims=(1, nx),
    ylims=(1, nx))

#RESULT: n_relax = 4 and n_relax = 8 were correct. I just need to rerun n=2 and n = 16 to get L2 error.

#%% compare different levels of smoothing
perc = "50p"
nx = 128
boundary = "neumann"
smooth1 = 1
smooth2 = 16
phi_smooth_1 = initialization_from_file("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(smooth1)_$(perc)_$(boundary)_v2.csv", nx, nx)
phi_smooth_2 = initialization_from_file("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(smooth2)_$(perc)_$(boundary)_v2.csv", nx, nx)

gr()
heatmap(phi_smooth_1 - phi_smooth_2,
    c=cgrad([:blue, :white, :red]),
    title="n = $(smooth1) - n = $(smooth2)",
    aspect_ratio=1,
    clim=(-1, 1),
    xlims=(1, nx),
    ylims=(1, nx))

#%%
# build smaller grid-size ICs from 512 gridsize
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/"
using DelimitedFiles
using DataFrames
function restrict_ch(uf, vf, nxc, nyc)
    uc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    vc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    for i in 1:nxc
        for j in 1:nyc
            uc[i, j] = 0.25 * (uf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i - 1), round(Int, 2 * j)] + uf[round(Int, 2 * i), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i), round(Int, 2 * j)])
            vc[i, j] = 0.25 * (vf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i - 1), round(Int, 2 * j)] + vf[round(Int, 2 * i), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i), round(Int, 2 * j)])
        end
    end
    return uc, vc
end
function initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
    phi = readdlm(file, delim, Float64)
    if size(phi) != (nx, ny)
        print("Warning: phi from file is wrong size: $(size(phi)) Expected: $(nx), $(ny)")
    end
    if transpose_matrix
        phi = transpose(phi)
    end
    return phi
end

nx = 512
ny = nx
n_relax = 4
phi = initialization_from_file("$(indir)initial_phi_$(nx)_smooth_n_relax_$(n_relax)_25p.csv", nx, nx)
mu = zeros(Float64, nx, ny)

new_nx = 256
phi_small, mu_small = restrict_ch(phi, mu, new_nx, new_nx)
writedlm("$(indir)initial_phi_$(new_nx)_smooth_n_relax_$(n_relax)_from512_25p.csv", phi_small, ',')
new_nx = 128
phi_small, mu_small = restrict_ch(phi_small, mu_small, new_nx, new_nx)
writedlm("$(indir)initial_phi_$(new_nx)_smooth_n_relax_$(n_relax)_from512_25p.csv", phi_small, ',')
new_nx = 64
phi_small, mu_small = restrict_ch(phi_small, mu_small, new_nx, new_nx)
writedlm("$(indir)initial_phi_$(new_nx)_smooth_n_relax_$(n_relax)_from512_25p.csv", phi_small, ',')

#%%
using DelimitedFiles
using Dates
using Plots

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output/"

include("../../CH_multigrid_solver.jl")
for nx in [64, 128, 256]
    for dt in [1e-6, 1e-5, 1e-4, 1e-7]
        ny = nx
        n_relax = 8
        tol = 1e-6
        dt = 6.25e-6
        m = 8
        epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
        total_time = 0.01
        max_it = Int.(total_time / dt)
        date_time = now()
        phi = initialization_from_file("$(indir)initial_phi_$(nx)_smooth_n_relax_$(n_relax).csv", nx, nx)
        final_phi = multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=false, suffix="dt_$(dt)_smoothed_$(n_relax)", check_dir=false)
        gr()
        heatmap(1:size(final_phi, 1),
            1:size(final_phi, 2), final_phi,
            c=cgrad([:blue, :white, :red]),
            title="Final value t = $(total_time) \n for phi_$(nx)_smooth_n_relax_$(n_relax)_dt_$(dt)")
    end
end

#%%
using DelimitedFiles
using Dates
using Plots

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output/"

include("../../CH_multigrid_solver.jl")
for nx in [128]
    for dt = [1e-2 / (nx^2), 1e-1 / (nx^2)]
        for n_relax = [2, 4, 8]
            ny = nx
            tol = 1e-6
            m = 8
            epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
            # total_time = 0.01
            total_time = (0.1 / 64^2) * 500
            max_it = Int.(total_time / dt)
            date_time = now()
            name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_n_relax_$(n_relax)_eps_$(epsilon)"
            phi = initialization_from_file("$(indir)initial_phi_$(nx)_smooth_n_relax_$(n_relax).csv", nx, nx)
            final_phi = multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=false, print_phi=true, ns=1000, suffix=name, check_dir=false)
            writedlm("$(outdir)/$(name)_final_phi.csv", final_phi, ',')
        end
    end
end

#%%
using DelimitedFiles
using Dates
using Plots

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output/"

nx = 128
dt = 5.5e-6
n_relax = 4
ny = nx
tol = 1e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))

max_it = 2000
total_time = max_it * dt
date_time = now()
name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_n_relax_$(n_relax)_eps_$(epsilon)"
phi = initialization_from_file("$(indir)initial_phi_$(nx)_smooth_n_relax_$(n_relax).csv", nx, nx)
date_time = now()
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, ns=10, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, print_phi=false, suffix=name, check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_smoothed_no_print" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end
# writedlm("$(outdir)/$(name)_final_phi.csv", final_phi, ',')
