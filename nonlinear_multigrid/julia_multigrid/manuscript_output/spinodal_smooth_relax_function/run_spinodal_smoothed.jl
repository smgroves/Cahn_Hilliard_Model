#%% FIGURE 2 GENERATING DIFFERENT IC FOR SPINODAL DECOMPOSITION
import Dates
function relax(c_new, mu_new, su, sw, nxt, nyt, c_relax, xright, xleft, yright, yleft, dt, epsilon2, boundary)
    ht2 = ((xright - xleft) / nxt)^2
    a = MVector{4,Float64}(undef)
    f = MVector{2,Float64}(undef)
    for iter in 1:c_relax
        for i in 1:nxt
            for j in 1:nyt
                if boundary == "neumann"
                    if i > 1 && i < nxt
                        x_fac = 2.0
                    else
                        x_fac = 1.0
                    end
                    if j > 1 && j < nyt
                        y_fac = 2.0
                    else
                        y_fac = 1.0
                    end
                elseif boundary == "periodic"
                    x_fac = 2.0
                    y_fac = 2.0
                end
                a[1] = 1 / dt
                a[2] = (x_fac + y_fac) / ht2
                a[3] = -(x_fac + y_fac) * epsilon2 / ht2 - 3 * (c_new[i, j])^2
                a[4] = 1.0

                f[1] = su[i, j]
                f[2] = sw[i, j] - 2 * (c_new[i, j])^3

                if i > 1
                    f[1] += mu_new[i-1, j] / ht2
                    f[2] -= epsilon2 * c_new[i-1, j] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[nxt, j] / ht2
                    f[2] -= epsilon2 * c_new[nxt, j] / ht2
                end
                if i < nxt
                    f[1] += mu_new[i+1, j] / ht2
                    f[2] -= epsilon2 * c_new[i+1, j] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[1, j] / ht2
                    f[2] -= epsilon2 * c_new[1, j] / ht2
                end
                if j > 1
                    f[1] += mu_new[i, j-1] / ht2
                    f[2] -= epsilon2 * c_new[i, j-1] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[i, nyt] / ht2
                    f[2] -= epsilon2 * c_new[i, nyt] / ht2
                end
                if j < nyt
                    f[1] += mu_new[i, j+1] / ht2
                    f[2] -= epsilon2 * c_new[i, j+1] / ht2
                elseif boundary == "periodic"
                    f[1] += mu_new[i, 1] / ht2
                    f[2] -= epsilon2 * c_new[i, 1] / ht2
                end
                det = a[1] * a[4] - a[2] * a[3]
                c_new[i, j] = (a[4] * f[1] - a[2] * f[2]) / det
                mu_new[i, j] = (-a[3] * f[1] + a[1] * f[2]) / det

            end
        end
    end
    return c_new, mu_new
end

function source(c_old, nx, ny, dt)
    src_mu = zeros(Float64, nx, ny)
    src_c = zeros(Float64, nx, ny)
    ct = laplace(c_old, nx, ny)
    for i in 1:nx
        for j in 1:ny
            src_c[i, j] = c_old[i, j] / dt - ct[i, j]
            src_mu[i, j] = 0
        end
    end
    return src_c, src_mu
end

#%%
using DelimitedFiles
using Dates

# indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_+1_-1_IC/output/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC_FIGURE1_FIGURE2/"

include("../../CH_multigrid_solver.jl")
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
phi = biased_spin_matrix(nx, p_plus)
writedlm("$(outdir)$(perc)/initial_phi_$(nx)_$(perc).csv", phi, ',')
#%%

boundary = "periodic"


sc, smu = source(phi, nx, ny, dt)

mu = zeros(Float64, nx, ny)
Cahn = epsilon^2  # ϵ^2 defined as a global variable

domain_left = 0
domain_right = 1
n_relax = 16
tol = 1e-6
dt = 6.25e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
total_time = 0.06
max_it = Int.(total_time / dt)
phi_smooth, mu_smooth = relax(phi, mu, sc, smu, nx, nx, n_relax, domain_right, domain_left, domain_right, domain_left, dt, Cahn, boundary)
writedlm("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary).csv", phi, ',')


#%% visualize the initial condition from file
# import Pkg;
# Pkg.add("Plots");
using Plots
using DelimitedFiles
perc = "50p"
nx = 128
boundary = "periodic"
n_relax = 4
phi_smooth = initialization_from_file("$(outdir)$(perc)/initial_phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary).csv", nx, nx)
gr()
heatmap(1:size(phi_smooth, 1),
    1:size(phi_smooth, 2), phi_smooth,
    c=cgrad([:blue, :white, :red]),
    title="Initial condition for phi_$(nx)_smooth_n_relax_$(n_relax)_$(perc)_$(boundary)",
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

include("../../CH_multigrid_solver.jl")
nx = 128
dt = 5.5e-6
n_relax = 4
ny = nx
tol = 1e-5
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
# total_time = (0.1 / 64^2) * 500

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
