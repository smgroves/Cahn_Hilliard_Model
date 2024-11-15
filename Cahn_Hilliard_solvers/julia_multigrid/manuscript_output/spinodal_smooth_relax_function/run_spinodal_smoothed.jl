#%%

function relax(c_new, mu_new, su, sw, nxt, nyt, c_relax, xright, xleft, dt, Cahn)
    ht2 = ((xright - xleft) / nxt)^2
    # a = MVector{4,Float64}(undef)
    # f = MVector{2,Float64}(undef)
    a = zeros(Float64, 4)
    f = zeros(Float64, 2)
    for iter in 1:c_relax
        for i in 1:nxt
            for j in 1:nyt
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
                a[1] = 1 / dt
                a[2] = (x_fac + y_fac) / ht2
                #a[2] = -(x_fac + y_fac) * Cahn / ht2 - d2f(c_new[i][j]);
                a[3] = -(x_fac + y_fac) * Cahn / ht2 - 3 * (c_new[i, j])^2
                cnew_val = (c_new[i, j])
                d2f = -3 * (c_new[i, j])^2

                a[4] = 1.0

                f[1] = su[i, j]
                f[2] = sw[i, j] - 2 * (c_new[i, j])^3

                if i > 1
                    f[1] += mu_new[i-1, j] / ht2
                    f[2] -= Cahn * c_new[i-1, j] / ht2
                end
                if i < nxt
                    f[1] += mu_new[i+1, j] / ht2
                    f[2] -= Cahn * c_new[i+1, j] / ht2
                end
                if j > 1
                    f[1] += mu_new[i, j-1] / ht2
                    f[2] -= Cahn * c_new[i, j-1] / ht2
                end
                if j < nyt
                    f[1] += mu_new[i, j+1] / ht2
                    f[2] -= Cahn * c_new[i, j+1] / ht2
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

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_+1_-1_IC/output/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/"

include("../../CH_multigrid_solver.jl")
# tol = 1e-6
nx = 512
ny = nx
tol = 1e-6
dt = 6.25e-6
m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
total_time = 0.06
max_it = Int.(total_time / dt)
date_time = now()
phi = initialization_from_file("$(indir)/initial_phi_$(nx).csv", nx, nx)

sc, smu = source(phi, nx, ny, dt)

mu = zeros(Float64, nx, ny)
Cahn = epsilon^2  # ϵ^2 defined as a global variable

domain_left = 0
domain_right = 1
n_relax = 8
phi_smooth, mu_smooth = relax(phi, mu, sc, smu, nx, nx, n_relax, domain_right, domain_left, dt, Cahn)
writedlm("$(outdir)initial_phi_$(nx)_smooth_n_relax_$(n_relax).csv", phi, ',')

#%%
using DelimitedFiles
using Dates
using Plots

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output/"

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
