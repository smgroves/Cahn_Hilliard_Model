#%%

using DelimitedFiles
using Dates

nx = 128
phi = fill(-1.0, nx, nx)
width = 20
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/stripe_IC/IC/"
for i in 1:nx
    for j in 1:nx
        if abs(j - nx / 2) < width / 2
            phi[i, j] = 1
        end
    end
end
writedlm("$(outdir)initial_phi_$(nx)_width_$(width).csv", phi, ',')

#%%
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/stripe_IC/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/stripe_IC/output/"
width = 20
nx = 128
include("../../CH_multigrid_solver.jl")

for dt = [1e-1 / (nx^2), 1 / (nx^2)]
    ny = nx
    tol = 1e-6
    m = 8
    epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
    # total_time = 0.01
    total_time = (1 / 64^2) * 500
    max_it = Int.(total_time / dt)
    date_time = now()
    name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_eps_$(epsilon)"
    phi = initialization_from_file("$(indir)initial_phi_$(nx)_width_$(width).csv", nx, nx)
    final_phi = multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=false, print_phi=false, suffix=name, check_dir=false)
    writedlm("$(outdir)/$(name)_final_phi.csv", final_phi, ',')
end