#%%

using DelimitedFiles
using Dates

include("../../CH_multigrid_solver.jl")
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/flat_interface/IC"

function initialization_tanh_transition(nx, ny, h; gam=0.1, interface_location=nx / 2)
    phi = zeros(Float64, nx, ny)            # Initialize the phi array to zeros
    x = h .* (0:nx-1)                       # x-coordinates: scaled by step size h
    y = h .* (0:ny-1)                       # y-coordinates: scaled by step size h
    xx, yy = meshgrid(x, y)                 # Create a meshgrid of x and y

    # Create a transition along the x-axis
    phi = @.tanh((xx .- interface_location * h) / (sqrt(2) * gam))

    return phi
end

m = 8
epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
nx = 128
phi = initialization_tanh_transition(nx, nx, 1 / 128, gam=epsilon)
writedlm("$(indir)/initial_phi_$(nx).csv", phi, ',')

#%%
total_time = 0.01
dt = 6.25e-6
tol = 1e-6
max_it = Int.(total_time / dt)
date_time = now()
name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_eps_$(epsilon)"
phi = initialization_from_file("$(indir)/initial_phi_128.csv", nx, nx)
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/flat_interface/output"
time_passed = @elapsed multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=true, print_phi=true, suffix=name, check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time "spinodal_t_0.01_IC" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end

