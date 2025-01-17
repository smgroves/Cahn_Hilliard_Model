#%%
using Plots
using DelimitedFiles
using Dates
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/IC/"

function create_matrix_with_circles(radius_1, radius_2, space, nx)
    # Initialize the matrix with -1
    matrix = -ones(Int, nx, nx)

    # Define the center row (central horizontal axis)
    center_row = nx / 2

    # Calculate the centers of the two circles
    # center1_col = nx / 2 - div(space, 2) - radius_1
    # center2_col = nx / 2 + div(space, 2) + radius_2
    dist_from_edge = floor(0.5 * (nx - 2 * (radius_1) - space - 2 * (radius_2)))
    # Position the left circle 2 grid points from the left edge
    center1_col = radius_1 + dist_from_edge

    # Calculate the center of the second circle based on the space
    center2_col = center1_col + radius_1 + space + radius_2

    # Draw the first circle
    for i in 1:nx
        for j in 1:nx
            if (i - center_row)^2 + (j - center1_col)^2 ≤ radius_1^2
                matrix[i, j] = 1
            end
        end
    end

    # Draw the second circle
    for i in 1:nx
        for j in 1:nx
            if (i - center_row)^2 + (j - center2_col)^2 ≤ radius_2^2
                matrix[i, j] = 1
            end
        end
    end

    return matrix
end

nx = 128
radius_1 = 10
radius_2 = 15
space = 10
phi = create_matrix_with_circles(radius_1, radius_2, space, nx)
writedlm("$(outdir)initial_phi_$(nx)_r1_$(radius_1)_r2_$(radius_2)_space_$(space).csv", phi, ',')

#%%
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/output/"
nx = 128
include("../../CH_multigrid_solver.jl")

for dt = [1e-4, 1e-3]
    ny = nx
    tol = 1e-6
    m = 8
    epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
    total_time = 0.2
    # total_time = (1 / 64^2) * 500
    max_it = Int.(total_time / dt)
    date_time = now()
    name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_eps_$(epsilon)_r1_$(radius_1)_r2_$(radius_2)_space_$(space)"
    phi = initialization_from_file("$(indir)initial_phi_$(nx)_r1_$(radius_1)_r2_$(radius_2)_space_$(space).csv", nx, nx)
    final_phi = multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=false, print_phi=true, suffix=name, check_dir=false)
    writedlm("$(outdir)/$(name)_final_phi.csv", final_phi, ',')
end

#%%
using Plots
using DelimitedFiles
using Dates
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/IC/"
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/output/"
nx = 128
include("../../CH_multigrid_solver.jl")
radius_1 = 20
radius_2 = 30
space = 10
for dt = [1e-3, 1e-4]
    ny = nx
    tol = 1e-6
    m = 8
    epsilon = m * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
    total_time = 5
    # total_time = (1 / 64^2) * 500
    max_it = Int.(total_time / dt)
    date_time = now()
    name = "MG_$(max_it)_dt_$(dt)_Nx_$(nx)_eps_$(epsilon)_r1_$(radius_1)_r2_$(radius_2)_space_$(space)"
    phi = initialization_from_file("$(indir)initial_phi_$(nx)_r1_$(radius_1)_r2_$(radius_2)_space_$(space).csv", nx, nx)
    final_phi = multigrid_solver(phi, nx, tol, outdir, dt=dt, epsilon=epsilon, max_it=max_it, print_mass=true, print_e=true, overwrite=false, print_r=false, print_phi=false, suffix=name, check_dir=false)
    writedlm("$(outdir)/$(name)_final_phi.csv", final_phi, ',')
end