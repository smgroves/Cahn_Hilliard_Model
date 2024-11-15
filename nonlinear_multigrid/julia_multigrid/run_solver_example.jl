# uncomment the following line to install the packages if needed
# using Pkg
# Pkg.add("DataFrames")
# Pkg.add("BenchmarkTools")
# Pkg.add("StaticArrays")
# Pkg.add("ProfileView")
# Pkg.add("DelimitedFiles")
# Pkg.add("LinearAlgebra")
# Pkg.add("Printf")

include("./CH_multigrid_solver.jl")
nx = 128
tol = 1e-6

#####################################
# EXAMPLE 1: SPINODAL DECOMPOSITION
#####################################

#initialize phi as a 128x128 grid of -1 and +1 (spinodal decomposition, default)
#other options are: "random" (-1 to 1 range), "geometric" (CPC initial conditions from paper), "droplet" (single droplet from paper), "file" (import from file)
phi = initialization(nx, nx)
#run the CH solver with tol = 1e-5 for 10,000 time steps. Note option to calculate mass and energy at each time step.
#note that the default is to print out every 10 time steps to a csv for phi/psi.
outdir = "./example"
main(phi, nx, tol, outdir, max_it=10000, print_mass=true, print_e=true)

#####################################
# EXAMPLE 2: INITIALIZATION FROM FILE
#####################################
# note that if your imported file is concentrations (>=0), they will be rescaled to the range 0 to 1 for the solver, 
# making the assumption that the maximum concentration in the data has phase = 1. The output will be psi in [0,1] rather than phi [-1,1].
outdir = "./example_from_text"
psi = initialization(nx, nx, method="file", initial_file="./input/example.csv", delim=",")
main(psi, nx, tol, outdir, max_it=10000)

