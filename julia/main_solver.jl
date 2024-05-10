include("solver_type_stable_v4.jl")

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/obsolete_tests"
suffix = "random_ic"
@time main(100, 10000, 1e-5, outdir, suffix)

