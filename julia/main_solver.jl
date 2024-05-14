################ solve v4 #####################
# version 4 includes global variables (const). 
using Plots
using DelimitedFiles

include("solver_type_stable_v4.jl")
# To execute code blocks in VS code, use shortcut  shift opt enter

# %% 
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/mass_cons_1000"
suffix = "random_ic"
nx = 64
max_it = 1000
max_it_CH = 10000
tol = 1.0e-5
@time main(max_it, max_it_CH, tol, outdir, suffix=suffix, initialize="random")

#%%

v = readdlm("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:1000, v, title="Discrete normalized energy \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Discrete Energy")
ylims!(0, 1)
savefig(p, "$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

v = readdlm("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:1000, v, title="Average Mass \n $(suffix)", linewidth=3, yformatter=:plain)
xlabel!("Time step")
ylabel!("Average Mass")
ylims!(v[1] - 0.000001, v[1] + 0.000001)
savefig(p, "$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

#%%
# TANGENT INITIAL CONDITIONS

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/mass_cons_1000_tan_IC_v4"
suffix = "tan_IC_v4"
nx = 128
max_it = 1000
max_it_CH = 10000
tol = 1.0e-5
@time main(max_it, max_it_CH, tol, outdir, suffix=suffix, initialize="function")
#%%
#Plotting tangent IC
v = readdlm("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:1000, v, title="Discrete normalized energy \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Discrete Energy")
ylims!(0, maximum(v))
savefig(p, "$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

v = readdlm("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:1000, v, title="Average Mass \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Average Mass")
ylims!(v[1] - 0.000001, v[1] + 0.000001)
savefig(p, "$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)


#%%
################ solve v5 #####################
# version 5 does not include global variables (const). All variables that are changed are defined locally in main_v5. 
using Plots
using DelimitedFiles
include("solver_type_stable_v5.jl")
#%%

outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/mass_cons_1000_tan_IC"
suffix = "tan_IC"
nx = 128
max_it = 1000
max_it_CH = 10000
tol = 1.0e-5
@time main_v5(nx, max_it, max_it_CH, tol, outdir, suffix=suffix, initialize="function")
#%%
#Plotting tangent IC
v = readdlm("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:1000, v, title="Discrete normalized energy \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Discrete Energy")
ylims!(0, maximum(v))
savefig(p, "$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

v = readdlm("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:1000, v, title="Average Mass \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Average Mass")
ylims!(v[1] - 0.000001, v[1] + 0.000001)
savefig(p, "$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

#%%
# Plotting ave mass and energy for 15000 time steps of tangent IC, grid sizes 128 <<<<<
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/tan_IC_ave_mass_128"
suffix = "tan_IC"
nx = 128
max_it = 15000
max_it_CH = 10000
@time main_v5(nx, max_it, max_it_CH, 0.01, outdir, suffix=suffix, initialize="function", overwrite=true)
@time main_v5(nx, max_it, max_it_CH, 0.001, outdir, suffix=suffix, initialize="function", overwrite=false)
@time main_v5(nx, max_it, max_it_CH, 0.0001, outdir, suffix=suffix, initialize="function", overwrite=false)
@time main_v5(nx, max_it, max_it_CH, 0.00001, outdir, suffix=suffix, initialize="function", overwrite=false)
@time main_v5(nx, max_it, max_it_CH, 0.000001, outdir, suffix=suffix, initialize="function", overwrite=false)


#%%
# Plotting ave mass and energy for 15000 time steps of tangent IC, grid sizes 256 <<<<<<
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/tan_IC_ave_mass_256"
suffix = "tan_IC"
nx = 256
max_it = 15000
max_it_CH = 10000
@time main_v5(nx, max_it, max_it_CH, 0.01, outdir, suffix=suffix, initialize="function", overwrite=true)
@time main_v5(nx, max_it, max_it_CH, 0.001, outdir, suffix=suffix, initialize="function", overwrite=false)
@time main_v5(nx, max_it, max_it_CH, 0.0001, outdir, suffix=suffix, initialize="function", overwrite=false)
@time main_v5(nx, max_it, max_it_CH, 0.00001, outdir, suffix=suffix, initialize="function", overwrite=false)

#%%
# Plotting ave mass and energy for 15000 time steps of geometric CPC, grid sizes 256 <<<<<<
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/CPC_256"
suffix = "CPC-geometric"
nx = 256
max_it = 10000
max_it_CH = 10000
@time main_v5(nx, max_it, max_it_CH, 0.00001, outdir, suffix=suffix, initialize="geometric", overwrite=true)
#%%
tol = 0.00001
v = readdlm("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:max_it, v, title="Discrete normalized energy \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Discrete Energy")
ylims!(0, maximum(v))
savefig(p, "$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

v = readdlm("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:max_it, v .- v[1], title="Average Mass \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Average Mass - Initial")
# ylims!(v[1] - 0.000001, v[1] + 0.000001)
savefig(p, "$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

#%%
# Plotting ave mass and energy for 15000 time steps of geometric CPC, grid sizes 256, M = 8 (to match C) <<<<<<
outdir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/CPC_256_M=8"
suffix = "CPC-geometric_M=8"
nx = 256
max_it = 15000
max_it_CH = 10000
@time main_v5(nx, max_it, max_it_CH, 0.00001, outdir, suffix=suffix, initialize="geometric", overwrite=true, M=8, ns=50)
#%%
tol = 0.00001
v = readdlm("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:max_it, v, title="Discrete normalized energy \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Discrete Energy")
ylims!(0, maximum(v))
savefig(p, "$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)

v = readdlm("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt")
p = plot(1:max_it, v .- v[1], title="Average Mass \n $(suffix)", linewidth=3)
xlabel!("Time step")
ylabel!("Average Mass - Initial")
# ylims!(v[1] - 0.000001, v[1] + 0.000001)
savefig(p, "$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).pdf")
display(p)