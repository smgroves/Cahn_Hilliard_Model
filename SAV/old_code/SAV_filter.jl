import Pkg
# Pkg.add("FFTW")
using FFTW
using Plots
using LinearAlgebra
# using Printf
# using Pkg; Pkg.add("CSV")
using CSV
using Pkg;
Pkg.add("Tables");
using Tables
# Pkg.add("PythonPlot")
# using PythonPlot
# pythonplot()
# using Pkg;
# Pkg.add(url="https://github.com/RobBlackwell/Parula.jl");
# using Parula

##############################
######### Functions ##########
##############################

#MATLAB to JULIA
# fft2 --> fft
#real(fft2(x)) --> real.(fft(x))
# ifft2 --> ifft

function f(phi)
    a = 0.5
    psi = (phi .+ 1) ./ 2
    fphi = 4 .* psi .^ 2 .* (6 * a .- 4 * (1 + a) * psi .+ 3 * psi .^ 2) ./ (3 * epsilon^2)
    return fphi
end

function r0_fun(phi0)
    r0 = sqrt(hx * hy * sum(f(phi0)) + C0)
    return r0
end

function fft2_filtered(x)
    return real.(fft(x))
end

function Lap(phi)
    Lphi = real.(ifft(k2 .* fft2_filtered(phi)))
    return Lphi
end

function A_inv_CN(phi)
    Ai = real.(ifft(fft2_filtered(phi) ./ (1 .+ dt / 2 * M * k4 - dt / 2 * gamma0 / epsilon^2 * M * k2)))
    return Ai
end


function df(phi)
    #  dfphi = (phi.^3-(1+gamma0)*phi)/epsilon^2;
    a = 0.5
    dpsidphi = 1 / 2
    psi = (phi .+ 1) ./ 2
    dfphi = 16 * (1 .- psi) .* (a .- psi) .* psi / epsilon^2 * dpsidphi
    return dfphi
end

function b_fun(phi)
    E1 = fft2_filtered(f(phi))
    b = df(phi) ./ sqrt(E1[1, 1] * hx * hy + C0)
    return b
end

function g_fun_CN(phi0, r0, b)
    bphi0 = fft2_filtered(b .* phi0)
    bphi0 = hx * hy * bphi0[1, 1]
    E1 = fft2_filtered(f(phi0))
    g = phi0 - dt / 2 * M * Lap(Lap(phi0)) + dt / 2 * gamma0 / epsilon^2 * M * Lap(phi0) + dt * M .* Lap(b) * (r0 - 1 / 4 * bphi0 - 1 / 2 * Beta * dt * r0 * (r0 - sqrt(E1[1, 1] * hx * hy + C0)))
    return g
end

function r_fun(phi, phi0, r0, b)
    bphi0 = fft2_filtered(b .* phi0)
    bphi0 = hx * hy * bphi0[1, 1]
    bphi = fft2_filtered(b .* phi)
    bphi = hx * hy * bphi[1, 1]

    E1 = fft2_filtered(f(phi0))
    r = r0 + 1 / 2 * (bphi - bphi0) - Beta * dt * r0 * (r0 - sqrt(E1[1, 1] * hx * hy + C0))
    return r
end


function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function print_mat(file, matrix)
    open(file, "a", lock=false) do f
        for i = axes(matrix, 1)
            for j = axes(matrix, 2)
                Printf.@printf(f, "%16.15f ", matrix[i, j])
            end
            println(f)
        end
    end
end
##############################
######## Parameters ##########
##############################

## Space parameters
global Lx = 1
global Ly = Lx
Nx = 2^7
Ny = Nx
global hx = Lx / Nx
global hy = Ly / Ny
para = Dict{String,Any}("Lx" => Lx)
para["Ly"] = Ly
para["Nx"] = Nx
para["Ny"] = Nx

x = hx .* (0:Nx-1)
y = hy .* (0:Ny-1)
xx, yy = meshgrid(x, y)

## Sarah's parameters
h = hx
h2 = h^2
gam = 8 * h / (2 * sqrt(2) * atanh(0.9)) # gam = epsilon here
max_it = 15000

## Interface and energy parameters
global epsilon = gam / 2^0
global C0 = 0
global M = epsilon^2
para["epsilon"] = epsilon
para["C0"] = C0
para["M"] = M

## Relaxation and stabilization
global Beta = 0
global gamma0 = 0
para["Beta"] = Beta
para["gamma0"] = gamma0

## Time parameters
global dt = 1e-1 * h2
global dtout = 50
global T = max_it * dt
para["dt"] = dt
para["dtout"] = dtout
para["T"] = T

##############################
###### Initialize data #######
##############################

#Initial data by importing
#phi0 = read from matrix

## Initial data by function
R = @.sqrt((xx - 0.5)^2 + (yy - 0.5)^2)
R0 = 0.1
eps_c = 0.01
psi0 = 0.5 * (1 .+ @.tanh((R0 .- R) / (2 * eps_c)))


# % phi0 = 0.05*sin(xx).*sin(yy) para.phi0=phi0
# % phi0 = 0.05*sin(x)'.*sin(y) para.phi0=phi0
# % r0 = r0_fun(phi0) para.r0=r0

# % C_cond=10
# % C_ref =10
# % psi0=psi0 * C_ref/C_cond

phi0 = 2 .* psi0 .- 1    # psi0=(phi0+1)/2
para["phi0"] = phi0


## Plot initial figure
p = Plots.contour(x, y, real.(psi0),
    fill=true,
    xlabel="X",
    ylabel="Y",
    title="Initial Conditions",
    size=(600, 550),
    c=:viridis)

Plots.display(p)

##############################
#### Crank Nicolson scheme ###
##############################

function CH2d_SAV_CN(para)
    # global dt epsilon k2 k4 C0 hx hy Lx Ly gamma0 Beta M
    T = para["T"]
    dt = para["dt"]
    Lx = para["Lx"]
    Ly = para["Ly"]
    Nx = para["Nx"]
    Ny = para["Ny"]
    epsilon = para["epsilon"]
    C0 = para["C0"]
    M = para["M"]
    gamma0 = para["gamma0"]
    Beta = para["Beta"]

    phi0 = para["phi0"]

    t = 0
    Nt = Int.(T / dt)
    E = zeros(Int.(Nt + 1))
    D = zeros(Int.(Nt + 1))
    O = zeros(Int.(Nt + 1))


    x = hx * (0:Nx-1)
    y = hy * (0:Ny-1)
    k_x = 1im * vcat(0:Nx÷2, (-Nx÷2+1:-1)) * (2 * pi / Lx)
    k_y = 1im * vcat(0:Ny÷2, (-Ny÷2+1:-1)) * (2 * pi / Ly)
    # % [kx,ky] = meshgrid(k_x,k_y)
    k_xx = k_x .^ 2
    k_yy = k_y .^ 2
    kxx, kyy = meshgrid(k_xx, k_yy)
    global k2 = kxx + kyy
    # % k2 = k_xx+k_yy
    global k4 = k2 .^ 2

    r0 = r0_fun(phi0)
    Phi = phi0
    # on each iteration, Phi = vcat(Phi, phi)

    #initial modified energy
    Lphi = Lap(phi0)
    E[1] = hx * hy * sum(1 / 2 .* phi0 .* (-Lphi)) + r0^2

    #Error
    D[1] = (r0 - sqrt(hx * hy * sum(f(phi0)))) / sqrt(hx * hy * sum(f(phi0)))
    O[1] = norm(imag.(fft(phi0)))

    # Time iteration
    for nt in 1:Nt
        t = t + dt
        phi_bar = A_inv_CN(phi0 + dt / 2 * M * Lap(df(phi0)))

        # Step 1
        b = b_fun(phi_bar)

        g = g_fun_CN(phi0, r0, b)

        AiLb = A_inv_CN(M * Lap(b))
        Aig = A_inv_CN(g)

        gamma = -fft2_filtered(b .* AiLb)
        gamma = gamma[1, 1] * hx * hy

        # Step 2
        bphi = fft2_filtered(b .* Aig)
        bphi = bphi[1, 1] * hx * hy / (1 + dt / 4 * gamma)

        # Step 3
        global phi = dt / 4 * bphi .* AiLb + Aig
        global r = r_fun(phi, phi0, r0, b)

        # Update phi0, phi1, r1
        phi0 = phi
        r0 = r

        # print out
        if mod(nt, dtout) == 0
            ll = Int.(nt / dtout)
            tt = Nx*ll+1:Nx*(ll+1)
            println(tt)
            Phi = cat(Phi, phi, dims=1)
        end


        #Modified energy
        E[1+nt] = hx * hy * sum(1 / 2 .* phi .* (-Lap(phi))) + r0^2


        # Error
        D[1+nt] = (r0 - sqrt(hx * hy * sum(f(phi)))) / sqrt(hx * hy * sum(f(phi)))
        O[1+nt] = norm(imag.(fft(phi)))
    end
    return phi, r, E, D, Phi, O
end


##############################
######## Run Solver ##########
##############################

# CN Solver (Phi collected every dtout results)
# [phi, ~, E, ~, Phi, O] = 
phi, r, E, D, Phi, O = CH2d_SAV_CN(para)



psi = (phi .+ 1) ./ 2
Psi = (Phi .+ 1) ./ 2

Plots.contour(x, y, real.(psi),
    fill=true,
    xlabel="X",
    ylabel="Y",
    title="Psi at time $T",
    size=(600, 550),
    c=:viridis)


CSV.write("./SAV/Julia/Phi_julia.csv", Tables.table(Phi), writeheader=false)

CSV.write("./SAV/Julia/D.csv", Tables.table(D), writeheader=false)

CSV.write("./SAV/Julia/E.csv", Tables.table(epsilon .* E), writeheader=false)