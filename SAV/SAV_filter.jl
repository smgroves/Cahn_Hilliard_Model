import Pkg
# Pkg.add("FFTW")
using FFTW


# Pkg.add("PythonPlot")
using Plots
# using PythonPlot
# pythonplot()


function f(phi)
    a = 0.5
    psi = (phi + 1) / 2
    fphi = 4 * (psi .^ 2) .* (6 * a - 4 * (1 + a) * psi + 3 * psi .^ 2) / (3 * epsilon^2)
    return fphi
end

function r0_fun(phi0)
    r0 = @.sqrt(hx * hy * sum(sum(f(phi0))) .+ C0)
    return r0
end

function fft2_filtered(x)
    return real.(fft(x))
end

function Lap(phi)
    Lphi = real.(ifft(k2 .* fft2_filtered(phi)))
    return Lphi
end

## Space parameters
global Lx = 1
global Ly = Lx
para = Dict{String,Any}("Lx" => Lx)
para["Ly"] = Ly
Nx = 2^4
para["Nx"] = Nx
Ny = Nx
para["Ny"] = Nx
global hx = Lx / Nx
global hy = Ly / Ny
x = hx .* (0:Nx-1)
y = hy .* (0:Ny-1)


function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

xx, yy = meshgrid(x, y)

## Sarah's parameters
h = hx
h2 = h^2
gam = 8 * h / (2 * sqrt(2) * atanh(0.9)) # gam is our epsilon
max_it = 15000

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
Plots.contour(x, y, real.(psi0),
    fill=true,
    xlabel="X",
    ylabel="Y",
    title="Initial Conditions",
    size=(600, 550))

# display(p)

## Interface and energy parameters
global epsilon = gam / 2^0
para["epsilon"] = epsilon
global C0 = 0
para["C0"] = C0
global M = epsilon^2
para["M"] = M

## Relaxation and stabilization
global Beta = 0
para["Beta"] = Beta
global gamma0 = 0
para["gamma0"] = gamma0


## Time parameters
global dt = 1e-1 * h2
para["dt"] = dt
dtout = 50
para["dtout"] = dtout
T = max_it * dt
para["T"] = T

## CN Solver (phi collected every dtout results)
# [phi, ~, E, ~, Phi, O] = CH2d_SAV_CN(para)
# psi = (phi + 1) / 2
# Psi = (Phi + 1) / 2
# Plots.contour(x, y, real.(psi),
#     fill=true,
#     xlabel="X",
#     ylabel="Y",
#     title="psi")

## Crank Nicolson scheme
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
    Nt = round(T / dt)
    E = zeros(Int.(Nt + 1))

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

    # % phi  = phi0
    # % r0 = para.r0
    r0 = r0_fun(phi0)
    # % r  = r0
    Phi = phi0
    # on each iteration, Phi = vcat(Phi, phi)

    #initial modified energy
    Lphi = Lap(phi0)
    E[1] = hx * hy * sum(sum(1 / 2 * phi0 .* (-Lphi))) .+ r0^2
end

CH2d_SAV_CN(para)