# This script was used to test different alternative ICs for the critical radius simulations on Rivanna.
using Dates

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function initialize_two_halves(nx, ny, h; R0=0.1, gam=0.01)
    phi = zeros(Float64, nx, ny)
    x = h .* (0:nx-1)
    y = h .* (0:ny-1)
    xx, yy = meshgrid(x, y)
    x_center = 0.5
    y_center1 = 0.0
    y_center2 = 1.0

    R1 = @. sqrt((xx - x_center)^2 + (yy - y_center1)^2)
    R2 = @. sqrt((xx - x_center)^2 + (yy - y_center2)^2)

    #the +1 is necessary because we are superimposing two phase fields and the background of both is -1, so everything with be 1 less than we want. 
    phi = @.tanh((R0 .- R1) / (sqrt(2) * gam)) .+ @.tanh((R0 .- R2) / (sqrt(2) * gam)) .+ 1

    return phi
end

using Plots

nx, ny = 100, 100
h = 0.01
R0 = 0.1
gam = 0.01

phi = initialize_top(nx, ny, h; R0=R0, gam=gam)
heatmap(phi, color=:viridis)
