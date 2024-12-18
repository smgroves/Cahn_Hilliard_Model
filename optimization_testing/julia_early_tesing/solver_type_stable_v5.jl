# Julia implementation of a Multigrid solver for the Cahn-Hilliard equation
# use this version for printing out variables from relax and vcycle functions for troubleshooting Julia error
# Author: Sarah Groves
# December 17, 2023

# uncomment the following line to install the DataFrames package
# using Pkg
# Pkg.add("DataFrames")
# Pkg.add("BenchmarkTools")
# Pkg.add("StaticArrays")
# Pkg.add("ProfileView")
# Pkg.add("DelimitedFiles")
# Pkg.add("LinearAlgebra")
# Pkg.add("Printf")


# using ProfileView
using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random
Random.seed!(1234) #note that when using a random seed, you must RESTART the REPL or run in a new instance for the runs to be the same. 

global version = "v5" #undef -> 0
const c_relax::Int = 2  # number of SMOOTH relaxation operations defined as a global variable
const xleft = 0.0  # left x-coordinate defined as a global variable
const xright = 1.0  # right x-coordinate defined as a global variable
const yleft = 0.0  # left y-coordinate defined as a global variable
const yright = 1.0  # right y-coordinate defined as a global variable

#dmatrix function: m = Array{Float64}(undef, 0, 0)

# double d2f(double c) { return 3.0 * c * c; } 

# laplacian function: laplacian(m, nx, ny, h2)
# tested with same results as Python implementation
function laplace(a, nxt, nyt)

    lap_a = zeros(Float64, nxt, nyt)
    ht2 = ((xright - xleft) / nxt)^2
    for i in 1:nxt
        for j in 1:nyt
            if i > 1
                dadx_L = (a[i, j] - a[i-1, j])
            else
                dadx_L = 0
            end
            if i < nxt
                dadx_R = (a[i+1, j] - a[i, j])
            else
                dadx_R = 0
            end
            if j > 1
                dady_B = (a[i, j] - a[i, j-1])
            else
                dady_B = 0
            end
            if j < nyt
                dady_T = (a[i, j+1] - a[i, j])
            else
                dady_T = 0
            end
            lap_a[i, j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2
        end
    end
    return lap_a
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
                # open("$(outdir)f_full_initial_$(suffix).txt", "a", lock=false) do file
                #     for i in 1:2
                #         Printf.@printf(file, "%16.15f ", f[i])
                #     end
                #     println(file)
                # end

                # open("$(outdir)sw_$(suffix).txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", sw[i, j])
                # end

                # open("$(outdir)su_$(suffix).txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", su[i, j])
                # end

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
                # open("$(outdir)cnew_before_update_$(suffix).txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", c_new[i, j])
                # end

                # open("$(outdir)mu_new_before_update_$(suffix).txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", mu_new[i, j])
                # end

                c_new[i, j] = (a[4] * f[1] - a[2] * f[2]) / det
                mu_new[i, j] = (-a[3] * f[1] + a[1] * f[2]) / det
                # open("$(outdir)a_v5_full_$(suffix).txt", "a", lock=false) do file
                #     for i in 1:4
                #         Printf.@printf(file, "%16.15f ", a[i])
                #     end
                #     println(file)
                # end
                # open("$(outdir)f_v5_full_$(suffix).txt", "a", lock=false) do file
                #     for i in 1:2
                #         Printf.@printf(file, "%16.15f ", f[i])
                #     end
                #     println(file)
                # end

                # open("$(outdir)cnew_after_update_$(suffix).txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", c_new[i, j])
                # end

                # open("$(outdir)mu_new_after_update_$(suffix).txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", mu_new[i, j])
                # end

                #         open("$(outdir)d2f_$(suffix).txt", "a", lock=false) do file
                #             Printf.@printf(file, "%16.15f ", d2f)
                #         end

            end
            #     open("$(outdir)d2f_$(suffix).txt", "a", lock=false) do file
            #         println(file)
            #     end

            #     open("$(outdir)cnew_$(suffix).txt", "a", lock=false) do file
            #         println(file)
            #     end
            # open("$(outdir)cnew_after_update_$(suffix).txt", "a", lock=false) do file
            #     println(file)
            # end

            # open("$(outdir)cnew_before_update_$(suffix).txt", "a", lock=false) do file
            #     println(file)
            # end
            # open("$(outdir)mu_new_after_update_$(suffix).txt", "a", lock=false) do file
            #     println(file)
            # end

            # open("$(outdir)mu_new_before_update_$(suffix).txt", "a", lock=false) do file
            #     println(file)
            # end
            # open("$(outdir)sw_$(suffix).txt", "a", lock=false) do file
            #     println(file)
            # end
            # open("$(outdir)su_$(suffix).txt", "a", lock=false) do file
            #     println(file)
            # end
        end
    end

    # open("$(outdir)a_v4_$(suffix).txt", "a", lock=false) do file
    #     for i in 1:4
    #         Printf.@printf(file, "%16.15f ", a[i])
    #     end
    #     println(file)
    # end

    # open("$(outdir)f_v4_$(suffix).txt", "a", lock=false) do file
    #     for i in 1:2
    #         Printf.@printf(file, "%16.15f ", f[i])
    #     end
    #     println(file)
    # end



    return c_new, mu_new
end

function restrict_ch(uf, vf, nxc, nyc)
    uc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    vc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    for i in 1:nxc
        for j in 1:nyc
            uc[i, j] = 0.25 * (uf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i - 1), round(Int, 2 * j)] + uf[round(Int, 2 * i), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i), round(Int, 2 * j)])
            vc[i, j] = 0.25 * (vf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i - 1), round(Int, 2 * j)] + vf[round(Int, 2 * i), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i), round(Int, 2 * j)])
        end
    end
    return uc, vc
end

function nonL(c_new, mu_new, nxt, nyt, dt, Cahn)
    ru = zeros(Float64, nxt, nyt)
    rw = zeros(Float64, nxt, nyt)
    lap_c = laplace(c_new, nxt, nyt)
    lap_mu = laplace(mu_new, nxt, nyt)
    for i in 1:nxt
        for j in 1:nyt
            ru[i, j] = c_new[i, j] / dt - lap_mu[i, j]
            rw[i, j] = mu_new[i, j] / dt - c_new[i, j]^3 + Cahn * lap_c[i, j]
        end
    end
    return ru, rw
end

# df(c) function: c.^3
# d2f(c) function: 3*c.^2
function defect(uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, Cahn)
    ruc, rwc = nonL(uc_new, wc_new, nxc, nyc, dt, Cahn)
    ruf, rwf = nonL(uf_new, wf_new, nxf, nyf, dt, Cahn)
    ruf = suf - ruf
    rwf = swf - rwf
    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = ruc + rruf
    dwc = rwc + rrwf
    return duc, dwc
end

function prolong_ch(uc, vc, nxc, nyc)
    uf = zeros(Float64, 2 * nxc, 2 * nyc)
    vf = zeros(Float64, 2 * nxc, 2 * nyc)
    for i in 1:nxc
        for j in 1:nyc
            uf[2*i-1, 2*j-1] = uc[i, j]
            uf[2*i-1, 2*j] = uc[i, j]
            uf[2*i, 2*j-1] = uc[i, j]
            uf[2*i, 2*j] = uc[i, j]
            vf[2*i-1, 2*j-1] = vc[i, j]
            vf[2*i-1, 2*j] = vc[i, j]
            vf[2*i, 2*j-1] = vc[i, j]
            vf[2*i, 2*j] = vc[i, j]
        end
    end
    return uf, vf
end

function vcycle(uf_new, wf_new, su, sw, nxf, nyf, ilevel, c_relax, xright, xleft, dt, Cahn, n_level)
    print_mat("$(outdir)/uf_new_1_$(version)_$(suffix).csv", uf_new)
    # open("$(outdir)a_v4_full_$(suffix).txt", "a", lock=false) do file
    #     Printf.@printf(file, "First relax %d %d \n", ilevel, nxf)
    # end
    # open("$(outdir)mu_new_after_relax0_$(suffix).txt", "a", lock=false) do file
    # Printf.@printf(file, "First relax %d %d \n", ilevel, nxf)
    # end
    uf_new, wf_new = relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright,
        xleft, dt, Cahn)
    print_mat("$(outdir)/uf_new_2_$(version)_$(suffix).csv", uf_new)

    if ilevel < n_level
        nxc = trunc(Int64, nxf / 2)
        nyc = trunc(Int64, nyf / 2)
        uc_new, wc_new = restrict_ch(uf_new, wf_new, nxc, nyc)
        print_mat("$(outdir)/uf_new_3_$(version)_$(suffix).csv", uf_new)

        duc, dwc = defect(uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, Cahn)
        print_mat("$(outdir)/duc.csv", duc)
        print_mat("$(outdir)/dwc.csv", dwc)

        uc_def = copy(uc_new)
        wc_def = copy(wc_new)

        uc_def, wc_def = vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1, c_relax, xright, xleft, dt, Cahn, n_level)


        uc_def = uc_def - uc_new
        print_mat("$(outdir)/coarse_grid_correction_u.csv", uc_def)
        wc_def = wc_def - wc_new
        print_mat("$(outdir)/coarse_grid_correction_w.csv", wc_def)

        uf_def, wf_def = prolong_ch(uc_def, wc_def, nxc, nyc)
        print_mat("$(outdir)/coarse_grid_correction_interpolated_u.csv", uf_def)
        print_mat("$(outdir)/coarse_grid_correction_interpolated_w.csv", wf_def)

        uf_new = uf_new + uf_def
        wf_new = wf_new + wf_def
        print_mat("$(outdir)/uf_new_4_$(version)_$(suffix).csv", uf_new)
        # open("$(outdir)mu_new_after_relax0_$(suffix).txt", "a", lock=false) do file
        # Printf.@printf(file, "Second relax %d %d \n", ilevel, nxf)
        # end
        # open("$(outdir)f_v4_full_$(suffix).txt", "a", lock=false) do file
        #     Printf.@printf(file, "Second relax %d %d \n", ilevel, nxf)
        # end
        uf_new, wf_new = relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright,
            xleft, dt, Cahn)
        print_mat("$(outdir)/uf_new_5_$(version)_$(suffix).csv", uf_new)

    end

    print_mat("$(outdir)/uf_new_6_$(version)_$(suffix).csv", uf_new)

    return uf_new, wf_new
end

function error2(c_old, c_new, mu, nxt, nyt, dt)
    rr = zeros(Float64, nxt, nyt)
    x = 0.0
    for i in 1:nxt
        for j in 1:nyt
            rr[i, j] = mu[i, j] - c_old[i, j]
        end
    end
    # print_mat("$(outdir)error2/mu_$(version)_$(suffix).csv", mu)
    # print_mat("$(outdir)error2/c_old_$(version)_$(suffix).csv", c_old)
    # print_mat("$(outdir)error2/c_new_$(version)_$(suffix).csv", c_old)
    # print_mat("$(outdir)error2/rr_1_$(version)_$(suffix).csv", rr)

    sor = laplace(rr, nxt, nyt)
    for i in 1:nxt
        for j in 1:nyt
            rr[i, j] = sor[i, j] - (c_new[i, j] - c_old[i, j]) / dt
        end
    end

    # print_mat("$(outdir)error2/rr_2_$(version)_$(suffix).csv", rr)

    for i in 1:nxt
        for j in 1:nyt
            x += rr[i, j]^2
        end
    end
    res2 = sqrt(x / (nxt * nyt))
    # print_mat("$(outdir)res2_$(version)_$(suffix).csv", res2)

    return res2
end

function initialize_geometric_CPC(nx, ny)
    phi = zeros(Float64, nx, ny)
    CPC_width = 20
    cohesin_width = 4
    for i in 1:nx
        for j in 1:ny
            if i > round(nx / 2) - CPC_width && i < round(nx / 2) + CPC_width
                if j > round(ny / 2) - CPC_width && j < round(ny / 2) + CPC_width
                    phi[i, j] = 1
                elseif i > round(nx / 2) - cohesin_width && i < round(nx / 2) + cohesin_width
                    phi[i, j] = 1
                else
                    phi[i, j] = -1
                end
            else
                phi[i, j] = -1
            end
        end
    end
    return phi
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function initialization_from_function(nx, ny, h)
    phi = zeros(Float64, nx, ny)
    x = h .* (0:nx-1)
    y = h .* (0:ny-1)
    xx, yy = meshgrid(x, y)
    R = @.sqrt((xx - 0.5)^2 + (yy - 0.5)^2)
    R0 = 0.1
    eps_c = 0.01
    psi0 = 0.5 * (1 .+ @.tanh((R0 .- R) / (2 * eps_c)))
    phi = 2 .* psi0 .- 1    # psi0=(phi0+1)/2
    return phi
end

function initialization_random(nx, ny)
    return (1 .- 2 .* rand(nx, ny))
end

function initialization_spinodal(nx, ny)
    return (rand([-1.0, 1.0], nx, ny))
end

function initialization_from_file(file, nx, ny, delim=',', transpose_matrix=true)
    phi = readdlm(file, delim, Float64)
    if size(phi) != (nx, ny)
        print("Warning: phi from file is wrong size: $(size(phi)) Expected: $(nx), $(ny)")
    end
    if transpose_matrix
        phi = transpose(phi)
    end
    return phi
end

function cahn(c_old, c_new, mu, nx, ny, dt, max_it_CH, tol, c_relax, xright, xleft, Cahn, n_level)
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx, ny, dt)
    # print_mat("$(outdir)sc_$(version)_$(suffix).csv", sc)

    while resid2 > tol && it_mg2 < max_it_CH

        c_new, mu = vcycle(c_new, mu, sc, smu, nx, ny, 1, c_relax, xright, xleft, dt, Cahn, n_level)
        # print_mat("$(outdir)mu_$(version)_$(suffix).csv", mu)
        # print_mat("$(outdir)c_new_$(version)_$(suffix).csv", c_new)

        resid2 = error2(c_old, c_new, mu, nx, ny, dt)
        it_mg2 += 1
    end
    # println(resid2)
    # println(it_mg2)
    return c_new
end

#will need to ensure matrices are floats not ints
function convert_matrix(matrix)
    if eltype(matrix) == Float64
        return matrix  # Matrix is already of type Float64, no need to convert
    elseif eltype(matrix) == Int
        return convert(Matrix{Float64}, matrix)  # Convert to Matrix{Float64}
    else
        throw(ArgumentError("Unsupported matrix element type"))
    end
end

#print_data function: writedlm(filename, m, " ")

function calculate_mass(phi, h2, nx, ny)
    ave_mass = sum(phi) / (h2 * nx * ny)
    return ave_mass
end

function f(phi)
    fphi = (1 / 4) .* ((phi .^ 2) .- 1) .^ 2
    return fphi
end

function calculate_discrete_energy(phi, h2, nx, ny, Cahn)
    a = h2 * sum(f(phi))
    sum_i = 0.0
    for i in 1:(nx-1)
        for j in 1:ny
            sum_i += (phi[i+1, j] - phi[i, j])^2

        end
    end
    b = (Cahn / 2) * sum_i
    sum_j = 0.0
    for i in 1:nx
        for j in 1:(ny-1)
            sum_j += (phi[i, j+1] - phi[i, j])^2
        end
    end
    c = (Cahn / 2) * sum_j
    E = a + b + c
    return E
end

function calculate_discrete_norm_energy(phi, phi0, h2, nx, ny, Cahn)
    E0 = calculate_discrete_energy(phi0, h2, nx, ny, Cahn)
    E = calculate_discrete_energy(phi, h2, nx, ny, Cahn)
    return E / E0
end

function main_v5(nx, max_it, max_it_CH, tol, outdir, ; suffix="", overwrite=true, print_phi=true, print_mass=true, print_e=true, initialize="function",
    dt=0, M=4, ns=50, gam=0, initial_file="", delim=',')
    while true
        ny = nx
        n_level::Int = trunc(log(nx) / log(2.0) + 0.1)  # original c code uses natural log too
        # todo: check if these are needed; it appears that only ht2 (temp h^2) is used in the code
        h = xright / nx  # space step size defined as a global variable
        h2 = h^2 #space step size squared defined as a global variable
        if dt == 0
            dt = 0.1 * h2  # ∆t defined as a global variable
        end

        if gam == 0
            gam = M * h / (2 * sqrt(2) * atanh(0.9))
        end
        Cahn = gam^2  # ϵ^2 defined as a global variable
        epsilon = gam / 2^0
        if isdir(outdir)
            if !isempty(outdir)
                if overwrite == false
                    println("Warning: Directory is not empty. Results may be appended to existing files. Are you sure you want to continue? (Y/N)")
                    input = readline()
                    if input == "Y" || input == "y"
                        println("Appending to any existing files.")
                    else
                        println("End.")
                        break
                    end
                else
                    println("Warning: overwriting directory with new files. Are you sure you want to continue? (Y/N)")
                    input = readline()
                    if input == "Y" || input == "y"
                        rm(outdir, recursive=true)
                        mkdir(outdir)
                    else
                        println("End.")
                        break
                    end
                end
            end
        else
            mkdir(outdir)
        end
        println("nx = $nx, ny = $ny, dt = $dt, epsilon = $gam, max_it = $max_it,max_it_CH= $max_it_CH, ns = $ns, n_level = $n_level")
        mu = zeros(Float64, nx, ny)
        # oc = initialization(nx, ny)
        if initialize == "random"
            oc = initialization_random(nx, ny)
        elseif initialize == "function"
            oc = initialization_from_function(nx, ny, h)
        elseif initialize == "geometric"
            oc = initialize_geometric_CPC(nx, ny)
        elseif initialize == "file"
            oc = initialization_from_file(initial_file, nx, ny, delim)
        elseif initialize == "spinodal"
            oc = initialization_spinodal(nx, ny)
        else
            println("Warning: initialize must be one of [random, function, geometric, file].")
        end
        nc = copy(oc)
        oc0 = copy(oc)
        if print_phi
            open("$(outdir)/phi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "w", lock=false) do f
                writedlm(f, oc, " ")
            end
        end
        for it in 1:max_it
            nc = cahn(oc, nc, mu, nx, ny, dt, max_it_CH, tol, c_relax, xright, xleft, Cahn, n_level)

            if print_mass
                print_mat("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt", calculate_mass(oc, h2, nx, ny))
            end
            # oc_psi = (oc .+ 1) ./ 2
            # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/ave_mass/ave_mass_psi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", calculate_mass(oc_psi))
            if print_e
                print_mat("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt", calculate_discrete_norm_energy(oc, oc0, h2, nx, ny, Cahn))
            end
            oc = copy(nc)
            if it % ns == 0
                if print_phi
                    open("$(outdir)/phi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "a", lock=false) do f
                        writedlm(f, oc, " ")
                    end
                end
                println(it)
            end
        end
        break
    end
end
