# Julia implementation of a Multigrid solver for the Cahn-Hilliard equation
# use this version for printing out variables from relax and vcycle functions for troubleshooting Julia error
# Author: Sarah Groves
# December 17, 2023

# uncomment the following line to install the DataFrames package
# using Pkg
# Pkg.add("DataFrames")
# Pkg.add("BenchmarkTools")
# Pkg.add("StaticArrays");
# Pkg.add("ProfileView")

# using ProfileView
using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf


const nx = 256
const ny = 256
const n_level::Int = trunc(log(nx) / log(2.0) + 0.1)  # original c code uses natural log too
const c_relax::Int = 2  # number of SMOOTH relaxation operations defined as a global variable
const xleft = 0.0  # left x-coordinate defined as a global variable
const xright = 1.0  # right x-coordinate defined as a global variable
const yleft = 0.0  # left y-coordinate defined as a global variable
const yright = 1.0  # right y-coordinate defined as a global variable

# todo: check if these are needed; it appears that only ht2 (temp h^2) is used in the code
const h = xright / nx  # space step size defined as a global variable
const h2 = h^2 #space step size squared defined as a global variable
const dt = 0.1 * h2  # ∆t defined as a global variable
const gam = 4 * h / (2 * sqrt(2) * atanh(0.9))
const Cahn = gam^2  # ϵ^2 defined as a global variable
const tol = 1e-6
const ns = 10
global version = "v4" #undef -> 0

#dmatrix function: m = Array{Float64}(undef, 0, 0)

# double d2f(double c) { return 3.0 * c * c; } 

# laplacian function: laplacian(m, nx, ny, h2)
# tested with same results as Python implementation
function laplace(a, nxt, nyt)
    global xright
    global xleft
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

function source(c_old; nx=nx, ny=ny, dt=dt)
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

function relax(c_new, mu_new, su, sw, nxt, nyt; c_relax=c_relax, xright=xright, xleft=xleft, dt=dt, Cahn=Cahn)
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
                open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/f_full_initial.txt", "a", lock=false) do file
                    for i in 1:2
                        Printf.@printf(file, "%16.15f ", f[i])
                    end
                    println(file)
                end

                open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/sw.txt", "a", lock=false) do file
                    Printf.@printf(file, "%16.15f ", sw[i, j])
                end

                open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/su.txt", "a", lock=false) do file
                    Printf.@printf(file, "%16.15f ", su[i, j])
                end

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
                # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/cnew_before_update.txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", c_new[i, j])
                # end

                # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu_new_before_update.txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", mu_new[i, j])
                # end

                c_new[i, j] = (a[4] * f[1] - a[2] * f[2]) / det
                mu_new[i, j] = (-a[3] * f[1] + a[1] * f[2]) / det
                # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/a_v5_full.txt", "a", lock=false) do file
                #     for i in 1:4
                #         Printf.@printf(file, "%16.15f ", a[i])
                #     end
                #     println(file)
                # end
                # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/f_v5_full.txt", "a", lock=false) do file
                #     for i in 1:2
                #         Printf.@printf(file, "%16.15f ", f[i])
                #     end
                #     println(file)
                # end

                # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/cnew_after_update.txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", c_new[i, j])
                # end

                # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu_new_after_update.txt", "a", lock=false) do file
                #     Printf.@printf(file, "%16.15f ", mu_new[i, j])
                # end

                #         open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/d2f.txt", "a", lock=false) do file
                #             Printf.@printf(file, "%16.15f ", d2f)
                #         end

            end
            #     open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/d2f.txt", "a", lock=false) do file
            #         println(file)
            #     end

            #     open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/cnew.txt", "a", lock=false) do file
            #         println(file)
            #     end
            # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/cnew_after_update.txt", "a", lock=false) do file
            #     println(file)
            # end

            # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/cnew_before_update.txt", "a", lock=false) do file
            #     println(file)
            # end
            # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu_new_after_update.txt", "a", lock=false) do file
            #     println(file)
            # end

            # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu_new_before_update.txt", "a", lock=false) do file
            #     println(file)
            # end
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/sw.txt", "a", lock=false) do file
                println(file)
            end
            open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/su.txt", "a", lock=false) do file
                println(file)
            end
        end
    end

    # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/a_v4.txt", "a", lock=false) do file
    #     for i in 1:4
    #         Printf.@printf(file, "%16.15f ", a[i])
    #     end
    #     println(file)
    # end

    # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/f_v4.txt", "a", lock=false) do file
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

function nonL(c_new, mu_new, nxt, nyt; dt=dt, Cahn=Cahn)
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
function defect(uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc)
    ruc, rwc = nonL(uc_new, wc_new, nxc, nyc)
    ruf, rwf = nonL(uf_new, wf_new, nxf, nyf)
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

function vcycle(uf_new, wf_new, su, sw, nxf, nyf, ilevel)
    global n_level
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/wf_new_1_$(version).csv", wf_new)
    # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/a_v4_full.txt", "a", lock=false) do file
    #     Printf.@printf(file, "First relax %d %d \n", ilevel, nxf)
    # end
    # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu_new_after_relax0.txt", "a", lock=false) do file
    # Printf.@printf(file, "First relax %d %d \n", ilevel, nxf)
    # end
    uf_new, wf_new = relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax=c_relax, xright=xright,
        xleft=xleft, dt=dt, Cahn=Cahn)
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/wf_new_2_$(version).csv", wf_new)

    if ilevel < n_level
        nxc = trunc(Int64, nxf / 2)
        nyc = trunc(Int64, nyf / 2)
        uc_new, wc_new = restrict_ch(uf_new, wf_new, nxc, nyc)
        # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/wf_new_3_$(version).csv", wf_new)

        duc, dwc = defect(uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new, nxc, nyc)

        uc_def = copy(uc_new)
        wc_def = copy(wc_new)

        uc_def, wc_def = vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1)

        uc_def = uc_def - uc_new
        wc_def = wc_def - wc_new

        uf_def, wf_def = prolong_ch(uc_def, wc_def, nxc, nyc)

        uf_new = uf_new + uf_def
        wf_new = wf_new + wf_def
        # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/wf_new_4_$(version).csv", wf_new)
        # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu_new_after_relax0.txt", "a", lock=false) do file
        # Printf.@printf(file, "Second relax %d %d \n", ilevel, nxf)
        # end
        # open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/f_v4_full.txt", "a", lock=false) do file
        #     Printf.@printf(file, "Second relax %d %d \n", ilevel, nxf)
        # end
        uf_new, wf_new = relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax=c_relax, xright=xright,
            xleft=xleft, dt=dt, Cahn=Cahn)
        # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/wf_new_5_$(version).csv", wf_new)

    end
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/wf_new_6_$(version).csv", wf_new)

    return uf_new, wf_new
end

function error2(c_old, c_new, mu, nxt, nyt; dt=dt)
    rr = zeros(Float64, nxt, nyt)
    x = 0.0
    for i in 1:nxt
        for j in 1:nyt
            rr[i, j] = mu[i, j] - c_old[i, j]
        end
    end
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu.csv", mu)
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/c_old.csv", c_old)
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/rr_1.csv", rr)

    sor = laplace(rr, nxt, nyt)
    for i in 1:nxt
        for j in 1:nyt
            rr[i, j] = sor[i, j] - (c_new[i, j] - c_old[i, j]) / dt
        end
    end

    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/rr_2.csv", rr)

    for i in 1:nxt
        for j in 1:nyt
            x += rr[i, j]^2
        end
    end
    res2 = sqrt(x / (nxt * nyt))
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/res2_$(version).csv", res2)

    return res2
end

function initialization(nx, ny)
    phi = zeros(Float64, nx, ny)
    CPC_width = 5
    cohesin_width = 1
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

function cahn(c_old, c_new, mu; nx=nx, ny=ny, dt=dt, max_it_CH=10000, tol=1e-10)
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx=nx, ny=ny, dt=dt)
    # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/sc_$(version).csv", sc)

    while resid2 > tol && it_mg2 < max_it_CH
        c_new, mu = vcycle(c_new, mu, sc, smu, nx, ny, 1)
        # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/mu_7_$(version).csv", mu)
        # print_mat("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/julia/test_256/c_new_$(version).csv", c_new)

        resid2 = error2(c_old, c_new, mu, nx, ny, dt=dt)
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

function main(max_it, max_it_CH)
    println("nx = $nx, ny = $ny, dt = $dt, Cahn = $Cahn, max_it = $max_it,max_it_CH= $max_it_CH, ns = $ns, n_level = $n_level")
    mu = zeros(Float64, nx, ny)
    oc = initialization(nx, ny)
    nc = copy(oc)
    open("output_$(nx)_$(max_it).txt", "w", lock=false) do f
        writedlm(f, oc, " ")
    end
    for it in 1:max_it
        nc = cahn(oc, nc, mu, nx=nx, ny=ny, dt=dt, max_it_CH=max_it_CH, tol=tol)
        oc = copy(nc)
        if it % ns == 0
            open("output_$(nx)_$(max_it).txt", "a", lock=false) do f
                writedlm(f, oc, " ")
            end
            # println(it)
        end
    end
end



@time main(1, 1) #ignore first one with compile time

function write(max_it, max_it_CH)
    time_passed = @elapsed main(max_it, max_it_CH)
    open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Job_specs_all_py_c_julia.csv", "a", lock=false) do f
        writedlm(f, [c_relax Cahn "Julia" dt max_it max_it_CH n_level ns nx ny time_passed tol], ",")
    end
end

# write(10, 10000)
# write(100, 1000)
# write(1000, 1000)
# write(10000, 1000)

# write(10, 10000)
# write(100, 10000)
# write(1000, 10000)
# write(10000, 10000)

# write(10, 100000)
# write(100, 100000)
# write(1000, 100000)
# write(10000, 100000)



# @time main(100, 10000)

# @profview main(10, 10000)

# N=16: 4.647 ms (12239 allocations: 4.62 MiB)
# N=32: 12.360 ms (30829 allocations: 16.66 MiB)
