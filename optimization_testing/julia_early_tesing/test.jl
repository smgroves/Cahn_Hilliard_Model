using Pkg
Pkg.add("Plots")
julia
using Plots
x = 1:10;
y = rand(10); # These are the plotting data
plot(x, y, label="my label")

for i in range(1, 8)
    print(i, "\n")
end

a = 3
b = a
b += 1
print(b)
print(a)

g(x=[rand(), 2]) = push!(x, 3)

print(g())

m = [1 2; 4 5; 7 8]
print(m[3, 2])
print(m)

m = Matrix{Float64}(undef, 3, 4)

#%%
function test_recursion(ilevel, n_level)
    # print_mat("$(outdir)/uf_new_1_$(version)_$(suffix).csv", uf_new)
    # open("$(outdir)a_v4_full_$(suffix).txt", "a", lock=false) do file
    #     Printf.@printf(file, "First relax %d %d \n", ilevel, nxf)
    # end
    # open("$(outdir)mu_new_after_relax0_$(suffix).txt", "a", lock=false) do file
    # Printf.@printf(file, "First relax %d %d \n", ilevel, nxf)
    # end
    println("Outside if statement")
    # print_mat("$(outdir)/uf_new_2_$(version)_$(suffix).csv", uf_new)

    if ilevel < n_level
        println("inside if", ilevel)
        test_recursion(ilevel + 1, n_level)
        println("after recursion")
    end
end

test_recursion(1, 4)