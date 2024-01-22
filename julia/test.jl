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