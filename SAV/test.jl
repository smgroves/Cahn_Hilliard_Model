Lx = 1
Ly = Lx
Nx = 2^4
Ny = Nx
hx = Lx / Nx
hy = Ly / Ny
x = hx .* (0:Nx-1)
y = hy .* (0:Ny-1)
A = x
A = vcat(A, y)
print(A)
print(x)
