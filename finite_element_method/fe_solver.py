"""This is adapted from a demo by 2009 Garth N. Wells. Adapted by Sarah Groves Nov 15 2024.  """

# https://olddocs.fenicsproject.org/dolfin/1.4.0/python/demo/documented/cahn-hilliard/python/documentation.html
# install with conda install -c conda-forge fenics-dolfinx mpich pyvista
import random
from dolfinx import *


# Class representing the intial conditions
class InitialConditions(Expression):
    def __init__(self):
        random.seed(2 + MPI.rank(mpi_comm_world()))

    def eval(self, values, x):
        values[0] = 0.63 + 0.02 * (0.5 - random.random())
        values[1] = 0.0

    def value_shape(self):
        return (2,)


# Class for interfacing with the Newton solver
class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a

    def F(self, b, x):
        assemble(self.L, tensor=b)

    def J(self, A, x):
        assemble(self.a, tensor=A)


# Model parameters
lmbda = 1.0e-02  # surface parameter
dt = 5.0e-06  # time step
theta = 0.5  # time stepping family, e.g. theta=1 -> backward Euler, theta=0.5 -> Crank-Nicolson

# Form compiler options
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "quadrature"

# Create mesh and define function spaces
mesh = UnitSquareMesh(96, 96)
V = FunctionSpace(mesh, "Lagrange", 1)
ME = V * V

# Define trial and test functions
du = TrialFunction(ME)
q, v = TestFunctions(ME)

# Define functions
u = Function(ME)  # current solution
u0 = Function(ME)  # solution from previous converged step

# Split mixed functions
dc, dmu = split(du)
c, mu = split(u)
c0, mu0 = split(u0)

# Create intial conditions and interpolate
u_init = InitialConditions()
u.interpolate(u_init)
u0.interpolate(u_init)

# Compute the chemical potential df/dc
c = variable(c)
f = 100 * c**2 * (1 - c) ** 2
dfdc = diff(f, c)

# mu_(n+theta)
mu_mid = (1.0 - theta) * mu0 + theta * mu

# Weak statement of the equations
L0 = c * q * dx - c0 * q * dx + dt * dot(grad(mu_mid), grad(q)) * dx
L1 = mu * v * dx - dfdc * v * dx - lmbda * dot(grad(c), grad(v)) * dx
L = L0 + L1

# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(L, u, du)

# Create nonlinear problem and Newton solver
problem = CahnHilliardEquation(a, L)
solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6

# Output file
file = File("output.pvd", "compressed")

# Step in time
t = 0.0
T = 50 * dt
while t < T:
    t += dt
    u0.vector()[:] = u.vector()
    solver.solve(problem, u.vector())
    file << (u.split()[0], t)

plot(u.split()[0])
interactive()
