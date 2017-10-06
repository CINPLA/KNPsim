from fenics import *
import math
from newton import *

error = []
h = []
N_list = [3, 5, 10, 15]
for N in N_list:
    # Mesh and spatial coordinates
    mesh = UnitCubeMesh(N, N, N)
    x = SpatialCoordinate(mesh)

    # Function space
    W = VectorFunctionSpace(mesh, "CG", 1)

    # Define functions
    #u = TrialFunction(W)
    u = Function(W)
    v = TestFunction(W)
    u_sol = Function(W)

    # Define exact solution
    u_e = Expression((
                    "sin(pow(x[0], 4))",            # x-direction
                    "cos(pow(x[1], 4))",            # y-direction
                    "cos(pow(x[2], 4))*sin(x[2])"   # z-direction
                     ), degree=4)
    u_x = sin(x[0]**4)
    u_y = cos(x[1]**4)
    u_z = cos(x[2]**4)*sin(x[2])
    u_vec = as_vector([u_x, u_y, u_z])

    # Create right hand side f
    f = div(grad(u_vec))

    # Solve for f and exact bc
    a = -inner(grad(u), grad(v))*dx
    L = inner(f,v)*dx

    F = a - L

    u_sol = Function(W)
    bcs = [DirichletBC(W, u_e, "on_boundary")]

    # Built-in solver
    # solve(F == 0, u, bcs)


    # Manual newton
    w = TrialFunction(W)
    J = derivative(F, u, w)
    Newton_manual(J, F, u, u_sol, bcs=bcs)

    # plot(interpolate(u_e,W))
    # plot(u_sol, interactive=True)

    error.append(errornorm(u_e, u, norm_type="l2", degree_rise=3))
    h.append(mesh.hmin())

# FIXME: Visualize u_e and u_sol to see where the error is. On the boundary?
#        The built-in solver is working, and gives a perfect convergence
#        (although N is very small)

print "Convergence rate:"
print error
print h
for i in range(len(N_list) - 1):
    print math.log(error[i] / error[i+1]) / math.log(h[i] / h[i+1])
