from fenics import *
import math

def run_mms(dt, N):
    # Mesh and spatial coordinates
    mesh = UnitCubeMesh(N, N, N)
    x = SpatialCoordinate(mesh)

    # Function space
    W = VectorFunctionSpace(mesh, "CG", 1)

    # Define functions
    u = TrialFunction(W)
    v = TestFunction(W)
    u_sol = Function(W)

    # Set time step and define variable t
    k = Constant(dt)
    t = Constant(dt)

    # Define exact solution
    u_e = Expression((
                    "sin(pow(x[0], 4))*sin(t)",            # x-direction
                    "cos(pow(x[1], 4))*cos(t)",            # y-direction
                    "cos(pow(x[2], 4))*sin(x[2])*sin(t)"   # z-direction
                     ), degree=4, t=0)
    u_x = sin(x[0]**4)*sin(t)
    u_y = cos(x[1]**4)*cos(t)
    u_z = cos(x[2]**4)*sin(x[2])*sin(t)
    u_vec = as_vector([u_x, u_y, u_z])

    # Create right hand side f
    f = diff(u_vec, t) - div(grad(u_vec))

    # Initial condition
    u_1 = interpolate(u_e, W)

    # Solve for f and exact bc
    F = 1/k * inner(u - u_1, v)*dx + inner(grad(u), grad(v))*dx - inner(f, v)*dx
    a = lhs(F)
    L = rhs(F)

    # Function to store solution
    u_sol = Function(W)

    # Boundary condition
    bcs = [DirichletBC(W, u_e, "on_boundary")]

    # Built-in solver
    t_ = 0
    for i in range(10):
        solve(a == L, u_sol, bcs)

        u_1.vector().zero()
        u_1.vector().axpy(1, u_sol.vector())

        # plot(u_1.sub(0))

        t_ += dt
        t.assign(t_)

    u_e.t = t_
    return errornorm(u_e, u_sol, norm_type="l2", degree_rise=3), mesh.hmin()


if __name__ == "__main__":
    error = []
    h = []
    N_list = [3, 5, 10, 15]
    dt = 1e-6

    # Spatial convergence test
    for N in N_list:
        error_h, hmin = run_mms(dt, N)
        h.append(hmin)
        error.append(error_h)

    print "Spatial convergence rate:"
    print h
    for i in range(len(N_list) - 1):
        print math.log(error[i] / error[i+1]) / math.log(h[i] / h[i+1])
