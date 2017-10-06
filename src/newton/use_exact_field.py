from fenics import *
from newton import Newton_manual
import math

"""
In this script we aim to solve the Poisson-Nernst-Planck equations, i.e. the
concentration dynamics of two ions, assuming they act by electrodiffusion.
Variables are
c1 - concentration of ion type 1
c2 - concentration of ion type 2
phi - the electric field
The dynamics equations are
dc1/dt = D(nabla^2(c1) + (1/psi)*nabla(c1*nabla(phi)))
dc2/dt = D(nabla^2(c1) + (1/psi)*nabla(c1*nabla(phi)))
nabla^2 phi = F(z1*c1 + z2*c2)/eps
D, psi, F, eps are physical constants (see below).
z1, z2 are the valencies of the ions.
In this example, the ions start with a high concentration on the left half of
an interval, and a low concentration on the right half. The initial
concentrations are equal.
We use dirichlet boundary conditions for the ions, and a pure von neuman
boundary for the electric field.
"""

set_log_level(99)



def run_mms(dt, N, end_time, theta=0):
    mesh = UnitIntervalMesh(N)

    # Defining functions and FEniCS stuff:
    degree = 1
    V = FunctionSpace(mesh, 'CG', degree)
    R = FunctionSpace(mesh, 'R', 0)
    W = MixedFunctionSpace([V, V, V, R])

    # Params:
    F = 9.648e4   # Faradays constant, C/mol
    T = 300       # Temperature, Kelvin
    R = 8.314     # Rayleighs constant, J/(K*mol)
    psi = R*T/F
    eps_0 = 8.854 # Vacuum permitivity, pF/m
    eps_r = 80    # Relative permitivity of water, no dimension
    eps = eps_0*eps_r

    D1 = 2.0 # diffusion coefficient
    D2 = 1.0 # diffusion coefficient
    z1 = 1   # valency
    z2 = -1  # valency
    time = 0

    # dt = 1e-3 # time step, ms



    x = SpatialCoordinate(mesh)
    t = Constant(0)
    phi_cc = ((sin(pi*x[0]))**2 - 0.5)*cos(t)**2
    # # phi_cc = 0
    c1_cc = cos(x[0])**3*cos(t)
    c2_cc = 1/z2*(-eps/F*div(nabla_grad(phi_cc)) - z1*(c1_cc))
    f1 = diff(c1_cc,t) - D1*div(nabla_grad(c1_cc) + (1.0/psi)*z1*c1_cc*nabla_grad(phi_cc))
    f2 = diff(c2_cc,t) - D2*div(nabla_grad(c2_cc) + (1.0/psi)*z2*c2_cc*nabla_grad(phi_cc))


    phi_e = Expression("(pow(sin(pi*x[0]), 2) - 0.5) * pow(cos(t),2)", t=time, degree=degree, domain=mesh)
    # phi_e = Expression(0, t=time)
    c1_e = Expression("pow(cos(x[0]), 3) * cos(t)", t=time, degree=degree)
    c2_e = Expression("1.0/z2*(-eps/F*2*pi*pi*pow(cos(t),2)*cos(2*pi*x[0]) - z1*pow(cos(x[0]), 3) * cos(t))",  \
        z2=z2, z1=z1, eps=eps, F=F, degree=degree, t=time)

    u_0 = Expression(("pow(cos(x[0]), 3) * cos(t)", "1.0/z2*(-eps/F*2*pi*pi*pow(cos(t),2)*cos(2*pi*x[0]) - z1*pow(cos(x[0]), 3) * cos(t))", "(pow(sin(pi*x[0]), 2) - 0.5) * pow(cos(t),2)", "0"), z2=z2, z1=z1, eps=eps, F=F, degree=degree, t=time)

    # P1 = FiniteElement('P', triangle, 1)
    # R = FiniteElement('R', triangle, 0)
    # element = MixedElement([P1, P1, P1, R])
    # W = FunctionSpace(mesh, element)

    v_1, v_2, v_phi, d = TestFunctions(W)


    # u = project(u_0, W)
    u = Function(W)
    u_new = Function(W)

    c1, c2, phi, dummy = split(u)
    c1_new, c2_new, phi_new, dummy_new = split(u_new)
    theta_float = theta
    theta = Constant(theta)
    c1_theta = theta*c1 + (1-theta)*c1_new
    c2_theta = theta*c2 + (1-theta)*c2_new
    phi_theta = theta*phi + (1-theta)*phi_new
    dummy_theta = theta*d + (1-theta)*dummy_new



    # assigner = FunctionAssigner(W.sub(0), V)



    assign(u.sub(0), interpolate(c1_e, V))
    assign(u.sub(1), interpolate(c2_e, V))
    assign(u.sub(2), interpolate(phi_e, V))
    # u_0 = Expression()

    # boundary conditions
    def boundary(x, on_boundary):
        return on_boundary

    bcs = [DirichletBC(W.sub(0), c1_e, "on_boundary"), DirichletBC(W.sub(1), c2_e, "on_boundary")]

    rho = F*(z1*c1_new + z2*c2_new)
    k = Constant(dt)
    form = ((c1_new-c1)*v_1 + k*inner(D1*nabla_grad(c1_theta) + \
        D1*c1_theta*z1*nabla_grad(phi_cc)/psi, nabla_grad(v_1)) - k*f1*v_1)*dx + \
        ((c2_new-c2)*v_2 + k*inner(D2*nabla_grad(c2_theta) + \
        D2*c2_theta*z2*nabla_grad(phi_cc)/psi, nabla_grad(v_2)) - k*f2*v_2)*dx + \
        (phi_new*v_phi + dummy_new*v_phi + phi_new*d - phi_e*v_phi)*dx

    dw = TrialFunction(W)
    Jac = derivative(form, u_new, dw)
    u_res = Function(W)

    tv = 0
    n_iter = int(end_time / dt)
    error_plot = Function(V)
    error_plot2 = Function(V)
    error_plot3 = Function(V)



    for i in range(n_iter):
        tv += (1-theta_float)*dt
        t.assign(tv)
        tv += theta_float*dt
        # f1.t = tv
        # f2.t = tv
        c1_e.t = tv
        c2_e.t = tv
        phi_e.t = tv
        # Newton_manual(Jac, form, u_new, u_res,bcs=bcs, max_it=100, atol = 1e-12, rtol=1e-12)

        solve(form==0, u_new, bcs)
        assign(error_plot,u_new.sub(0))
        # assign(error_plot,error_plot-project(c1_e, V))
        # plot(error_plot - project(c1_e,V))
        # interactive()

        # c1, c2, phi, dummy = u.split()
        # c1_new, c2_new, phi_new, dummy_new = u_new.split()
        # assign(c1, c1_new)
        # assign(c2, c2_new)
        # assign(phi, phi_new)
        # assign(dummy, dummy_new)

        u.assign(u_new)

        # assign(u.sub(0), u_new.sub(0))
        # assign(u.sub(1), u_new.sub(1))
        # assign(u.sub(2), u_new.sub(2))
        # assign(c2, c2_new)
        # assign(phi, phi_new)

    interactive()
    c1_e_f = project(c1_e, V)
    c1_sol = Function(V)
    assign(c1_sol, u.sub(0))

    c2_e_f = project(c2_e, V)
    c2_sol = Function(V)
    assign(c2_sol, u.sub(1))

    phi_e_f = project(phi_e, V)
    phi_sol = project(phi_cc, V)

    error_c1 = errornorm(c1_e_f, c1_sol, norm_type="l2", degree_rise=0)
    error_c2 = errornorm(c2_e_f, c2_sol, norm_type="l2", degree_rise=0)
    error_phi = errornorm(phi_e_f, phi_sol, norm_type="l2", degree_rise=0)

    # print "norms:"
    # print norm(u.sub(0), norm_type="l2")
    # print norm(u.sub(1), norm_type="l2")
    # print norm(u.sub(2), norm_type="l2")

    return error_c1, error_c2, error_phi, mesh.hmin()


def run_convergence(N_list, dt_list, theta=0):
    errors_c1 = []
    errors_c2 = []
    errors_phi = []
    h = []
    type_of_convergence = "Spatial" if len(N_list) > len(dt_list) else "Temporal"
    print "="*15, type_of_convergence, "="*15
    # end_time = 0.05
    end_time = dt_list[0]*10 if type_of_convergence == "Spatial" else max(dt_list)*4
    for N in N_list:
        for dt in dt_list:
            print N, dt
            error_c1, error_c2, error_phi, hmin = run_mms(dt, N, end_time,
                theta=theta)
            h.append(hmin)
            errors_c1.append(error_c1)
            errors_c2.append(error_c2)
            errors_phi.append(error_phi)

    N = max(len(N_list), len(dt_list)) - 1
    h = h if len(N_list) > len(dt_list) else dt_list

    print errors_c1
    print errors_c2
    print errors_phi

    print "Spatial convergence C1"
    for i in range(N):
        print math.log(errors_c1[i] / errors_c1[i+1]) / math.log(h[i] / h[i+1])

    print "Spatial convergence C2"
    for i in range(N):
        print math.log(errors_c2[i] / errors_c2[i+1]) / math.log(h[i] / h[i+1])

    print "Spatial convergence phi"
    for i in range(N):
        print math.log(errors_phi[i] / errors_phi[i+1]) / math.log(h[i] / h[i+1])

    print "\n"

if __name__ == '__main__':
        theta = 0
        run_convergence([10, 20, 40, 80], [1e-6], theta = theta)
        run_convergence([1000], [1e-2, 0.5e-2, 1e-3, 0.5e-3], theta=theta)
