from fenics import *
from newton import Newton_manual
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



N = 1000
mesh = IntervalMesh(N, 0, 100)


# Defining functions and FEniCS stuff:
V = FunctionSpace(mesh, 'CG', 1)
R = FunctionSpace(mesh, 'R', 0)
W = MixedFunctionSpace([V, V, V, R])
v_1, v_2, v_phi, d = TestFunctions(W)



u = Function(W)
u_new = Function(W)

init_cond = Expression('10 + (10*(x[0] < 50))', degree=1)
assign(u.sub(0), interpolate(init_cond, V))
assign(u.sub(1), interpolate(init_cond, V))

# c1, c2, phi, dummy = split(u)
c1, c2, phi, dummy = u[0],u[1], u[2], u[3]
c1_new, c2_new, phi_new, dummy_new = u_new[0],u_new[1],u_new[2], u_new[3]
# c1_new, c2_new, phi_new, dummy_new = split(u_new)


# init_cond = Expression('10 + (100*(x[0] < 0.5))', degree=1)
# c1.assign(interpolate(init_cond, V))
# c2.assign(interpolate(init_cond, V))

# boundary conditions
def boundary(x, on_boundary):
    return on_boundary

bcs = [DirichletBC(W.sub(0), init_cond, boundary), DirichletBC(W.sub(1), init_cond, boundary)]



# Params:
F = 9.648e4 # Faradays constant, C/mol
T = 300 # Temperature, Kelvin
R = 8.314 # Rayleighs constant, J/(K*mol)
psi = R*T/F
eps_0 = 8.854 # Vacuum permitivity, pF/m
eps_r = 80 # Relative permitivity of water, no dimension
eps = eps_0*eps_r

D1 = 2.0 # diffusion coefficient
D2 = 1.0 # diffusion coefficient
z1 = 1 # valency
z2 = -1 # valency

dt = 1e-4 # time step, ms

# Form:
rho = F*(z1*c1_new + z2*c2_new)
# rho = Constant(0)
form = ((c1_new-c1)*v_1 + dt*inner(D1*nabla_grad(c1_new) + \
    D1*c1_new*z1*nabla_grad(phi_new)/psi, nabla_grad(v_1)))*dx + \
    ((c2_new-c2)*v_2 + dt*inner(D2*nabla_grad(c2_new) + \
    D2*c2_new*z2*nabla_grad(phi_new)/psi, nabla_grad(v_2)))*dx + \
    (eps*inner(nabla_grad(phi_new),nabla_grad(v_phi)) + dummy_new*v_phi + phi_new*d - rho*v_phi)*dx


dw = TrialFunction(W)
Jac = derivative(form, u_new, dw)
u_res = Function(W)

func_plot = Function(W)
c1plot = Function(V)
c2plot = Function(V)
phiplot = Function(V)

t = 0
for i in range(2000):
    t += dt
    Newton_manual(Jac, form, u_new, u_res,bcs=bcs, max_it=1000, atol = 1e-9, rtol=1e-6)
    # solve(form==0, u_new, bcs)
    # u.assign(u_new)
    assign(u, u_new)


    assign(c1plot, u_new.sub(0))
    assign(c2plot, u_new.sub(1))
    assign(phiplot, u_new.sub(2))

    plot(c1plot, title=str(t))
    plot(c2plot)
    plot(phiplot)


interactive()
