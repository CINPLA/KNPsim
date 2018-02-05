from knpsim import *
from dolfin import *
import time

x0 = 0
x1 = 100e-6
xmid = 50e-6

mesh = IntervalMesh(10000, x0, x1)
geometry = Geometry(mesh)
simulator = Simulator(geometry)

print "loaded mesh and made spaces!"


def boundary(x, on_boundary):
    return on_boundary


lambda_o = 1.6  # ECS tortuousity, Chen & Nicholson 2000;

init_cond_Na = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)
z_Na = 1
D_Na = 1.33e-9/lambda_o**2
init_Na = init_cond_Na
c_boundary_Na = init_cond_Na
ion_Na = Ion(simulator, z_Na, D_Na, init_Na, c_boundary_Na, None, "Na")

init_cond_Cl = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)
z_Cl = -1
D_Cl = 2.03e-9/lambda_o**2
init_Cl = init_cond_Cl
c_boundary_Cl = init_cond_Cl
ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, None, "Cl")

dt = 1e-3
time_solver = Time_solver(simulator, dt, t_stop=10, rtol=1e-3)
potential = PoissonPotential(simulator)

print "initializing"
simulator.initialize_simulator()
print "initialized!"

# live_plotter = Live_plotter(simulator)

fname = dirname + "/pnp_binary_long.h5"
notes = """
This simulation considers a step concentration profile in 1D, solved with PNP.
"""
state_saver = State_saver(fname, simulator, notes)

time_solver.solve()
