import sys
import os
dirname, filename = os.path.split(os.path.abspath(__file__))
rel_path = "/../../../"
print dirname+rel_path
sys.path.append(dirname+rel_path)
from impKNP import *
from dolfin import *
import time

x0 = 49.9e-6
x1 = 50.1e-6
xmid = 50e-6

mesh = IntervalMesh(10000, x0, x1)
geometry = Geometry(mesh)
simulator = Simulator(geometry)

print "loaded mesh and made spaces!"

def boundary(x, on_boundary):
    return on_boundary

lambda_o = 1.6 # ECS tortuousity, Chen & Nicholson 2000;

init_cond_Na = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)
z_Na = 1
D_Na = 1.33e-9/lambda_o**2
init_Na = init_cond_Na
c_boundary_Na = init_cond_Na
ion_Na = Ion(simulator, z_Na, D_Na, init_Na, c_boundary_Na, boundary, "Na")

init_cond_Cl = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)
z_Cl = -1
D_Cl = 2.03e-9/lambda_o**2
init_Cl = init_cond_Cl
c_boundary_Cl = init_cond_Cl
ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, boundary, "Cl")

dt = 1e-10
time_solver = Time_solver(simulator, dt, t_stop=1e-7, rtol=1e-5)
potential = PoissonPotential(simulator)

print "initializing"
simulator.initialize_simulator()
print "initialized!"

# live_plotter = Live_plotter(simulator)

fname = dirname + "/pnp_binary.h5"
notes = "This simulation considers a step concentration profile in 1D, solved with PNP."
state_saver = State_saver(fname,simulator, notes)

time_solver.solve()
