import sys
import os
dirname, filename = os.path.split(os.path.abspath(__file__))
rel_path = "/../../../"
print dirname+rel_path
sys.path.append(dirname+rel_path)
from impKNP import *
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

lambda_o = 1.6 # ECS tortuousity, Chen & Nicholson 2000;

init_cond_X = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)

D_Na = 1.33e-9/lambda_o**2
D_Cl = 2.03e-9/lambda_o**2

D_X = 2*D_Na*D_Cl/(D_Na + D_Cl)

z_X = 0
init_X = init_cond_X
c_boundary_X = init_cond_X
ion_X = Ion(simulator, z_X, D_X, init_X, c_boundary_X, boundary, "X")

dt = 1e-3
time_solver = Time_solver(simulator, dt, t_stop=1)
potential = ZeroPotential(simulator)

print "initializing"
simulator.initialize_simulator()
print "initialized!"

# live_plotter = Live_plotter(simulator)

fname = dirname + "/zero_binary_long.h5"
notes = "This simulation considers a step concentration profile in 1D, solved with no field."
state_saver = State_saver(fname,simulator, notes)

time_solver.solve()
