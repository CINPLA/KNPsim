"""
This script simulates a system of a single ion species, starting with a step-
function concentration profile, using the diffusion formalism. It uses a time
step of 1.0 ms, and runs from 0 to 10 seconds. The results are stored in the
file knp_zoom.h5.
"""
from knpsim import *
from dolfin import *
import time

# set up mesh
x0 = 0
x1 = 100e-6
xmid = 50e-6
mesh = IntervalMesh(10000, x0, x1)

# initialize geometry and simulator
geometry = Geometry(mesh)
simulator = Simulator(geometry)

print("loaded mesh and made spaces!")


def boundary(x, on_boundary):
    return on_boundary


# Set up the ion species
lambda_o = 1.6  # ECS tortuousity, Chen & Nicholson 2000;

init_cond_X = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)
D_Na = 1.33e-9/lambda_o**2
D_Cl = 2.03e-9/lambda_o**2
D_X = 2*D_Na*D_Cl/(D_Na + D_Cl)
z_X = 0
init_X = init_cond_X
c_boundary_X = init_cond_X
ion_X = Ion(simulator, z_X, D_X, init_X, c_boundary_X, None, "X")

# Set up the time solver
dt = 1e-3
time_solver = Time_solver(simulator, dt, t_stop=10)

# Set potential type
potential = ZeroPotential(simulator)

# Initialize simulator
simulator.initialize_simulator()

# Set up state saver
fname = "modified_diffusion.h5"
notes = "This simulation considers a step concentration profile in 1D, \
    solved with no field."
state_saver = State_saver(fname, simulator, notes)

# Run simulation
time_solver.solve()
