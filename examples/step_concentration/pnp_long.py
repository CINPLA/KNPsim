"""
This script simulates a system of two ion species, starting with a step-
function concentration profile, using the PNP formalism. It uses a time step
of 1.0 ms, and runs from 0 to 10 seconds. The results are stored in the file
pnp_long.h5
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

# Set up the time solver
dt = 1e-3
time_solver = Time_solver(simulator, dt, t_stop=10, rtol=1e-3)

# set potential type
potential = PoissonPotential(simulator)

# Initialize simulator
print "initializing"
simulator.initialize_simulator()
print "initialized!"

# Set up state saver
fname = "pnp_long.h5"
notes = "This simulation considers a step concentration profile in 1D, solved \
    with PNP."

state_saver = State_saver(fname, simulator, notes)

# Run simulation
time_solver.solve()
