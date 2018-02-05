from knpsim import *
from dolfin import *
import time

# Set up mesh
p0 = Point(0, 0, 0)
p1 = Point(400e-6, 400e-6, 20e-6)
mesh = BoxMesh(p0, p1, 30, 30, 5)

# initialize geometry and simulator
geometry = Geometry(mesh)
simulator = Simulator(geometry)


def boundary(x, on_boundary):
    return on_boundary


# Set up the ion species
lambda_o = 1.6  # ECS tortuousity, Chen & Nicholson 2000;

init_cond_Ca = Expression('1.4', degree=4)
z_Ca = 2
D_Ca = 0.71e-9/lambda_o**2
init_Ca = init_cond_Ca
c_boundary_Ca = init_cond_Ca
ion_Ca = Ion(simulator, z_Ca, D_Ca, init_Ca, c_boundary_Ca, boundary, "Ca")

init_cond_Na = Expression('150', degree=4)
z_Na = 1
D_Na = 1.33e-9/lambda_o**2
init_Na = init_cond_Na
c_boundary_Na = init_cond_Na
ion_Na = Ion(simulator, z_Na, D_Na, init_Na, c_boundary_Na, boundary, "Na")

init_cond_Cl = Expression('155.8', degree=4)
z_Cl = -1
D_Cl = 2.03e-9/lambda_o**2
init_Cl = init_cond_Cl
c_boundary_Cl = init_cond_Cl
ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, boundary, "Cl")

init_cond_K = Expression('3', degree=4)
z_K = 1
D_K = 1.96e-9/lambda_o**2
init_K = init_cond_K
c_boundary_K = init_cond_K
ion_K = Ion(simulator, z_K, D_K, init_K, c_boundary_K, boundary, "K")


# Define current sources
def mag_func(t):
    ecsfrac = 0.2
    return 0.1*(t < 1e0)*1e-9/ecsfrac


def neg_mag_func(t):
    return -mag_func(t)


p1 = Point(120e-6, 200e-6, 10e-6)  # position of point source
p2 = Point(280e-6, 200e-6, 10e-6)  # position of point sink

current_1 = Current(mag_func, ion_K)
currents = [current_1]
delta = Delta(p1, currents)
simulator.add_point_source(delta)

current_2 = Current(neg_mag_func, ion_K)
currents = [current_2]
delta = Delta(p2, currents)
simulator.add_point_source(delta)

# Set up time solver
dt = 2e-3
time_solver = Time_solver(simulator, dt, t_stop=2e0)

# Set potential type
potential = KirchoffPotential(simulator)

# Initialize simulator
print("initializing")
simulator.initialize_simulator()
print("initialized!")


fname = "knp_point_source.h5"
notes = "This simulation considers a point source in a 3d grid, with knp"
state_saver = State_saver(fname, simulator, notes)

time_solver.solve()
