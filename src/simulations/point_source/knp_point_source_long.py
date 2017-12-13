import sys
import os
dirname, filename = os.path.split(os.path.abspath(__file__))
rel_path = "/../../"
print dirname+rel_path
sys.path.append(dirname+rel_path)
from impKNP import *
from dolfin import *
import time

p0 = Point(0,0,0)
p1 = Point(400e-6,400e-6,40e-6)
mesh = BoxMesh(p0, p1, 50, 50, 5)
geometry = Geometry(mesh)
simulator = Simulator(geometry)

print "loaded mesh and made spaces!"

def boundary(x, on_boundary):
    return on_boundary

lambda_o = 1.6 # ECS tortuousity, Chen & Nicholson 2000;


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
# f_Cl = Expression("x[0]*x[0]*(100-x[0])*(100-x[0])*t", t=0)
ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, boundary, "Cl")

init_cond_K = Expression('3', degree=4)
z_K = 1
D_K = 1.96e-9/lambda_o**2
init_K = init_cond_K
c_boundary_K = init_cond_K
# f_Cl = Expression("x[0]*x[0]*(100-x[0])*(100-x[0])*t", t=0)
ion_K = Ion(simulator, z_K, D_K, init_K, c_boundary_K, boundary, "K")



def mag_func(t):
    return 1.0*(t<1e0)*1e-9

def neg_mag_func(t):
    return -mag_func(t)


p1 = Point(280e-6, 200e-6, 20e-6)
p2 = Point(120e-6, 200e-6, 20e-6)

current_1 = Current(mag_func, ion_K)
currents = [current_1]
delta = Delta(p1, currents)
simulator.add_point_source(delta)

current_2 = Current(neg_mag_func, ion_K)
currents = [current_2]
delta = Delta(p2, currents)
simulator.add_point_source(delta)

dt = 1e-2
time_solver = Time_solver(simulator, dt, t_stop=2e0)
potential = KirchoffPotential(simulator)

print "initializing"

simulator.initialize_simulator()

print "initialized!"

# live_plotter = Live_plotter(simulator)

fname = dirname + "/knp_long_with_ca.h5"
notes = "This simulation considers a point source in a 3d grid, with knp"
state_saver = State_saver(fname,simulator, notes)

time_solver.solve()
