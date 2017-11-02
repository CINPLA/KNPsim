"""
This simulation sets up a 1-D system with a length of 2 um (starting at 49 um,
ending at 51 um). The system has a larger concentration of K+ of the right side
of the system. starting at t=0, the solutions on the left and right side are
allowed to mix.
"""

import sys
import os
dirname, filename = os.path.split(os.path.abspath(__file__))
rel_path = "/../../../"
print dirname+rel_path
sys.path.append(dirname+rel_path)
from impKNP import *
from dolfin import *
import time

K_list = [2.5, 5, 10, 20,40, 80]

for k_concentration in K_list:
    mesh = IntervalMesh(2000, 49e-6, 51e-6)
    geometry = Geometry(mesh)
    simulator = Simulator(geometry)

    print "loaded mesh and made spaces!"

    def boundary(x, on_boundary):
        return on_boundary

    lambda_o = 1.6 # ECS tortuousity, Chen & Nicholson 2000;

    init_cond_Na = Expression('150 - left_con*(x[0]<50e-6)', degree=4, left_con=k_concentration)
    z_Na = 1
    D_Na = 1.33e-9/lambda_o**2
    init_Na = init_cond_Na
    c_boundary_Na = init_cond_Na
    ion_Na = Ion(simulator, z_Na, D_Na, init_Na, c_boundary_Na, boundary, "Na")

    init_cond_Cl = Expression('153', degree=4)
    z_Cl = -1
    D_Cl = 2.03e-9/lambda_o**2
    init_Cl = init_cond_Cl
    c_boundary_Cl = init_cond_Cl
    # f_Cl = Expression("x[0]*x[0]*(100-x[0])*(100-x[0])*t", t=0)
    ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, boundary, "Cl")


    init_cond_K = Expression('3 + left_con*(x[0]<50e-6)', degree=4, left_con=k_concentration)
    z_K = 1
    D_K = 1.96e-9/lambda_o**2
    init_K = init_cond_K
    c_boundary_K = init_cond_K
    # f_Cl = Expression("x[0]*x[0]*(100-x[0])*(100-x[0])*t", t=0)
    ion_K = Ion(simulator, z_K, D_K, init_K, c_boundary_K, boundary, "K")


    dt = 1e-10
    time_solver = Time_solver(simulator, dt, t_stop=2e-7, max_iter=6)
    potential = PoissonPotential(simulator)

    # def mag_func(t):
    #     return 1
    #
    # def neg_mag_func(t):
    #     return -mag_func(t)

    # delta = Delta(ion_Na, Point(20), mag_func)
    # simulator.add_point_source(delta)
    #
    # delta = Delta(ion_Cl, Point(80), mag_func)
    # simulator.add_point_source(delta)

    print "initializing"

    simulator.initialize_simulator()

    print "initialized!"

    # live_plotter = Live_plotter(simulator)

    fname = "/media/andreavs/datadrive/knp_sims_SI/advanced_step_function/pnp_left_con_"+ str(k_concentration) +".h5"
    notes = "This simulation considers a step concentration profile in 1D, solved with PNP."
    state_saver = State_saver(fname,simulator, notes)

    # live_plotter.plot()
    time_solver.solve()
