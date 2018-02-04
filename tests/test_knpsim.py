import pytest
from dolfin import *
from knpsim import *

#TODO: test Time_solver.set_time_step_size(dt)

def test_nofield_constant_solution():
    """
    Tests if an ion species with a constant initial concentration with no field
    will remain constant.
    """
    EPS = 1e-14
    import time

    x0 = 0
    x1 = 100e-6
    xmid = 50e-6

    mesh = IntervalMesh(10000, x0, x1)
    geometry = Geometry(mesh)
    simulator = Simulator(geometry)

    def boundary(x, on_boundary):
        return on_boundary

    init_cond_X = Expression('100', degree=4, xmid=xmid)

    D_X = 1
    z_X = 0

    init_X = init_cond_X
    c_boundary_X = init_cond_X
    ion_X = Ion(simulator, z_X, D_X, init_X, c_boundary_X, None, "X")

    dt = 1e-2
    time_solver = Time_solver(simulator, dt, t_stop=10)
    potential = ZeroPotential(simulator)
    simulator.initialize_simulator()

    c = simulator.u.sub(0)
    start_c = project(c, geometry.V)

    for i in range(4):
        time_solver.solve_for_time_step()
        c = project(simulator.u.sub(0), geometry.V)
        err = project(c-start_c, geometry.V)
        err_norm = norm(err)
        assert abs(err_norm) < EPS
