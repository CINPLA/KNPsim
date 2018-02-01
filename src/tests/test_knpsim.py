import pytest
from dolfin import *
import knpsim

def test_constant_solution():
    from knpsim import *
    from dolfin import *
    import time

    x0 = 0
    x1 = 100e-6
    xmid = 50e-6

    mesh = IntervalMesh(10000, x0, x1)
    geometry = Geometry(mesh)
    simulator = Simulator(geometry)

    def boundary(x, on_boundary):
        return on_boundary

    lambda_o = 1.6  # ECS tortuousity, Chen & Nicholson 2000;

    init_cond_X = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)

    D_X = 1
    z_X = 0

    init_X = init_cond_X
    c_boundary_X = init_cond_X
    ion_X = Ion(simulator, z_X, D_X, init_X, c_boundary_X, None, "X")

    dt = 1e-2
    time_solver = Time_solver(simulator, dt, t_stop=10)
    potential = ZeroPotential(simulator)
    simulator.initialize_simulator()

    for i in range(4):
        time_solver.solve_for_time_step()
