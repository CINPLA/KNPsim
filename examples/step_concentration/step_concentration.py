from knpsim import *
from dolfin import *
from argparse import ArgumentParser
import time


def str2bool(arg):
    if arg.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif arg.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def step_concentration(zoom, potential, dt, T, x0, x1, relax, modified, rtol):
    # Set up mesh
    xmid = (x0 + x1) / 2
    mesh = IntervalMesh(10000, x0, x1)

    # Initialize geometry and simulator
    geometry = Geometry(mesh)
    simulator = Simulator(geometry)

    # Define boundary
    boundary = (lambda x, on_boundary: on_boundary) if zoom else None

    # Initial condition
    init_cond = Expression('(140 + 10*(x[0]>=xmid))', degree=4, xmid=xmid)

    # Define constants
    lambda_o = 1.6  # ECS tortuousity, Chen & Nicholson 2000;
    z_Cl = -1
    D_Cl = 2.03e-9/lambda_o**2
    z_Na = 1
    D_Na = 1.33e-9/lambda_o**2
    D_X = 2*D_Na*D_Cl/(D_Na + D_Cl)
    z_X = 0

    # Include ion species
    if not modified:
        ion_Na = Ion(simulator, z_Na, D_Na, init_cond, init_cond, boundary,
                     "Na")
        ion_Cl = Ion(simulator, z_Cl, D_Cl, init_cond, init_cond, boundary,
                     "Cl")
    else:
        ion_X = Ion(simulator, z_X, D_X, init_cond, init_cond, boundary, "X")

    # Set up time solver
    time_solver = Time_solver(simulator, dt, t_stop=T, rtol=rtol, relax=relax)

    # Set potential type
    print(potential)
    if potential == "kirchoff":
        potential = KirchoffPotential(simulator)
        solver = "KNP"
    elif potential == "poisson":
        potential = PoissonPotential(simulator)
        solver = "PNP"
    else:
        potential = ZeroPotential(simulator)
        solver = "zero_potential"

    # Initialize simulator
    simulator.initialize_simulator()

    # Set up state saver
    zoom = "zoom" if zoom else "long"
    fname = solver.lower() + "_" + zoom + ".h5"
    fname = fname if not modified else fname = "modified_diffusion.h5"
    notes = "This simulation considers a step concentration profile in 1D," + \
        " solved with %s." % solver.replace("_", " ")

    state_saver = State_saver(fname, simulator, notes)

    # Run simulation
    time_solver.solve()


if __name__ == "__main__":
    description = """This script simulates a system of two ion species, starting with a step-
function concentration profile, using either PNP or KNP formalism."""
    parser = parser = ArgumentParser(description=description)

    parser.add_argument("-z", "--zoom",
                        default=False,
                        dest="zoom",
                        type=str2bool,
                        help="?")

    parser.add_argument("-p", "--potential",
                        default="kirchoff",
                        dest="potential",
                        choices=["kirchoff", "poisson", "zero"],
                        type=str,
                        help="?")

    parser.add_argument("-d", "--dt",
                        default=1e-3,
                        dest="dt",
                        type=float,
                        help="?")

    parser.add_argument("-T", "--time",
                        default=10,
                        dest="T",
                        type=float,
                        help="?")

    parser.add_argument("-s", "--start",
                        default=-50e-6,
                        dest="x0",
                        type=float,
                        help="?")

    parser.add_argument("-e", "--end",
                        default=50e-6,
                        dest="x1",
                        type=float,
                        help="?")

    parser.add_argument("-w", "--relax",
                        default=0.95,
                        dest="relax",
                        type=float,
                        help="?")

    parser.add_argument("-m", "--modified",
                        default=False,
                        dest="modified",
                        type=bool,
                        help="?")

    parser.add_argument("-r", "--rtol",
                        default=5e-3,
                        dest="rtol",
                        type=float,
                        help="?")

    args = parser.parse_args()

    step_concentration(args.zoom, args.potential, args.dt, args.T, args.x0,
                       args.x1, args.relax, args.modified, args.rtol)
