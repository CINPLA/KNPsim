import sys
import random
from dolfin import Expression, Constant

# TODO: Update docstring
class Ion:
    """
    This class keeps track of mesh and the function space of the electro
    diffusion solver.
    """
    def __init__(self, simulator, charge, diffusion_coefficient,
                 initial_condition, boundary_concentration, boundary, name,
                 f=Expression("0", degree=2, t=0)):
        self.z = charge
        self.D = diffusion_coefficient
        self.name = name
        self.simulator = simulator
        self.id = random.randint(0, 100000000)
        self.boundary_condition = boundary_concentration
        simulator.add_ion(self)
        self.boundary = boundary
        self.f = f

        if isinstance(initial_condition, Expression):
            self.initial_condition = initial_condition
        elif isinstance(initial_condition, Constant):
            self.initial_condition = initial_condition
        else:
            print("initial_condition is unexpected type! Exiting...")
            sys.exit(1)
        # TODO: Accept more types...
