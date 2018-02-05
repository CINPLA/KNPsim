import sys
import random
from dolfin import Expression, Constant


class Ion:
    """
    This class keeps track of mesh and the function space of the electro
    diffusion solver. It also holds the boundary condition, and ion parameters.

    Args:
        simulator (Simulator): The simulator isinstance.
        charge (int): The valence of the ion species.
        diffusion_coefficient (float): The diffusion coefficient of the ion
            species.
        initial_condition (Expression): A FEniCS Expression of the initial
            condition for the ion species.
        boundary_concentration (Expression): The boundary values for the ion
            species
        boundary (function): The boundary of the domain
        name (str): The name of the ion species (typically 'K', 'Ca', etc)
        f (Expression, optional): A source term for the ion species.

        Note:
            f is not used for point source terms. Use the Delta class for this.
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
        elif isinstance(initial_condition, Function):
            self.initial_condition = initial_condition
        else:
            print("initial_condition is unexpected type! Exiting...")
            sys.exit(1)
        # TODO: Accept more types...
