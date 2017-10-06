import sys
import random
from dolfin import *

class Ion:
    def __init__(self, simulator, charge, diffusion_coefficient, initial_condition, boundary_concentration, boundary, name, f=Expression("0", t=0)):
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
        else:
            print "initial_condition is unexpected type! Exiting..."
            sys.exit(1)
        # TODO: Accept more types...
