from dolfin import *


class Delta:
    """
    A class to hold the position and currents at a point source.
    """
    def __init__(self, point, currents):
        self.point = point
        self.currents = currents


class Current:
    """
    A class to hold the currents, used by Delta. magnitude_function should be
    callable.
    """
    def __init__(self, magnitude_function, ion=None):
        self.magnitude_function = magnitude_function
        self.ion = ion
