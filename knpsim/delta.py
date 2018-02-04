from dolfin import *


class Delta:
    """
    A class to hold the position and currents at a point source.

    Args:
        point (Point): A FEniCS Point object holding the location of the point
            source.
        currents (list): a list of currents.
    """
    def __init__(self, point, currents):
        self.point = point
        self.currents = currents


class Current:
    """
    A class to hold the currents, used by Delta.

    Args:
        magnitude_function (function): a function that takes exactly one
            argument and returns a float (i.e., the current amplitude as a
            function of time).
        ion (Ion, optional): which ion species the current belongs to.
    """
    def __init__(self, magnitude_function, ion=None):
        self.magnitude_function = magnitude_function
        self.ion = ion
