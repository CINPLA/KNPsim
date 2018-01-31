from dolfin import *


class Delta:
    def __init__(self, point, currents):
        self.point = point
        self.currents = currents


class Current:
    def __init__(self, magnitude_function, ion=None):
        self.magnitude_function = magnitude_function
        self.ion = ion
