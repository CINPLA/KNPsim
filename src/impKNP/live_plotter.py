from dolfin import *
import numpy as np

class Live_plotter:
    def __init__(self, simulator):
        simulator.live_plotter = self
        self.simulator = simulator
        self.total_charge = Function(self.simulator.geometry.V)
        self.plot_functions = []
        for i in range(len(self.simulator.ion_list)+1):
            self.plot_functions.append(Function(self.simulator.geometry.V))



    def plot(self):
        rho = Function(self.simulator.geometry.V)

        for (idx, ion) in enumerate(self.simulator.ion_list):
            rho += ion.z*ion.c
            assign(self.plot_functions[idx],self.simulator.u.sub(idx))
            plot(self.plot_functions[idx], title=ion.name + ", t=" + str(self.simulator.time_solver.t))

        assign(self.plot_functions[self.simulator.N],self.simulator.u.sub(self.simulator.N))
        plot(self.plot_functions[self.simulator.N], title="Potential")
        self.total_charge.assign(project(rho,self.simulator.geometry.V))
        plot(self.total_charge, title="total charge")
