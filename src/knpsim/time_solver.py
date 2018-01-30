from dolfin import *
from potential import KirchoffPotential, PoissonPotential
import sys
import time

from newton.newton import *
import numpy as np

class Time_solver:
    def __init__(self, simulator, dt, t_start=0., t_stop=1.0, atol=1e-9, rtol = 1e-6, max_iter=10):
        self.simulator = simulator
        self.dt = dt
        self.t_start = t_start
        self.t = t_start
        self.t_stop = t_stop
        simulator.set_time_solver(self)
        self.t_list = [t_start]
        self.atol = atol
        self.rtol = rtol
        self.max_iter = max_iter

    def set_time_step_size(self, dt):
        self.dt = dt

    def solve_for_time_step(self):
        t_new = self.t + self.dt
        self.t_list.append(t_new)

        # update source densities
        for ion in self.simulator.ion_list:
            ion.f.t = t_new

        pointsources = []
        for delta in self.simulator.deltas:
            I = 0
            for current in delta.currents:
                I_i = current.magnitude_function(self.t)
                I += I_i
                if current.ion != None:
                    pointsources.append(PointSource(self.simulator.geometry.W.sub(current.ion.index), delta.point, I_i/(self.simulator.F*current.ion.z)))
            if isinstance(self.simulator.potential,KirchoffPotential):
                pointsources.append(PointSource(self.simulator.geometry.W.sub(self.simulator.N+1), delta.point, I))


        # call solver:
        Newton_manual(self.simulator.Jac, self.simulator.form, \
            self.simulator.u_new, self.simulator.u_res,bcs=self.simulator.bcs, \
            deltas=pointsources, max_it=self.max_iter, atol = self.atol, rtol=self.rtol)
        # solve(self.simulator.form==0, self.simulator.u_new, self.simulator.bcs)

        # update old solution
        assign(self.simulator.u, self.simulator.u_new)
        if MPI.rank(mpi_comm_world()) == 0:
            print "Current time in simulation: " + str(self.t)
        self.t += self.dt
        if self.simulator.live_plotter:
            self.simulator.live_plotter.plot()
        if self.simulator.state_saver:
            self.simulator.state_saver.save_state()


    def solve(self):
        while self.t < self.t_stop:
            sim_t0 = time.clock()
            self.solve_for_time_step()
            sim_t1 = time.clock()
            if MPI.rank(mpi_comm_world()) == 0:
                print "The time step was solved in " + str(sim_t1-sim_t0) + " seconds."
        if self.simulator.state_saver:
            self.simulator.state_saver.finalize()
