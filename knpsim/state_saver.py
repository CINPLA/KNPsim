from dolfin import *
import numpy as np
import h5py
import time


class State_saver:
    """
    This class saves the result of an electrodiffusion simulation. It is
    initialized with a system. It should be initialized after
    simulator.initialize_simulator(), but before time_solver.solve().

    Args:
        filename (str): Relative path to where the result should be saved.
        simulator (Simulator): The Simulator instance for the simulation.
        notes (str): Any notes on the simulation.
        save_step (int, optional): The result is saved every save_step-th step.
    """
    def __init__(self, filename, simulator, notes, save_step=1):
        self.simulator = simulator
        self.filename = filename
        simulator.state_saver = self
        self.hdf = HDF5File(simulator.geometry.mesh.mpi_comm(), filename, 'w')
        self.hdf.write(simulator.geometry.mesh, "geometry/mesh")
        self.hdf.write(Function(
            simulator.geometry.W),
            "geometry/MixedFunctionSpace")
        self.hdf.write(Function(
            simulator.geometry.V),
            "geometry/FunctionSpace")
        self.hdf.write(
            vertex_to_dof_map(simulator.geometry.V).astype(float),
            'vertex_to_dof_map')

        attribute_holder = Function(simulator.geometry.V)
        attribute_holder = self.hdf.write(attribute_holder, 'attributes')
        self.attributes = self.hdf.attributes('attributes')

        self.attributes['notes'] = notes
        self.attributes['number_of_ions'] = self.simulator.N
        ion_names = np.array([ion.name for ion in self.simulator.ion_list])

        # in case we only want every n-th step to be saved:
        self.save_step = save_step
        self.index = 0
        self.save_times = []

        self.hdf.write(self.simulator.u, '/initial_state')
        self.simtime_start = time.time()
        self.save_state()

    def save_state(self):
        """
        This function saves the potential and the concentrations at the current
        time step. It is called automatically be ´Time_solver.solve()´.
        """
        self.index += 1
        if self.index == self.save_step:
            self.hdf.write(self.simulator.u, "/solution",
                           self.simulator.time_solver.t)
            self.save_times.append(self.simulator.time_solver.t)
            self.index = 0

    def finalize(self):
        """
        This function should be called at the end of a simulation where
        ´State_saver´ has been used. It stores metadata such as the simulation
        duration. It is called automatically by `Time_solver.solve()`
        """
        self.simtime_end = time.time()
        simtime = self.simtime_end - self.simtime_start
        n_procs = MPI.size(mpi_comm_world())
        self.hdf.close()
        f = h5py.File(self.filename, 'r+')
        if MPI.rank(mpi_comm_world()) == 0:
            grp = f.create_group("simulation stats")
            grp.attrs['n_procs'] = n_procs
            grp.attrs['simtime'] = simtime
            grp = f.create_group("FunctionSpace")
            grp.attrs['space'] = self.simulator.geometry.space
            grp.attrs['order'] = self.simulator.geometry.order
            dset_time = f.create_dataset(
                "time",
                (np.size(self.simulator.time_solver.t_list), ),
                dtype='f')
            dset_time[...] = np.array(self.simulator.time_solver.t_list)
            dset_time = f.create_dataset(
                "save_time",
                (np.size(self.save_times), ),
                dtype='f')
            dset_time[...] = np.array(self.save_times)

            grp_ions = f.create_group("ions")
            for idx, ion in enumerate(self.simulator.ion_list):
                grp_ion = grp_ions.create_group(ion.name)
                grp_ion.attrs['diffusion constant'] = ion.D
                grp_ion.attrs['charge'] = ion.z
                grp_ion.attrs['index'] = ion.index

            f.close()
