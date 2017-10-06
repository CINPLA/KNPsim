import h5py
import numpy as np
from dolfin import *
from simulator import Simulator
from delta import *

class Neuron_source:
    def __init__(self, simulator, h5path):
        h5_file = h5py.File(h5path, 'r')
        self.h5_file = h5_file
        self.simulator = simulator


        cell = h5_file['cell']
        time = np.array(cell['t'])
        data_dt = time[1] - time[0]
        ion_input = self.h5_file['ion_currents']

        for i in range(cell['x'].size):
            xmid = cell['x'][i]
            ymid = cell['y'][i]
            zmid = cell['z'][i]

            p = Point(xmid, ymid, zmid)
            cap_current_array = cell['icap'][i,:]

            def cap_current(t, i=i):
                if self.simulator.time_solver.dt <= data_dt:
                    tt = np.abs(time - t)
                    idx = np.argmin(tt)
                    cc = cap_current_array[idx]
                else:
                    curr_time = self.simulator.time_solver.t
                    new_time = curr_time + self.simulator.time_solver.dt
                    idx_list = np.where(np.logical_and(time >= curr_time, time <= new_time))
                    cc = np.mean(cap_current_array[idx_list])
                return cc

            currents = [Current(cap_current)]
            for ion in self.simulator.ion_list:
                def ion_current(t, i=i, ion=ion):
                    ion_current_array = ion_input['i_' + ion.name][i,:]
                    if self.simulator.time_solver.dt <= data_dt:
                        tt = np.abs(time - t)
                        idx = np.argmin(tt)
                        cc = ion_current_array[idx]
                    else:
                        curr_time = self.simulator.time_solver.t
                        new_time = curr_time + self.simulator.time_solver.dt
                        idx_list = np.where(np.logical_and(time >= curr_time, time <= new_time))
                        cc = np.mean(ion_current_array[idx_list])
                    return cc

                currents.append(Current(ion_current, ion))

            delta = Delta(p, currents)
            simulator.add_point_source(delta)



class Compartment:
    def __init__(self, neuron, idx):
        self.index = idx
        self.neuron = neuron
        cell = neuron.h5_file['cell']
        self.rmid = [cell['x'][idx], cell['y'][idx], cell['z'][idx]]
        self.rstart = [cell['xstart'][idx], cell['ystart'][idx], cell['zstart'][idx]]
        self.rend = [cell['xend'][idx], cell['yend'][idx], cell['zend'][idx]]

if __name__=='__main__':
    h5path = 'active.h5'
    ns = Neuron_source(h5path)
