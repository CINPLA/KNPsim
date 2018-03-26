import os.path
from knpsim import *
from dolfin import *
import time
import h5py
import numpy as np
# parameters['form_compiler']['optimize'] = True


assert os.path.isfile('neuron_input_1.h5'), "neuron data doesn't exist!"
assert os.path.isfile('mesh_cylinder.hdf5'), "mesh doesn't exist!"
assert os.path.isfile('knp_column_neuron.h5'), "KNP data doesn't exist!"

ECSfrac = 0.2  # Fraction of tissue being extracellular space

d = h5py.File('neuron_input_1.h5', 'r')
mesh = Mesh()
f = HDF5File(mesh.mpi_comm(), "mesh_cylinder.hdf5", 'r')
f.read(mesh, 'mesh', False)

geometry = Geometry(mesh)
simulator = Simulator(geometry)

def boundary(x, on_boundary):
    return on_boundary



lambda_o = 1.6 # ECS tortuousity, Chen & Nicholson 2000;

init_cond_Ca = Expression('1.4', degree=4)
z_Ca = 2
D_Ca = 0.71e-9/lambda_o**2
init_Ca = init_cond_Ca
c_boundary_Ca = init_cond_Ca
ion_Ca = Ion(simulator, z_Ca, D_Ca, init_Ca, c_boundary_Ca, boundary, "Ca")


init_cond_Na = Expression('150', degree=4)
z_Na = 1
D_Na = 1.33e-9/lambda_o**2
init_Na = init_cond_Na
c_boundary_Na = init_cond_Na
ion_Na = Ion(simulator, z_Na, D_Na, init_Na, c_boundary_Na, boundary, "Na")

init_cond_Cl = Expression('155.8', degree=4)
z_Cl = -1
D_Cl = 2.03e-9/lambda_o**2
init_Cl = init_cond_Cl
c_boundary_Cl = init_cond_Cl
ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, boundary, "Cl")

init_cond_K = Expression('3', degree=4)
z_K = 1
D_K = 1.96e-9/lambda_o**2
init_K = init_cond_K
c_boundary_K = init_cond_K
ion_K = Ion(simulator, z_K, D_K, init_K, c_boundary_K, boundary, "K")

factor = 1

n_points = len(np.transpose(d['x']))
x = np.array(d['x']).reshape(n_points)
y = np.array(d['y']).reshape(n_points)
z = np.array(d['z']).reshape(n_points)

current_arrays = ['ica', 'ina', 'ix', 'ik']

dt = 1e-4

for i in range(n_points):
    p = Point(x[i], y[i], z[i])
    currents = []

    def cap_current(t, space_idx=i, icap=d['icap']):
        time_idx = int(t/1e-4)
        time_stop_idx = int((t + dt)/1e-4)
        # print "icap start", time_idx
        # print "icap stop", time_stop_idx
        cc = np.mean(icap[space_idx, time_idx:time_stop_idx+1])/ECSfrac
        return cc

    currents.append(Current(cap_current))
    for idx, ion in enumerate(simulator.ion_list):
        ion_current_array = current_arrays[idx]

        def ion_current(t, idx=idx, space_idx=i, ion_current_array=d[current_arrays[idx]], ion=ion):
            time_idx = int(t/1e-4)
            time_stop_idx = int((t + dt)/1e-4)
            # print ion.name, time_idx
            # print ion.name, time_stop_idx
            cc = np.mean(ion_current_array[space_idx, time_idx:time_stop_idx+1])/ECSfrac
            # print cc
            return cc

        currents.append(Current(ion_current, ion))
    delta = Delta(p, currents)
    simulator.add_point_source(delta)

time_solver = Time_solver(simulator, dt, t_start=80, t_stop=80+200e-3)
potential = KirchoffPotential(simulator)
simulator.initialize_simulator()
t_n = 800
g = HDF5File(mesh.mpi_comm(), "knp_column_neuron.h5", 'r')
u = Function(geometry.W)
g.read(u, '/solution/vector_'+str(t_n))

assign(simulator.u.sub(0), project(u.sub(0), geometry.V))
assign(simulator.u.sub(1), project(u.sub(1), geometry.V))
assign(simulator.u.sub(2), project(u.sub(2), geometry.V))
assign(simulator.u.sub(3), project(u.sub(3), geometry.V))

fname = "late_column_neuron.h5"
notes = "This simulation considers Hay model neuron in a cylindrical column, using knp and accounting for ECSfrac. It starts at 80 s, and uses a small time step"
state_saver = State_saver(fname,simulator, notes)

time_solver.solve()
