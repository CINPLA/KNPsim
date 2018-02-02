import sys
# FIXME: Why did this come from?
sys.path.append("/home/andreavs/Dropbox/knpsim_backward_euler_plain_SI")
from impKNP import *
from dolfin import *
import time
import scipy.io as sio
import numpy as np
# parameters['form_compiler']['optimize'] = True

ECSfrac = 0.2 # Fraction of tissue being extracellular space

d = sio.loadmat('revdata_100fold.mat')

mesh = Mesh()
f = HDF5File(mesh.mpi_comm(), "mesh_co.hdf5", 'r')
f.read(mesh, 'mesh', False)

geometry = Geometry(mesh)
simulator = Simulator(geometry)

y_coor = mesh.coordinates()[:,1]
ymin = y_coor.min()
ymax = y_coor.max()

print "loaded mesh and made spaces!"

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
# f_Cl = Expression("x[0]*x[0]*(100-x[0])*(100-x[0])*t", t=0)
ion_Cl = Ion(simulator, z_Cl, D_Cl, init_Cl, c_boundary_Cl, boundary, "Cl")

init_cond_K = Expression('3', degree=4)
z_K = 1
D_K = 1.96e-9/lambda_o**2
init_K = init_cond_K
c_boundary_K = init_cond_K
# f_Cl = Expression("x[0]*x[0]*(100-x[0])*(100-x[0])*t", t=0)
ion_K = Ion(simulator, z_K, D_K, init_K, c_boundary_K, boundary, "K")

factor = ECSfrac*10

n_points = len(np.transpose(d['x']))
ina = d['jna']*simulator.F*ion_Na.z/factor
icl = d['jx']*simulator.F*ion_Cl.z/factor
ica = d['jca']*simulator.F*ion_Ca.z/factor
ik = d['jk']*simulator.F*ion_K.z/factor
icap = d['icap']/factor
x = d['x'].reshape(n_points)
y = d['y'].reshape(n_points)
z = d['z'].reshape(n_points)

i_ion = ion_Na.z*ina + ion_K.z*ik + ion_Ca.z*ica + ion_Cl.z*icl
current_arrays = [ica, ina, icl, ik]

# n_x = 13
# z_pos = np.linspace(150e-6, 1350e-6,n_x)
dt = 1e-1

for i in range(n_points):
    p = Point(x[i], y[i], z[i])
    currents = []

    def cap_current(t, space_idx=i, icap=icap):
        time_idx = int(t/1e-4)
        time_stop_idx = int((t + dt)/1e-4)
        print "icap start", time_idx
        print "icap stop", time_stop_idx
        cc = np.mean(icap[space_idx, time_idx:time_stop_idx+1])
        return cc

    currents.append(Current(cap_current))
    for idx, ion in enumerate(simulator.ion_list):
        ion_current_array = current_arrays[idx]

        def ion_current(t, idx=idx, space_idx=i, ion_current_array=ion_current_array, ion=ion):
            time_idx = int(simulator.time_solver.t/1e-4)
            time_stop_idx = int((t + dt)/1e-4)
            print ion.name, time_idx
            print ion.name, time_stop_idx
            cc = np.mean(ion_current_array[space_idx, time_idx:time_stop_idx+1])
            # print cc
            return cc

        currents.append(Current(ion_current, ion))
    delta = Delta(p, currents)
    simulator.add_point_source(delta)

time_solver = Time_solver(simulator, dt, theta=1,t_start=8.1, t_stop=9.9)

# ns = Neuron_source(simulator, 'active.h5')
potential = KirchoffPotential(simulator)

simulator.initialize_simulator()

live_plotter = Live_plotter(simulator)

# fname = "/media/andreavs/datadrive/knp_sims_SI/hay_model/knp_3.h5"
# notes = "This simulation considers Hay model neuron in a cylindrical column, using knp"
# state_saver = State_saver(fname,simulator, notes)

time_solver.solve()

interactive()
