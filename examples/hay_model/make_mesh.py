from fenics import *
from mshr import *
import scipy.io
import numpy as np

data = scipy.io.loadmat('neuron_input_1.hdf5')

if MPI.rank(MPI.comm_world) == 0:
    print(data.keys())

x = data['x']
y = data['y']
z = data['z']

xlen = x.max() - x.min()
ylen = y.max() - y.min()
zlen = z.max() - z.min()

if MPI.rank(MPI.comm_world) == 0:
    print(xlen, ylen, zlen)

padding_fraction = 0.1

x0_mesh = x.min() - padding_fraction*xlen
x1_mesh = x.max() + padding_fraction*xlen

y0_mesh = y.min() - padding_fraction*ylen
y1_mesh = y.max() + padding_fraction*ylen

z0_mesh = z.min() - padding_fraction*zlen
z1_mesh = z.max() + padding_fraction*zlen

xmid_mesh = (x.min() + x.max())/2.
ymid_mesh = (y.min() + y.max())/2.
zmid_mesh = (z.min() + z.max())/2.

r_mesh = np.sqrt((xlen/2.)**2 + (zlen/2.)**2)
r_mesh += 2*r_mesh


p_top = Point(xmid_mesh, y0_mesh, zmid_mesh)
p_bottom = Point(xmid_mesh, y1_mesh, zmid_mesh)
domain = Cylinder(p_top, p_bottom, r_mesh, r_mesh)
# domain = Box(Point(x0_mesh, y0_mesh, z0_mesh), Point(x1_mesh, y1_mesh, z1_mesh))

resolution = 40

mesh = generate_mesh(domain, resolution)

if MPI.rank(MPI.comm_world) == 0:
    print(mesh.coordinates().shape)

f = HDF5File(mesh.mpi_comm(), "mesh_co.hdf5", 'w')
f.write(mesh, 'mesh')
