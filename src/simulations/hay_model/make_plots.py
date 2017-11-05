import sys

import h5py
from fenics import *
sys.path.append("/home/andreavs/Dropbox/knpsim_backward_euler_plain_SI")
sys.path.append("/home/andreavs/Dropbox/knpsim_backward_euler_plain_SI/plot_analysis_impKNP")
#sys.path.append("/home/andreavs/Dropbox/knpsim_backward_euler/plot_analysis_impKNP")
from plot_analysis_impKNP import *
import numpy as np
import scipy.io
from matplotlib import animation, rc

# Load data
analysistools = AnalysisTools('knp_3.h5')
cell_morphology = scipy.io.loadmat('cell_morphology.mat')
data = scipy.io.loadmat('revdata_100fold.mat')



xstart = cell_morphology['xstart'][0] + (57.28 + 34.16)/2
xend = cell_morphology['xend'][0] + (57.28 + 34.16)/2
ystart = cell_morphology['ystart'][0] + (17.62 + 19.07)/2
yend = cell_morphology['yend'][0] + (17.62 + 19.07)/2
N = cell_morphology['N'][0][0]



plt.figure()
[plt.plot([xstart[idx], xend[idx]], [ystart[idx], yend[idx]], c='k', linewidth=0.4) for idx in range(N)]
plt.title("Neuron morphology")
plt.axis("image")
plt.xlabel("x position (um)")
plt.ylabel("y position (um)")



plt.figure()
[plt.plot([xstart[idx], xend[idx]], [ystart[idx], yend[idx]], c='k', linewidth=0.4) for idx in range(N)]
x_dat = data['x'][0]*1e6
y_dat = data['y'][0]*1e6
z_dat = data['z'][0]*1e6
plt.scatter(x_dat, y_dat)
plt.title("Neuron morphology, with point sources")
plt.axis("image")
plt.xlabel("x position (um)")
plt.ylabel("y position (um)")

f, (ax1, ax2) = plt.subplots(1,2)
[ax1.plot([xstart[idx], xend[idx]], [ystart[idx], yend[idx]], c='k', linewidth=0.4) for idx in range(N)]
ax1.axis("image")
ax1.set_xlabel("x position (um)")
ax1.set_ylabel("y position (um)")

[ax2.plot([xstart[idx], xend[idx]], [ystart[idx], yend[idx]], c='k', linewidth=0.4) for idx in range(N)]
x_dat = data['x'][0]*1e6
y_dat = data['y'][0]*1e6
z_dat = data['z'][0]*1e6
ax2.scatter(x_dat, y_dat)
ax2.axis("image")
ax2.set_xlabel("x position (um)")
ax2.set_ylabel("y position (um)")


plt.show()



# Set up arrays
xmin = analysistools.mesh.coordinates()[:,0].min()
xmax = analysistools.mesh.coordinates()[:,0].max()
xlen = xmax - xmin

ymin = analysistools.mesh.coordinates()[:,1].min()
ymax = analysistools.mesh.coordinates()[:,1].max()
ylen = ymax - ymin

zmin = analysistools.mesh.coordinates()[:,2].min()
zmax = analysistools.mesh.coordinates()[:,2].max()
zmid = (zmax + zmin)/2

N_array = 100
padding = 0.01

x_array = np.linspace(xmin + padding*xlen, xmax - padding*xlen, N_array)
y_array = np.linspace(ymin + padding*ylen, ymax - padding*ylen, N_array)
X,Y = np.meshgrid (x_array,y_array)


# convert to um
X = X*1e6
Y = Y*1e6


# make K plot
t_idx = 700
Z = np.zeros([N_array,N_array,t_idx])
ion_idx = 3

for k in range(t_idx):
    u = Function(analysistools.W)
    c = Function(analysistools.V)
    analysistools.hdf.read(u, '/solution/vector_'+str(k))
    c.assign(project(u.sub(ion_idx), analysistools.V))
    for i,x in enumerate(x_array):
        for j,y in enumerate(y_array):
            p = Point(x,y,zmid)
            Z[j,i,k] = c(p)

idx = 0
plt.figure()
mesh = plt.pcolormesh(X, Y, Z[:,:,idx], vmin=3, vmax=4, cmap='Reds')
[plt.plot([xstart[idx], xend[idx]], [ystart[idx], yend[idx]], c='k', linewidth=0.4) for idx in range(N)]
cb = plt.colorbar()
plt.axis("image")
plt.title("K concentration at time " + str(analysistools.time_series[0]) + " s")

def update(frame):
#         global mesh
#         if frame % 25 == 0:  # print only every 25th frame
#             print(frame, end="..")
        current = Z[:,:,frame]
        mesh.set_array(current[:-1,:-1].ravel())
        plt.title("K concentration at time " + str(analysistools.time_series[frame]) + " s")
        return [mesh, ]

anim = animation.FuncAnimation(
    plt.gcf(),
    update,
    interval=60,
    frames=t_idx
)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
anim.save('k_concentration2.mp4')


# Make pot plot
Z = np.zeros([N_array,N_array,t_idx])
ion_idx = 4

for k in range(t_idx):
    u = Function(analysistools.W)
    c = Function(analysistools.V)
    analysistools.hdf.read(u, '/solution/vector_'+str(k))
    c.assign(project(u.sub(ion_idx), analysistools.V))
    for i,x in enumerate(x_array):
        for j,y in enumerate(y_array):
            p = Point(x,y,zmid)
            Z[j,i,k] = c(p)


idx = 0
plt.figure()
mesh = plt.pcolormesh(X, Y, Z[:,:,idx], vmin=-0.00003, vmax=0.00003, cmap='RdBu')
[plt.plot([xstart[idx], xend[idx]], [ystart[idx], yend[idx]], c='k', linewidth=0.4) for idx in range(N)]
plt.title("Potential at time:" + str(analysistools.time_series[0]) + " s")
cb = plt.colorbar()
plt.axis("image")

def update(frame):
#         global mesh
#         if frame % 25 == 0:  # print only every 25th frame
#             print(frame, end="..")
        current = Z[:,:,frame]
        plt.title("Potential at time:" + str(analysistools.time_series[frame]) + " s")
        mesh.set_array(current[:-1,:-1].ravel())

        return [mesh, ]

anim = animation.FuncAnimation(
    plt.gcf(),
    update,
    interval=60,
    frames=t_idx
)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
anim.save('potential2.mp4')
