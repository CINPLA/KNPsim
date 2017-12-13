import sys

import h5py
from fenics import *
from analysistools import *
plt.rc('text', usetex=True)

a_pnp_zoom = AnalysisTools("fig_1_data/pnp_binary_zoom.h5")
a_knp_zoom = AnalysisTools("fig_1_data/knp_binary_zoom.h5")

a_pnp_long = AnalysisTools("fig_1_data/pnp_binary_long.h5")
a_knp_long = AnalysisTools("fig_1_data/knp_binary_long.h5")

a_nofield = AnalysisTools("fig_1_data/nofield_binary_long.h5")
a_mod_diff = AnalysisTools("fig_1_data/diffusion_x_binary_long.h5")

u = Function(a_pnp_zoom.W)
a_pnp_zoom.hdf.read(u, '/solution/vector_'+str(2))

uu = u[0]

print uu(Point(1e-9))
