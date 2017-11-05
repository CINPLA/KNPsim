from neuron import h
import numpy as np
import scipy.io
from pylab import *
import time
import sys
from os.path import exists
import matplotlib.pyplot as plt

morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"

import LFPy
cell_parameters = {
    'morphology' : morphology_file,
    # 'delete_sections' : [False]
}

cell = LFPy.Cell(**cell_parameters)


proximalpoint = 400
distalpoint = 620

nsegs = 20
dt = 0.025
syngmax = 0.000042
dtsave = 0.1

# Initialize the model
h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
load_file("import3d.hoc")
objref L5PC
load_file(\""""+biophys_file+"""\")
load_file(\""""+template_file+"""\")
L5PC = new L5PCtemplate(\""""+morphology_file+"""\")
dtsave = """+str(dtsave)+"""
objref st1, synlist
st1 = new IClamp(0.5)
st1.del = 700
L5PC.soma st1
synlist = new List()
objref isyn,tvec,sl
isyn = new Vector()
tvec = new Vector()
sl = new List()
double siteVec[2]
sl = L5PC.locateSites("apic","""+str(distalpoint)+""")
maxdiam = 0
for(i=0;i<sl.count();i+=1){
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {
    j = i
    maxdiam = dd
  }
}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
access L5PC.apic[siteVec[0]]
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
{vsoma = new Vector()}
{casoma = new Vector()}
{vdend = new Vector()}
{cadend = new Vector()}
{vdend2 = new Vector()}
{cadend2 = new Vector()}
access L5PC.soma
{vsoma.record(&v(0.5),dtsave)}
{casoma.record(&cai(0.5),dtsave)}
access L5PC.apic[siteVec[0]]
{vdend.record(&v(siteVec[1]),dtsave)}
{cadend.record(&cai(siteVec[1]),dtsave)}
access L5PC.soma
{isoma = new Vector()}
{isoma.record(&st1.i,dtsave)}
forsec L5PC.all {
  nseg = """+str(nsegs)+"""
}
dt = """+str(dt)+"""
""")



apical_nsec = 109
basal_nsec = 84
somatic_nsec = 1
axonal_nsec = 2

xxstart = []
xxend = []
yystart = []
yyend = []
for i in range(apical_nsec):
    h("access L5PC.apic["+str(i)+"]")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    xxstart.append(h.tmpvarx)
    xxend.append(h.tmpvarx2)
    yystart.append(h.tmpvary)
    yyend.append(h.tmpvary2)
        # print h.tmpvarx

for i in range(somatic_nsec):
    h("access L5PC.soma["+str(i)+"]")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    print h.tmpvarx, h.tmpvarx2, h.tmpvary, h.tmpvary2

[plt.plot([xxstart[idx], xxend[idx]], [yystart[idx], yyend[idx]], c='k', linewidth=0.4) for idx in range(len(xxstart))]



offsetx = (57.28 + 34.16)/2
offsety = (17.62 + 19.07)/2
#
[plt.plot([cell.xstart[idx]+offsetx, cell.xend[idx]+offsetx], [cell.ystart[idx]+offsety, cell.yend[idx]+offsety], linewidth=0.4, c='b') for idx in xrange(cell.totnsegs)]
plt.show()
