from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from neuron import h
import numpy as np
import scipy.io
from pylab import *
import time
import sys
import os
from os import path

myseed = int(sys.argv[1])
seed(myseed)

directory = 'simulation_data'
if not os.path.exists(directory):
    os.mkdir(directory)

directory = 'simulation_data/seed_' + str(myseed)
if not os.path.exists(directory):
    os.mkdir(directory)


morphology_file = "morphologies/cell1.asc"
biophys_file = "models/L5PCbiophys3.hoc"
template_file = "models/L5PCtemplate.hoc"

v0 = -80
ca0 = 0.0001

synlambda = 5.0  # frequency of synaptic inputs (Hz)
syntau = 2.0     # decay time (ms)
proximalpoint = 400
distalpoint = 620
BACdt = 5.0

nsegs = 20
dt = 0.025
syngmax = 0.000042
tstop = 10000
Nsynlocs = 10000
synloctype = 2

singleSimT = 200

Nsims = int(1.0*tstop/singleSimT + 0.9999)
dtsave = 0.1

# Initialize the model
h("""
load_file("stdlib.hoc")
load_file("stdrun.hoc")
load_file("import3d.hoc")
objref L5PC
load_file("{biophys_file}")
load_file("{template_file}")
L5PC = new L5PCtemplate("{morphology_file}")
dtsave = {dtsave}
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
sl = L5PC.locateSites("apic",{distalpoint})
maxdiam = 0
for(i=0;i<sl.count();i+=1){{
  dd1 = sl.o[i].x[1]
  dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
  if (dd > maxdiam) {{
    j = i
    maxdiam = dd
  }}
}}
siteVec[0] = sl.o[j].x[0]
siteVec[1] = sl.o[j].x[1]
access L5PC.apic[siteVec[0]]
objref vsoma, vdend, recSite, vdend2, isoma, cadend, cadend2, casoma
{{vsoma = new Vector()}}
{{casoma = new Vector()}}
{{vdend = new Vector()}}
{{cadend = new Vector()}}
{{vdend2 = new Vector()}}
{{cadend2 = new Vector()}}
access L5PC.soma
{{vsoma.record(&v(0.5),dtsave)}}
{{casoma.record(&cai(0.5),dtsave)}}
access L5PC.apic[siteVec[0]]
{{vdend.record(&v(siteVec[1]),dtsave)}}
{{cadend.record(&cai(siteVec[1]),dtsave)}}
access L5PC.soma
{{isoma = new Vector()}}
{{isoma.record(&st1.i,dtsave)}}
forsec L5PC.all {{
  nseg = {nsegs}
}}
dt = {dt}
""".format(**dict(dt=dt, dtsave=dtsave, biophys_file=biophys_file,
                  template_file=template_file, morphology_file=morphology_file,
                  distalpoint=distalpoint, nsegs=nsegs)))

# Initialize the variables where the transmembrane currents will be saved

# DO NOT CHANGE - These are numbers from the morphology file
apical_nsec = 109
basal_nsec = 84
somatic_nsec = 1
axonal_nsec = 2
all_nsec = apical_nsec + basal_nsec + somatic_nsec + axonal_nsec

apical_nsegs = apical_nsec*nsegs
basal_nsegs = basal_nsec*nsegs
somatic_nsegs = somatic_nsec*nsegs
axonal_nsegs = axonal_nsec*nsegs

all_nsegs = apical_nsegs + basal_nsegs + somatic_nsegs + axonal_nsegs

tmp_dict = dict(nsegs=nsegs, axonal_nsec_nsegs=axonal_nsec*nsegs,
                apical_nsec_nsegs=apical_nsec*nsegs,
                basal_nsec_nsegs=basal_nsec*nsegs)

h("""objref apicalina[{apical_nsec_nsegs}], \
apicalik[{apical_nsec_nsegs}], \
apicalica[{apical_nsec_nsegs}], \
apicalih[{apical_nsec_nsegs}], \
apicalil[{apical_nsec_nsegs}], \
apicalv[{apical_nsec_nsegs}], \
apicalicap[{apical_nsec_nsegs}], \
apicalimemb[{apical_nsec_nsegs}]
objref basalih[{basal_nsec_nsegs}], \
basalil[{basal_nsec_nsegs}], \
basalv[{basal_nsec_nsegs}], \
basalicap[{basal_nsec_nsegs}], \
basalimemb[{basal_nsec_nsegs}]
objref somaticina[{nsegs}], \
somaticik[{nsegs}], \
somaticica[{nsegs}], \
somaticih[{nsegs}], \
somaticil[{nsegs}], \
somaticv[{nsegs}], \
somaticicap[{nsegs}], \
somaticimemb[{nsegs}]
objref axonalil[{axonal_nsec_nsegs}], \
axonalv[{axonal_nsec_nsegs}], \
axonalicap[{axonal_nsec_nsegs}], \
axonalimemb[{axonal_nsec_nsegs}]""".format(**tmp_dict))

# Initialize the segment areas
print("objref complete")
A_apical = zeros(apical_nsegs)
A_basal = zeros(basal_nsegs)
A_somatic = zeros(somatic_nsegs)
A_axonal = zeros(axonal_nsegs)

Asec_apical = zeros(apical_nsec)
x_apical = zeros(apical_nsec)
y_apical = zeros(apical_nsec)
z_apical = zeros(apical_nsec)

Asec_basal = zeros(basal_nsec)
x_basal = zeros(basal_nsec)
y_basal = zeros(basal_nsec)
z_basal = zeros(basal_nsec)

Asec_somatic = zeros(somatic_nsec)
y_somatic = zeros(somatic_nsec)
x_somatic = zeros(somatic_nsec)
z_somatic = zeros(somatic_nsec)

Asec_axonal = zeros(axonal_nsec)
x_axonal = zeros(axonal_nsec)
y_axonal = zeros(axonal_nsec)
z_axonal = zeros(axonal_nsec)

# Set the recordings for the compartments in the apical dendrite.
# Calculate also the membrane areas of each segment
for i in range(apical_nsec):
    h("access L5PC.apic[" + str(i) + "]")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarz = z3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    h("tmpvarz2 = z3d(n3d()-1)")
    coord1 = [h.tmpvarx, h.tmpvary, h.tmpvarz]
    coord2 = [h.tmpvarx2, h.tmpvary2, h.tmpvarz2]
    thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]
    x_apical[i] = thispos3d[0]
    y_apical[i] = thispos3d[1]
    z_apical[i] = thispos3d[2]

    for j in range(nsegs):
        thispos = 0.5/nsegs + 1.0/nsegs*j
        if nsegs == 1:
            thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]
        else:
            thispos3d = [x + j*(y-x)/(nsegs-1) for x, y in zip(coord1, coord2)]

        nsegs_str = str(i*nsegs + j)
        h("{apicalina[" + nsegs_str + "] = new Vector()}")
        h("{apicalik[" + nsegs_str + "] = new Vector()}")
        h("{apicalica[" + nsegs_str + "] = new Vector()}")
        h("{apicalih[" + nsegs_str + "] = new Vector()}")
        h("{apicalil[" + nsegs_str + "] = new Vector()}")
        h("{apicalv[" + nsegs_str + "] = new Vector()}")
        h("{apicalicap[" + nsegs_str + "] = new Vector()}")
        h("{apicalimemb[" + nsegs_str + "] = new Vector()}")
        h("L5PC.apic[" + str(i) + "] insert extracellular")
        h("access L5PC.apic[" + str(i) + "]")
        h("{apicalina[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].ina(" + str(thispos) + "),dtsave)}")
        h("{apicalik[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].ik(" + str(thispos) + "),dtsave)}")
        h("{apicalica[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].ica(" + str(thispos) + "),dtsave)}")
        h("{apicalih[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].ihcn_Ih(" + str(thispos) + "),dtsave)}")
        h("{apicalil[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].i_pas(" + str(thispos) + "),dtsave)}")
        h("{apicalv[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].v(" + str(thispos) + "),dtsave)}")
        h("{apicalicap[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].i_cap(" + str(thispos) + "),dtsave)}")
        h("{apicalimemb[" + nsegs_str + "].record(&L5PC.apic[" + str(i) +
            "].i_membrane(" + str(thispos) + "),dtsave)}")
        h("L5PC.apic[" + str(i) + "].nseg = " + str(nsegs))
        h("tmpvar = area(" + str(thispos) + ")")
        A_apical[i*nsegs + j] = h.tmpvar
        Asec_apical[i] += A_apical[i*nsegs + j]

print("apical complete")

# Set the recordings for the compartments in the basal dendrite.
# Calculate also the membrane areas of each segment
for i in range(basal_nsec):
    h("access L5PC.dend[" + str(i) + "]")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarz = z3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    h("tmpvarz2 = z3d(n3d()-1)")
    coord1 = [h.tmpvarx, h.tmpvary, h.tmpvarz]
    coord2 = [h.tmpvarx2, h.tmpvary2, h.tmpvarz2]

    thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]

    x_basal[i] = thispos3d[0]
    y_basal[i] = thispos3d[1]
    z_basal[i] = thispos3d[2]

    for j in range(nsegs):
        thispos = 0.5/nsegs + 1.0/nsegs*j
        if nsegs == 1:
            thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]
        else:
            thispos3d = [x + j*(y-x)/(nsegs-1) for x, y in zip(coord1, coord2)]

        nsegs_str = str(i*nsegs + j)
        h("{basalih[" + nsegs_str + "] = new Vector()}")
        h("{basalil[" + nsegs_str + "] = new Vector()}")
        h("{basalv[" + nsegs_str + "] = new Vector()}")
        h("{basalicap[" + nsegs_str + "] = new Vector()}")
        h("{basalimemb[" + nsegs_str + "] = new Vector()}")
        h("L5PC.dend[" + str(i) + "] insert extracellular")
        h("access L5PC.dend[" + str(i) + "]")
        h("{basalih[" + nsegs_str + "].record(&L5PC.dend[" + str(i) +
            "].ihcn_Ih(" + str(thispos) + "),dtsave)}")
        h("{basalil[" + nsegs_str + "].record(&L5PC.dend[" + str(i) +
            "].i_pas(" + str(thispos) + "),dtsave)}")
        h("{basalv[" + nsegs_str + "].record(&L5PC.dend[" + str(i) +
            "].v(" + str(thispos) + "),dtsave)}")
        h("{basalicap[" + nsegs_str + "].record(&L5PC.dend[" + str(i) +
            "].i_cap(" + str(thispos) + "),dtsave)}")
        h("{basalimemb[" + nsegs_str + "].record(&L5PC.dend[" + str(i) +
            "].i_membrane(" + str(thispos) + "),dtsave)}")
        h("L5PC.dend[" + str(i) + "].nseg = " + str(nsegs))
        h("tmpvar = area(" + str(thispos) + ")")
        A_basal[i*nsegs + j] = h.tmpvar
        Asec_basal[i] += A_basal[i*nsegs + j]

print("basal complete")

# Set the recordings for the compartments in the soma.
# Calculate also the membrane area
for i in range(somatic_nsec):
    h("access L5PC.soma[" + str(i) + "]")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarz = z3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    h("tmpvarz2 = z3d(n3d()-1)")
    coord1 = [h.tmpvarx, h.tmpvary, h.tmpvarz]
    coord2 = [h.tmpvarx2, h.tmpvary2, h.tmpvarz2]

    thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]
    x_somatic[i] = thispos3d[0]
    y_somatic[i] = thispos3d[1]
    z_somatic[i] = thispos3d[2]

    for j in range(nsegs):
        thispos = 0.5/nsegs + 1.0/nsegs*j
        if nsegs == 1:
            thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]
        else:
            thispos3d = [x + j*(y-x)/(nsegs-1) for x, y in zip(coord1, coord2)]

        nsegs_str = str(i*nsegs + j)
        h("{somaticina[" + nsegs_str + "] = new Vector()}")
        h("{somaticik[" + nsegs_str + "] = new Vector()}")
        h("{somaticica[" + nsegs_str + "] = new Vector()}")
        h("{somaticih[" + nsegs_str + "] = new Vector()}")
        h("{somaticil[" + nsegs_str + "] = new Vector()}")
        h("{somaticv[" + nsegs_str + "] = new Vector()}")
        h("{somaticicap[" + nsegs_str + "] = new Vector()}")
        h("{somaticimemb[" + nsegs_str + "] = new Vector()}")
        h("L5PC.soma[" + str(i) + "] insert extracellular")
        h("access L5PC.soma[" + str(i) + "]")
        h("{somaticina[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].ina(" + str(thispos) + "),dtsave)}")
        h("{somaticik[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].ik(" + str(thispos) + "),dtsave)}")
        h("{somaticica[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].ica(" + str(thispos) + "),dtsave)}")
        h("{somaticih[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].ihcn_Ih(" + str(thispos) + "),dtsave)}")
        h("{somaticil[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].i_pas(" + str(thispos) + "),dtsave)}")
        h("{somaticv[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].v(" + str(thispos) + "),dtsave)}")
        h("{somaticicap[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].i_cap(" + str(thispos) + "),dtsave)}")
        h("{somaticimemb[" + nsegs_str + "].record(&L5PC.soma[" + str(i) +
            "].i_membrane(" + str(thispos) + "),dtsave)}")
        h("L5PC.soma[" + str(i) + "].nseg = " + str(nsegs))
        h("tmpvar = area(" + str(thispos) + ")")
        A_somatic[i*nsegs + j] = h.tmpvar
        Asec_somatic[i] += A_somatic[i*nsegs + j]

        if j == nsegs/2:
            somapos3d = thispos3d[:]

print("somatic complete")

# Set the recordings for the compartments in the axon initial segment.
# Calculate also the membrane areas of each segment
for i in range(axonal_nsec):
    thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]
    x_axonal[i] = thispos3d[0]
    y_axonal[i] = thispos3d[1]
    z_axonal[i] = thispos3d[2]

    for j in range(nsegs):
        if nsegs == 1:
            thispos3d = [(x + y)/2 for x, y in zip(coord1, coord2)]
        else:
            thispos3d = [x + j*(y-x)/(nsegs-1) for x, y in zip(coord1, coord2)]

        nsegs_str = str(i*nsegs + j)
        h("{axonalil[" + nsegs_str + "] = new Vector()}")
        h("{axonalv[" + nsegs_str + "] = new Vector()}")
        h("{axonalicap[" + nsegs_str + "] = new Vector()}")
        h("{axonalimemb[" + nsegs_str + "] = new Vector()}")
        h("L5PC.axon[" + str(i) + "] insert extracellular")
        h("access L5PC.axon[" + str(i) + "]")
        h("{axonalil[" + nsegs_str + "].record(&L5PC.axon[" + str(i) +
            "].i_pas(" + str(thispos) + "),dtsave)}")
        h("{axonalv[" + nsegs_str + "].record(&L5PC.axon[" + str(i) +
            "].v(" + str(thispos) + "),dtsave)}")
        h("{axonalicap[" + nsegs_str + "].record(&L5PC.axon[" + str(i) +
            "].i_cap(" + str(thispos) + "),dtsave)}")
        h("{axonalimemb[" + nsegs_str + "].record(&L5PC.axon[" + str(i) +
            "].i_membrane(" + str(thispos) + "),dtsave)}")
        h("L5PC.axon[" + str(i) + "].nseg = " + str(nsegs))
        h("tmpvar = area(" + str(thispos) + ")")
        A_axonal[i*nsegs + j] = h.tmpvar
        Asec_axonal[i] += A_axonal[i*nsegs + j]

print("axonal complete")

synbranch = [0]*Nsynlocs
syniseg = [0]*Nsynlocs
synx = [0.0]*Nsynlocs

# Calculate the probability of synapse being found in the basal dendrite.
if synloctype == 1:
    basalprob = 0.0
if synloctype == 2:
    basalprob = sum(A_basal)/(sum(A_basal) + sum(A_apical))
if synloctype == 3:
    basalprob = 1.0

print("Basal area:", str(sum(A_basal)))
print("Apical area:", str(sum(A_apical)))
print("basalprob =", str(basalprob))

# Calculate the probabilities for the synapse being in each segment
ps_basal = [1.0*x/sum(A_basal) for x in A_basal]
cumps_basal = cumsum(ps_basal)
ps_apical = [1.0*x/sum(A_apical) for x in A_apical]
cumps_apical = cumsum(ps_basal)

# Draw the random numbers, one to decide which branch, one to decide which
# segment, and one to determine the distance x from 0-end
rs_branch = rand(Nsynlocs)
rs_seg = rand(Nsynlocs)
rs_x = rand(Nsynlocs)

# print(Nsynlocs)
# asdas

ts_syn = []
seg_syn = [-1]*Nsynlocs
seg_syn_accurate = [-1]*Nsynlocs

# For each synapse, determine to which section it outputs the currents, and
# randomize the synapse activation times
# TODO: Might become faster if the set of AlphaSynapses at a single synaptic
#       location is replaced by a single point process that is activated at the
#       time instants drawn here.
allseg_syn = [-1]*Nsynlocs

for isyn in range(Nsynlocs):
    if rs_branch[isyn] <= basalprob:
        synbranch[isyn] = 1

    if synbranch[isyn] == 0:
        seg_syn_accurate[isyn] = next((i for i, x in enumerate(cumps_apical)
                                       if x > rs_seg[isyn]))
        seg_syn[isyn] = seg_syn_accurate[isyn]/nsegs
        segnum = seg_syn_accurate[isyn] % nsegs
        h("access L5PC.apic[" + str(seg_syn[isyn]) + "]")
        mystr = "L5PC.apic[" + str(seg_syn[isyn]) + "]"
        allseg_syn[isyn] = seg_syn_accurate[isyn]

    else:
        seg_syn_accurate[isyn] = next((i for i, x in enumerate(cumps_basal)
                                       if x > rs_seg[isyn]))
        #  seg_syn_accurate[isyn]
        seg_syn[isyn] = seg_syn_accurate[isyn]/nsegs
        segnum = seg_syn_accurate[isyn] % nsegs
        h("access L5PC.dend[" + str(seg_syn[isyn]) + "]")
        mystr = "L5PC.dend[" + str(seg_syn[isyn]) + "]"
        allseg_syn[isyn] = seg_syn_accurate[isyn] + apical_nsegs

    ts = []
    t = 0
    secx = 1.0*segnum/nsegs + 1.0/nsegs*rs_x[isyn]
    while t < tstop:
        t = t - 1000.0/synlambda*log(1-rand())
        ts.append(t)
        h("{synlist.append(new AlphaSynapse(" + str(secx) + "))}")
        h("syni = synlist.count()-1")
        h("synlist.o[syni].tau = " + str(syntau))
        h("synlist.o[syni].gmax = " + str(syngmax))
        h("synlist.o[syni].e = 0")
        h("synlist.o[syni].onset = " + str(t))
    ts_syn.append(ts[:])

Nsyns = h.syni + 1

h("""
tstop = """ + str(tstop) + """
v_init = """ + str(v0) + """
cai0_ca_ion = """ + str(ca0) + """
st1.amp = 0
st1.dur = 0
""")

print("Initializing...")
h("stdinit()")
print("Init complete")

tfin = 0

for isim in range(Nsims):
    h("{tvec.resize(0)}")
    h("{vsoma.resize(0)}")
    h("{vdend.resize(0)}")
    h("{casoma.resize(0)}")
    h("{cadend.resize(0)}")

    for i in range(apical_nsec*nsegs):
        h("{apicalina[" + str(i) + "].resize(0)}")
        h("{apicalik[" + str(i) + "].resize(0)}")
        h("{apicalica[" + str(i) + "].resize(0)}")
        h("{apicalih[" + str(i) + "].resize(0)}")
        h("{apicalil[" + str(i) + "].resize(0)}")
        h("{apicalv[" + str(i) + "].resize(0)}")
        h("{apicalicap[" + str(i) + "].resize(0)}")
        h("{apicalimemb[" + str(i) + "].resize(0)}")

    for i in range(basal_nsec*nsegs):
        h("{basalih[" + str(i) + "].resize(0)}")
        h("{basalil[" + str(i) + "].resize(0)}")
        h("{basalv[" + str(i) + "].resize(0)}")
        h("{basalicap[" + str(i) + "].resize(0)}")
        h("{basalimemb[" + str(i) + "].resize(0)}")

    for i in range(nsegs):
        h("{somaticina[" + str(i) + "].resize(0)}")
        h("{somaticik[" + str(i) + "].resize(0)}")
        h("{somaticica[" + str(i) + "].resize(0)}")
        h("{somaticih[" + str(i) + "].resize(0)}")
        h("{somaticil[" + str(i) + "].resize(0)}")
        h("{somaticv[" + str(i) + "].resize(0)}")
        h("{somaticicap[" + str(i) + "].resize(0)}")
        h("{somaticimemb[" + str(i) + "].resize(0)}")

    for i in range(axonal_nsec*nsegs):
        h("{axonalil[" + str(i) + "].resize(0)}")
        h("{axonalv[" + str(i) + "].resize(0)}")
        h("{axonalicap[" + str(i) + "].resize(0)}")
        h("{axonalimemb[" + str(i) + "].resize(0)}")

    tfin = tfin + singleSimT
    print("Starting run " + str(isim) + " until " + str(tfin) + " ms")
    h("continuerun(" + str(tfin) + ")")
    print("Run " + str(isim) + " complete, tfin = " + str(tfin))

    Vsoma = np.array(h.vsoma)
    Vdend = np.array(h.vdend)
    Casoma = np.array(h.casoma)
    Cadend = np.array(h.cadend)
    times = np.array([tfin-singleSimT + dtsave*i for i in range(len(Vsoma))])

    ina = zeros([all_nsec, len(times)])
    ik = zeros([all_nsec, len(times)])
    ica = zeros([all_nsec, len(times)])
    ih = zeros([all_nsec, len(times)])
    il = zeros([all_nsec, len(times)])
    v = zeros([all_nsec, len(times)])
    icap = zeros([all_nsec, len(times)])
    imemb = zeros([all_nsec, len(times)])
    A = zeros(all_nsec)
    x = zeros(all_nsec)
    y = zeros(all_nsec)
    z = zeros(all_nsec)

    a = np.array(h.apicalina[i])
    k = 0
    for i in range(apical_nsec):
        for j in range(nsegs):
            ina[k, :] += np.array(h.apicalina[i*nsegs + j])*A_apical[i*nsegs +
                                                                     j]
            ik[k, :] += np.array(h.apicalik[i*nsegs + j])*A_apical[i*nsegs + j]
            ica[k, :] += np.array(h.apicalica[i*nsegs + j])*A_apical[i*nsegs +
                                                                     j]
            ih[k, :] += np.array(h.apicalih[i*nsegs + j])*A_apical[i*nsegs + j]
            il[k, :] += np.array(h.apicalil[i*nsegs + j])*A_apical[i*nsegs + j]
            v[k, :] += np.array(h.apicalv[i*nsegs + j])*A_apical[i*nsegs + j]
            icap[k, :] += np.array(h.apicalicap[i*nsegs +
                                                j])*A_apical[i*nsegs + j]
            imemb[k, :] += np.array(h.apicalimemb[i*nsegs +
                                                  j])*A_apical[i*nsegs + j]
        A[k] = A_apical[i]
        x[k] = x_apical[i]
        y[k] = y_apical[i]
        z[k] = z_apical[i]
        k += 1

    for i in range(basal_nsec):
        for j in range(nsegs):
            ih[k, :] += np.array(h.basalih[i*nsegs + j])*A_basal[i*nsegs + j]
            il[k, :] += np.array(h.basalil[i*nsegs + j])*A_basal[i*nsegs + j]
            v[k, :] += np.array(h.basalv[i*nsegs + j])*A_basal[i*nsegs + j]
            icap[k, :] += np.array(h.basalicap[i*nsegs +
                                               j])*A_basal[i*nsegs + j]
            imemb[k, :] += np.array(h.basalimemb[i*nsegs +
                                                 j])*A_basal[i*nsegs + j]
        A[k] = A_basal[i]
        x[k] = x_basal[i]
        y[k] = y_basal[i]
        z[k] = z_basal[i]
        k += 1

    for i in range(somatic_nsec):
        for j in range(nsegs):
            ina[k, :] += np.array(h.somaticina[i*nsegs +
                                               j])*A_somatic[i*nsegs + j]
            ik[k, :] += np.array(h.somaticik[i*nsegs +
                                             j])*A_somatic[i*nsegs + j]
            ica[k, :] += np.array(h.somaticica[i*nsegs +
                                               j])*A_somatic[i*nsegs + j]
            ih[k, :] += np.array(h.somaticih[i*nsegs +
                                             j])*A_somatic[i*nsegs + j]
            il[k, :] += np.array(h.somaticil[i*nsegs +
                                             j])*A_somatic[i*nsegs + j]
            v[k, :] += np.array(h.somaticv[i*nsegs +
                                           j])*A_somatic[i*nsegs + j]
            icap[k, :] += np.array(h.somaticicap[i*nsegs +
                                                 j])*A_somatic[i*nsegs + j]
            imemb[k, :] += np.array(h.somaticimemb[i*nsegs +
                                                   j])*A_somatic[i*nsegs + j]
        A[k] = A_somatic[i]
        x[k] = x_somatic[i]
        y[k] = y_somatic[i]
        z[k] = z_somatic[i]
        k += 1

    for i in range(axonal_nsec):
        for j in range(nsegs):
            il[k, :] += np.array(h.axonalil[i*nsegs + j])*A_axonal[i*nsegs + j]
            v[k, :] += np.array(h.axonalv[i*nsegs + j])*A_axonal[i*nsegs + j]
            icap[k, :] += np.array(h.axonalicap[i*nsegs +
                                                j])*A_axonal[i*nsegs + j]
            imemb[k, :] += np.array(h.axonalimemb[i*nsegs +
                                                  j])*A_axonal[i*nsegs + j]
        A[k] = A_axonal[i]
        x[k] = x_axonal[i]
        y[k] = y_axonal[i]
        z[k] = z_axonal[i]
        k += 1

    # divide all currents by 100 to get currents in nA
    ina = ina/100
    ik = ik/100
    ica = ica/100
    ih = ih/100
    il = il/100
    icap = icap/100
    imemb = imemb/100

    # move soma to origin:
    h("access L5PC.soma[0]")
    h("tmpvarx = x3d(0)")
    h("tmpvary = y3d(0)")
    h("tmpvarz = z3d(0)")
    h("tmpvarx2 = x3d(n3d()-1)")
    h("tmpvary2 = y3d(n3d()-1)")
    h("tmpvarz2 = z3d(n3d()-1)")

    x_offset = (h.tmpvarx + h.tmpvarx2)/2.
    y_offset = (h.tmpvary + h.tmpvary2)/2.
    z_offset = (h.tmpvarz + h.tmpvarz2)/2.
    x = x - x_offset
    y = y - y_offset
    z = z - z_offset

    dict_ = {'ina': ina, 'ik': ik, 'ica': ica, 'ih': ih, 'il': il,
             'VtimesA': v, 'imemb': imemb, 'Vsoma': np.array(Vsoma),
             'icap': icap, 'times': np.array(times), 'A': A, 'x': x, 'y': y,
             'z': z, 'ts_syn': np.array(ts_syn),
             'sec_syn': np.array(allseg_syn)}

    scipy.io.savemat(path.join('simulation_data', 'seed_' + str(myseed),
                               'currsums_parts_' + str(Nsynlocs) +
                               'areagsynsmediumtau_fixeddt_type' +
                               str(synloctype) + '_amp'+str(syngmax) + '_tstop'
                               + str(tstop) + '_nseg' + str(nsegs) + '_dt' +
                               str(dt) + '_seed' + str(myseed) + '_sim' +
                               str(isim) + 'x' + str(singleSimT) + '.mat'),
                     dict_)
