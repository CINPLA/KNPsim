import scipy.io
import numpy as np
import sys

myseed = int(sys.argv[1])

synloctype = 2
syngmax = 0.000042
nsegs = 20
dt = 0.025
tstop = 10000
Nsynlocs = 10000
singleSimT = 200

Nsims = int(1.0*tstop/singleSimT+0.9999)

ina = []
ik = []
ica = []
ih = []
il = []
VtimesA = []
imemb = []
Vsoma = []
icap = []
times = []

for isim in range(Nsims):
    print('isim=', isim)
    print('Loading isim= ', isim, ', myseed=', myseed)
    A = scipy.io.loadmat('simulation_data/seed_' +
                         str(myseed) +
                         '/currsums_parts_' +
                         str(Nsynlocs) +
                         'areagsynsmediumtau_fixeddt_type' +
                         str(synloctype) +
                         '_amp' +
                         str(syngmax) +
                         '_tstop'+str(tstop) +
                         '_nseg' +
                         str(nsegs) +
                         '_dt' +
                         str(dt) +
                         '_seed' +
                         str(myseed) +
                         '_sim' +
                         str(isim) +
                         'x' +
                         str(singleSimT) +
                         '.mat')

    print A.keys()
    times.append(A['times'])
    ina.append(A['ina'])
    ik.append(A['ik'])
    ica.append(A['ica'])
    ih.append(A['ih'])
    il.append(A['il'])
    VtimesA.append(A['VtimesA'])
    imemb.append(A['imemb'])
    icap.append(A['icap'])
    Vsoma.append(A['Vsoma'])

    if isim == 0:
        ts_syn = A['ts_syn']
        sec_syn = A['sec_syn']
        x = A['x']
        y = A['y']
        z = A['z']
        area = A['A']

times = np.concatenate(times, axis=1)
times = times.reshape(np.size(times))

ina = np.concatenate(ina, axis=1)
ik = np.concatenate(ik, axis=1)
ica = np.concatenate(ica, axis=1)
ih = np.concatenate(ih, axis=1)
il = np.concatenate(il, axis=1)
VtimesA = np.concatenate(VtimesA, axis=1)
imemb = np.concatenate(imemb, axis=1)
icap = np.concatenate(icap, axis=1)
Vsoma = np.concatenate(Vsoma, axis=1)

print(ina.shape)
print(ik.shap)
print(ica.shape)
print(ih.shape)
print(il.shape)
print(VtimesA.shape)
print(imemb.shape)
print(icap.shape)
print(Vsoma.shape)

B = {
    'times': times,
    'ina': ina,
    'ik': ik,
    'ica': ica,
    'ih': ih,
    'il': il,
    'VtimesA': VtimesA,
    'imemb': imemb,
    'icap': icap,
    'Vsoma': Vsoma,
    'ts_syn': ts_syn,
    'sec_syn': sec_syn,
    'x': x,
    'y': y,
    'z': z,
    'area': area
}

scipy.io.savemat('simulation_data/currsums_parts_' + str(Nsynlocs) +
                 'areagsynsmediumtau_fixeddt_type' + str(synloctype) +
                 '_amp' + str(syngmax) + '_tstop' + str(tstop) +
                 '.0_nseg' + str(nsegs) + '_dt' + str(dt) + '_seed' +
                 str(myseed) + '_comb200.0.mat', B)
