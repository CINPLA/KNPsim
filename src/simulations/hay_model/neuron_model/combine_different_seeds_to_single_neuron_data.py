import scipy.io
import numpy as np
import sys
import os

iseed = int(sys.argv[1])

list_ica = [];
list_ina = [];
list_ik = [];
list_ix = [];
list_icap = [];
list_t = [];
list_imemb = [];
list_isyn = [];

# 'simulation_data/currsums_parts_10000areagsynsmediumtau_fixeddt_type2_amp4.2e-05_tstop10000.0_nseg20_dt0.025_seed1_com200.0.mat'
# 'simulation_data/currsums_parts_10000areagsynsmediumtau_fixeddt_type2_amp4.2e-05_tstop10000.0_nseg20_dt0.025_seed1_comb200.0'

for myseed in range(10):
    strseed = str(myseed*10 + iseed)
    filename = 'simulation_data/currsums_parts_10000areagsynsmediumtau_fixeddt_type2_amp4.2e-05_tstop10000.0_nseg20_dt0.025_seed' + strseed + '_comb200.0.mat'
    A = scipy.io.loadmat(filename)

    iion0 = A['ik'] + A['ina'] + A['ica'] + A['il'] + A['ih'];
    isyn = A['imemb'] - (iion0 + A['icap']);

    ix = A['il'] + A['ih'] + isyn; # all currents of unspecified ion species x
    ts_sim = A['times']
    ts_sim = np.reshape(ts_sim, np.size(ts_sim))

    # Remove first 1600 ms of each simulation, to give the neuron time to settle
    ts_sim2 = ts_sim[ts_sim >= 1600] + 8400*myseed

    ica = A['ica']
    ica = ica[:,np.where(ts_sim >= 1600)]
    ica = np.reshape(ica, (ica.shape[0], ica.shape[2]))

    ina = A['ina']
    ina = ina[:,np.where(ts_sim >= 1600)]
    ina = np.reshape(ina, (ina.shape[0], ina.shape[2]))

    ik = A['ik']
    ik = ik[:,np.where(ts_sim >= 1600)]
    ik = np.reshape(ik, (ik.shape[0], ik.shape[2]))

    ix = ix[:,np.where(ts_sim >= 1600)]
    ix = np.reshape(ix, (ix.shape[0], ix.shape[2]))

    icap = A['icap']
    icap = icap[:,np.where(ts_sim >= 1600)]
    icap = np.reshape(icap, (icap.shape[0], icap.shape[2]))

    imemb = A['imemb']
    imemb = imemb[:,np.where(ts_sim >= 1600)]
    imemb = np.reshape(imemb, (imemb.shape[0], imemb.shape[2]))

    isyn = isyn[:,np.where(ts_sim >= 1600)]
    isyn = np.reshape(isyn, (isyn.shape[0], isyn.shape[2]))

    list_t.append(ts_sim2)
    list_ica.append(ica)
    list_ina.append(ina)
    list_ik.append(ik)
    list_ix.append(ix)
    list_icap.append(icap)
    list_imemb.append(imemb)
    list_isyn.append(isyn)

    if myseed == 0:
        x = A['x']
        y = A['y']
        z = A['z']

times = np.concatenate(list_t)
ica = np.concatenate(list_ica, axis=1)
ina = np.concatenate(list_ina, axis=1)
ik = np.concatenate(list_ik, axis=1)
ix = np.concatenate(list_ix, axis=1)
icap = np.concatenate(list_icap, axis=1)
imemb = np.concatenate(list_imemb, axis=1)
isyn = np.concatenate(list_isyn, axis=1)

# Convert to base SI units:
times = times*1e-3  # s
ica = ica*1e-9      # nA
ina = ina*1e-9      # nA
ik = ik*1e-9        # nA
ix = ix*1e-9        # nA
icap = icap*1e-9    # nA
imemb = imemb*1e-9  # nA
isyn = isyn*1e-9    # nA

x = x*1e-6  # meters
y = y*1e-6  # meters
z = z*1e-6  # meters

B = {
    'times' : times,
    'ina'   : ina,
    'ik'    : ik,
    'ica'   : ica,
    'ix'    : ix,
    'imemb' : imemb,
    'icap'  : icap,
    'x' : x,
    'y' : y,
    'z' : z,
}

folder_name = 'simulation_data/neuron_' + str(iseed)
filename = folder_name + '/revdata_neuron_' + str(iseed) + '.mat'
os.mkdir(folder_name)
scipy.io.savemat(filename, B)
