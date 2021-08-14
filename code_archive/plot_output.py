import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pdb
from mpl_toolkits.mplot3d import Axes3D
import h5py
import matplotlib.pyplot as plt
from scipy.signal import welch
import scipy.signal as ss

tsim = 5000

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def zscore(x):
    return (x - np.mean(x))/np.std(x)

# Load data
config_file = "simulation_configECP.json"
lfp_file = "./outputECP/ecp.h5"
mem_pot_file = './outputECP/v_report.h5'
raster_file = './outputECP/spikes.h5'
conns_file = './network/BLA_BLA_edges.h5'


# load 
#f = h5py.File(mem_pot_file,'r')
#mem_potential = f['report']['exc_bg_bask']['data']

f = h5py.File(lfp_file,'r')
lfp = list(f['ecp']['data'])
lfp_arr = np.asarray(lfp)

f = h5py.File(raster_file,'r')
#gids = f['spikes']['exc_bg_bask']['node_ids']
#timestamps = f['spikes']['exc_bg_bask']['timestamps']
gids = f['spikes']['BLA']['node_ids']
timestamps = f['spikes']['BLA']['timestamps']

plt.figure()
plt.title('LFP')
lfp1 = zscore(lfp_arr[:,0])
#lfp2 = zscore(lfp_arr[:,2])
#lfp3 = zscore(lfp_arr[:,4])
#np.savetxt("lfp_ben.csv", lfp1, delimiter=",")
plt.plot(np.arange(0,tsim,0.1),lfp1)
#plt.plot(np.arange(0,tsim,0.1),lfp2)
#plt.plot(np.arange(0,tsim,0.1),lfp3)
plt.xlim(4000,5000)

freqs, psd = welch(lfp1,fs=10000)


plt.figure()
plt.semilogy(freqs, psd)
plt.title('LFP PSD')
plt.xlabel('Frequency')
plt.ylabel('Power')


[x,n] = np.unique(gids,return_counts=True)
plt.figure()
plt.title('INT FR')
plt.hist(n/(tsim/1000))

"""
gids2 = f['spikes']['mthalamus']['node_ids']
timestamps2 = f['spikes']['mthalamus']['timestamps']


[x,n] = np.unique(gids2,return_counts=True)
plt.figure()
plt.title('PN FR')
plt.hist(n/(tsim/1000))
"""

#f = h5py.File(conns_file,'r')
#src = f['edges']['SPWR_biophysical_SPWR_biophysical']['source_node_id']
#tgt = f['edges']['SPWR_biophysical_SPWR_biophysical']['target_node_id']
#wgt = f['edges']['SPWR_biophysical_SPWR_biophysical']['0']['syn_weight']
#
#id, pyr2baskconns = np.unique(src[(src[:]<=pyr[1]) & (tgt[:]>=bask[0]) \
#			& (tgt[:]<=bask[1])],return_counts=True)
#
#
#pyr2baskwgts = wgt[(src[:]<=pyr[1]) & (tgt[:]>=bask[0]) \
#			& (tgt[:]<=bask[1])]
#
#
#plt.figure()
#plt.subplot(2,1,1)
#plt.hist(pyr2baskconns,20)
#plt.title('pyr2baskconns')
#plt.subplot(2,1,2)
#plt.hist(pyr2baskwgts,20)
#plt.title('pyr2baskwgts')


# Plot data

#plt.figure()
#plt.plot(mem_potential[:,0])
#plt.plot(mem_potential[:,1])
#plt.plot(mem_potential[:,2])
#plt.plot(mem_potential[:,3])

#plt.figure()
#plt.plot(lfp_arr[:,0])
#plt.plot(lfp_arr[:,1])
#plt.plot(lfp_arr[:,2])

"""
plt.figure()
plt.plot(timestamps,gids,'r.')

plt.plot(timestamps2,gids2 + np.max(gids),'b.')
"""

#plt.figure()
#x = lfp_arr[:,0]
#fs = 10000
#f, Pxx_den = welch(x, fs, nperseg=2000)
#plt.semilogy(f, Pxx_den)
#plt.xlim(0,500)

plt.figure()
f, t, Sxx = ss.spectrogram(lfp_arr[:,0], fs=10000)
plt.pcolormesh(t, f, Sxx, shading='gouraud')
plt.xlabel('time')
plt.ylabel('freq')
plt.show()

