"""
A python implementation of matlab/analysis.m

TB - 8/4/21
"""

from scipy.signal import hanning,welch,decimate
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

def raster():
    pass

def raw_ecp(lfp):
    pass

def ecp_psd(ecp_h5,channel=0,skip_n=0,downsample=20,nfft=1024,fs=1000,noverlap=0):
    f = h5py.File(ecp_h5)
    data_raw = np.array(f['ecp']['data'])
    data = data_raw.T[channel] #flip verts and grab channel 0
    
    #skip_n first few
    data = data[skip_n:]

    #downsample the data to fit ms (steps used 20=1/.05 step)
    lfp_d = decimate(data,downsample)
    raw_ecp(lfp_d)
    win = hanning(nfft, True)

    f,pxx = welch(lfp_d,fs,window=win,noverlap=noverlap,nfft=nfft)
    
    plt.figure()
    #plt.yscale('symlog')
    #plt.xscale('log')
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(f, pxx*1000,linewidth=0.6)
    axes = plt.gca()
    axes.set_ylim([0.00001,0.1])
    plt.show()
    
    

def spike_frequency_histogram():
    pass


def run(config):
    
    dt = 0.05
    steps_per_ms = 1/dt
    skip_seconds = 5
    skip_ms = skip_seconds*1000
    skip_n = int(skip_ms * steps_per_ms)
    end_ms = 15000

    ecp_h5 = 'outputECP/ecp.h5'
    ecp_psd(ecp_h5, skip_n=skip_n)


if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('simulation_config.json')
