"""
Example Use:

python compare_theta.py --first-case case1 --first-case-runs run1,run2,run3,run4,run5 --first-case-label "Base Case" --second-case case_no_som --second-case-runs run1 --second-case-label "No SOM+"
"""

import argparse
import scipy
from scipy.signal import hanning,welch,decimate
import h5py
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def mean_ecp(ecps,skip_n=0,downsample=20,nfft=1024,fs=1000,noverlap=0):
    # Average all the runs
    #ecp = np.mean(np.array(ecps),axis=0)
    #skip_n first few
    pxxs = []
    for ecp in ecps:
        data = ecp[skip_n:]

        #downsample the data to fit ms (steps used 20=1/.05 step)
        lfp_d = decimate(data,downsample)
        win = hanning(nfft, True)

        f,pxx = welch(lfp_d,fs,window=win,noverlap=noverlap,nfft=nfft)

        pxxs.append(pxx)

    pxx = np.mean(np.array(pxxs),axis=0)

    return f,pxx

def ecp_psd(first_ecps, second_ecps, ax=None,label1='',label2='',skip_n=0):

    f1,pxx1 = mean_ecp(first_ecps,skip_n)    
    f2,pxx2 = mean_ecp(second_ecps,skip_n)

    fx1 = f1[np.where((f1>8) & (f1<12))]
    theta1 = pxx1[np.where((f1>8) & (f1<12))]
    fx2 = f2[np.where((f2>8) & (f2<12))]
    theta2 = pxx2[np.where((f2>8) & (f2<12))]

    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(fx1, theta1*1000,linewidth=0.6,label=label1)
    ax.plot(fx2, theta2*1000,linewidth=0.6,label=label2)
    #ax.set_ylim([0,0.1])

    for f,pxx,label in zip([f1, f2],[pxx1,pxx2],[label1,label2]):

        theta = pxx[np.where((f>8) & (f<12))]*1000
        gamma = pxx[np.where((f>50) & (f<60))]*1000
        mean_theta = theta.mean()
        peak_theta = theta.max() 
        mean_gamma = gamma.mean()
        peak_gamma = gamma.max()
        print('')
        print(label + " Mean theta (8Hz-12Hz)  : " + str(round(mean_theta,6)))
        #print("Mean gamma (50Hz-60Hz) : " + str(round(mean_gamma,6)))     
        #print('')
        print(label + " Peak theta (8Hz-12Hz)  : " + str(round(peak_theta,6)))
        #print("Peak gamma (50Hz-60Hz) : " + str(round(peak_gamma,6)))
        print('')

def get_ecp(ecp_h5_location):
     
    print("loading " + ecp_h5_location)
    ecp_channel = 0
    f = h5py.File(ecp_h5_location)
    data_raw = np.array(f['ecp']['data'])
    ecp = data_raw.T[ecp_channel] #flip verts and grab channel 0
    return ecp

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--first-case')
    parser.add_argument('--first-case-runs')
    parser.add_argument('--first-case-label')
    parser.add_argument('--second-case')
    parser.add_argument('--second-case-runs')
    parser.add_argument('--second-case-label')

    args = parser.parse_args()

    first_run_ecps = []
    second_run_ecps = []

    for run in args.first_case_runs.split(','):
        ecp = get_ecp(args.first_case + '/' + run + '/ecp.h5')
        first_run_ecps.append(ecp)

    for run in args.second_case_runs.split(','):
        ecp = get_ecp(args.second_case + '/' + run + '/ecp.h5')
        second_run_ecps.append(ecp)


    dt = 0.05
    steps_per_ms = 1/dt
    skip_seconds = 5
    skip_ms = skip_seconds*1000
    skip_n = int(skip_ms * steps_per_ms)
    end_ms = 15000

    fig, (ax1) = plt.subplots(1,1)#,figsize=(15,4.8))#6.4,4.8 default
    #fig.suptitle('Amygdala Analysis')    
    
    # AX1 - Base raster
    ax1.set_title("Theta Comparison")
    ax1.set_xlabel("Hz")
    ax1.set_ylabel("PSD [V^2/Hz]")
    print('plotting...') 
    ecp_psd(first_run_ecps, second_run_ecps, ax=ax1, label1=args.first_case_label, label2=args.second_case_label,skip_n=skip_n)
    ax1.legend()
    plt.show()
