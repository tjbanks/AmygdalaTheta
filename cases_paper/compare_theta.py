"""
Example Use (with one simple command [ha]):

python compare_theta.py --first-case case1 --first-case-runs run1,run2,run3,run4,run5 --first-case-label "Base Case" --second-case case_no_som --second-case-runs run1 --second-case-label "No SOM+"

python compare_theta.py --first-case case1 --first-case-runs run1 --first-case-label "Quiet Waking" --second-case case2 --second-case-runs run1 --second-case-label VPSI --third-case case3 --third-case-runs run1 --third-case-label "ACH High" --fourth-case case4 --fourth-case-runs run1 --fourth-case-label "ACH Medium" --fifth-case case5 --fifth-case-runs run1 --fifth-case-label "ACH High VPSI Off" --sixth-case case6 --sixth-case-runs run1 --sixth-case-label "ACH Medium VPSI Off"
"""

import argparse
import scipy
from scipy.signal.windows import hann as hanning
from scipy.signal import welch,decimate
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

def ecp_psd(first_ecps, second_ecps, ax=None,label1='',label2='',skip_n=0,third_ecps=None,label3='',fourth_ecps=None,label4='',fifth_ecps=None,label5='',sixth_ecps=None,label6=''):

    f1,pxx1 = mean_ecp(first_ecps,skip_n)    
    f2,pxx2 = mean_ecp(second_ecps,skip_n)

    fx1 = f1[np.where((f1>8) & (f1<12))]
    theta1 = pxx1[np.where((f1>8) & (f1<12))]
    fx2 = f2[np.where((f2>8) & (f2<12))]
    theta2 = pxx2[np.where((f2>8) & (f2<12))]

    label1 = label1 + " ("+ str(round(theta1.mean()*1000,6))+")"
    label2 = label2 + " ("+ str(round(theta2.mean()*1000,6))+")"
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(fx1, theta1*1000,linewidth=0.6,label=label1)
    ax.plot(fx2, theta2*1000,linewidth=0.6,label=label2)
    #ax.set_ylim([0,0.1])
    f_arr = [f1, f2]
    pxx_arr = [pxx1,pxx2]
    label_arr = [label1,label2]

    if third_ecps:
        f3,pxx3 = mean_ecp(third_ecps,skip_n)
        fx3 = f3[np.where((f3>8) & (f3<12))]
        theta3 = pxx3[np.where((f3>8) & (f3<12))]
        label3 = label3 + " ("+ str(round(theta3.mean()*1000,6))+")"
        ax.plot(fx3, theta3*1000,linewidth=0.6,label=label3)
        f_arr.append(f3)
        pxx_arr.append(pxx3)
        label_arr.append(label3)

    if fourth_ecps:
        f4,pxx4 = mean_ecp(fourth_ecps,skip_n)
        fx4 = f4[np.where((f4>8) & (f4<12))]
        theta4 = pxx4[np.where((f4>8) & (f4<12))]
        label4 = label4 + " ("+ str(round(theta4.mean()*1000,6))+")"
        ax.plot(fx4, theta4*1000,linewidth=0.6,label=label4)
        f_arr.append(f4)
        pxx_arr.append(pxx4)
        label_arr.append(label4)

    if fifth_ecps:
        f5,pxx5 = mean_ecp(fifth_ecps,skip_n)
        fx5 = f5[np.where((f5>8) & (f5<12))]
        theta5 = pxx5[np.where((f5>8) & (f5<12))]
        label5 = label5 + " ("+ str(round(theta5.mean()*1000,6))+")"
        ax.plot(fx5, theta5*1000,linewidth=0.6,label=label5)
        f_arr.append(f5)
        pxx_arr.append(pxx5)
        label_arr.append(label5)

    if sixth_ecps:
        f6,pxx6 = mean_ecp(sixth_ecps,skip_n)
        fx6 = f6[np.where((f6>8) & (f6<12))]
        theta6 = pxx6[np.where((f6>8) & (f6<12))]
        label6 = label6 + " ("+ str(round(theta6.mean()*1000,6))+")"
        ax.plot(fx6, theta6*1000,linewidth=0.6,label=label6)
        f_arr.append(f6)
        pxx_arr.append(pxx6)
        label_arr.append(label6)

    for f,pxx,label in zip(f_arr,pxx_arr,label_arr):

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
    
    parser.add_argument('--third-case')
    parser.add_argument('--third-case-runs')
    parser.add_argument('--third-case-label')
    
    parser.add_argument('--fourth-case')
    parser.add_argument('--fourth-case-runs')
    parser.add_argument('--fourth-case-label')

    parser.add_argument('--fifth-case')
    parser.add_argument('--fifth-case-runs')
    parser.add_argument('--fifth-case-label')

    parser.add_argument('--sixth-case')
    parser.add_argument('--sixth-case-runs')
    parser.add_argument('--sixth-case-label')

    args = parser.parse_args()

    first_run_ecps = []
    second_run_ecps = []
    third_run_ecps = []
    fourth_run_ecps = []
    fifth_run_ecps = []
    sixth_run_ecps = []

    for run in args.first_case_runs.split(','):
        ecp = get_ecp(args.first_case + '/' + run + '/ecp.h5')
        first_run_ecps.append(ecp)

    for run in args.second_case_runs.split(','):
        ecp = get_ecp(args.second_case + '/' + run + '/ecp.h5')
        second_run_ecps.append(ecp)

    if args.third_case_runs:
        for run in args.third_case_runs.split(','):
            ecp = get_ecp(args.third_case + '/' + run + '/ecp.h5')
            third_run_ecps.append(ecp)

    if args.fourth_case_runs:
        for run in args.fourth_case_runs.split(','):
            ecp = get_ecp(args.fourth_case + '/' + run + '/ecp.h5')
            fourth_run_ecps.append(ecp)

    if args.fifth_case_runs:
        for run in args.fifth_case_runs.split(','):
            ecp = get_ecp(args.fifth_case + '/' + run + '/ecp.h5')
            fifth_run_ecps.append(ecp)

    if args.sixth_case_runs:
        for run in args.sixth_case_runs.split(','):
            ecp = get_ecp(args.sixth_case + '/' + run + '/ecp.h5')
            sixth_run_ecps.append(ecp)

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
    ecp_psd(first_run_ecps, second_run_ecps, ax=ax1, label1=args.first_case_label, label2=args.second_case_label,skip_n=skip_n, third_ecps=third_run_ecps, label3=args.third_case_label, fourth_ecps=fourth_run_ecps, label4=args.fourth_case_label, fifth_ecps=fifth_run_ecps, label5=args.fifth_case_label, sixth_ecps=sixth_run_ecps, label6=args.sixth_case_label)
    ax1.legend()
    plt.show()
