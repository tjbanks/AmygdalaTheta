"""
A python implementation of matlab/analysis.m

TB - 8/4/21
"""
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

scale = 1
spikes = {}
spike_hist = {}
psds = {}

node_set = [
        {"name":"PN","start":0*scale,"end":799*scale,"color":"blue"},
        {"name":"PV","start":800*scale,"end":892*scale,"color":"red"},
        {"name":"SOM","start":893*scale,"end":943*scale,"color":"green"},
        {"name":"CR","start":944*scale,"end":999*scale,"color":"purple"}
    ]

def raster(spikes_df,node_set,skip_ms=0,ax=None,case=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms] 
    #spikes[case] = spikes_df.to_json()
    spikes[case] = {'timestamps':spikes_df['timestamps'].tolist(), 'node_ids':spikes_df['node_ids'].tolist()}
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        ax.scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels))
    ax.grid(True)

def raw_ecp(lfp):
    pass

def ecp_psd(ecps,skip_n=0,downsample=20,nfft=1024,fs=1000,noverlap=0,ax=None,case=None):
    
    # Average all the runs
    #ecp = np.mean(np.array(ecps),axis=0)
    #skip_n first few
    pxxs = []
    for ecp in ecps:
        data = ecp[skip_n:]

        #downsample the data to fit ms (steps used 20=1/.05 step)
        lfp_d = decimate(data,downsample)
        raw_ecp(lfp_d)
        win = hanning(nfft, True)

        f,pxx = welch(lfp_d,fs,window=win,noverlap=noverlap,nfft=nfft)

        pxxs.append(pxx)    

    pxx = np.mean(np.array(pxxs),axis=0)
    # check ach theta adjustment match up
    #if case == 1:
    #    pxx = pxx - 0.000000054
    #if case == 2:
    #    pxx = pxx - 0.000000053    

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(f, pxx*1000,linewidth=0.6)
    ax.set_ylim([0,0.1])
    
    theta = pxx[np.where((f>=4) & (f<=8))]*1000
    gamma = pxx[np.where((f>=50) & (f<=60))]*1000
    mean_theta = theta.mean()
    peak_theta = theta.max() 
    mean_gamma = gamma.mean()
    peak_gamma = gamma.max()
    print('')
    print("Mean theta (4Hz-8Hz)  : " + str(round(mean_theta,6)))
    #print("Mean gamma (50Hz-60Hz) : " + str(round(mean_gamma,6)))     
    #print('')
    print("Peak theta (4Hz-8Hz)  : " + str(round(peak_theta,6)))
    #print("Peak gamma (50Hz-60Hz) : " + str(round(peak_gamma,6)))
    print('')

    psds[case] = {'f':f.tolist(),'pxx':(pxx*1000).tolist()}

def spike_frequency_histogram(spikes_df,node_set,ms,skip_ms=0,ax=None,n_bins=10,case=None):
    print("Type : mean (std)")
    spike_hist[case] = {}
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        #skip the first few ms
        cell_spikes = cell_spikes[cell_spikes['timestamps']>skip_ms]
        spike_counts = cell_spikes.node_ids.value_counts()
        total_seconds = (ms-skip_ms)/1000
        spike_counts_per_second = spike_counts / total_seconds

        spikes_mean = spike_counts_per_second.mean()
        spikes_std = spike_counts_per_second.std()
        node_name = node['name']
        spike_hist[case][node_name] = {'mean':spikes_mean, 'std':spikes_std}

        label = "{} : {:.2f} ({:.2f})".format(node['name'],spikes_mean,spikes_std)
        print(label)
        c = "tab:" + node['color']
        if ax:
            ax.hist(spike_counts_per_second,n_bins,density=True,histtype='bar',label=label,color=c)
    if ax:
        ax.set_xscale('log')
        ax.legend() 
        
        
def run(case,show_plots=False,save_plots=False):
    

    dt = 0.05
    steps_per_ms = 1/dt
    skip_seconds = 5
    skip_ms = skip_seconds*1000
    skip_n = int(skip_ms * steps_per_ms)
    end_ms = 15000

    spikes_location =  'case' + str(case) + '/run1' + '/spikes.h5'
    
    print("loading " + spikes_location)
    f = h5py.File(spikes_location)
    spikes_df = pd.DataFrame({'node_ids':f['spikes']['BLA']['node_ids'],'timestamps':f['spikes']['BLA']['timestamps']}) 
    print("done")

    if show_plots or save_plots:
        case_path = 'case' + str(case)
        run_paths = []
        for files in os.listdir(case_path):
            f = os.path.join(case_path,files)
            if os.path.isdir(f) and not files.startswith('_'):
               run_paths.append(f)
        
        ecps = []
        for run_path in run_paths:
            ecp_h5_location = os.path.join(run_path,'ecp.h5')

            print("loading " + ecp_h5_location)
            ecp_channel = 0
            f = h5py.File(ecp_h5_location)
            data_raw = np.array(f['ecp']['data'])
            ecp = data_raw.T[ecp_channel] #flip verts and grab channel 0
            ecps.append(ecp)

    if show_plots or save_plots:
        print("plotting...")
        fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,4.8))#6.4,4.8 default
        fig.suptitle('Amygdala Theta Analysis Case ' + str(case))
        ecp_psd(ecps, skip_n=skip_n, ax=ax2,case=case)
        spike_frequency_histogram(spikes_df,node_set,end_ms,skip_ms=skip_ms,ax=ax3,case=case)
        raster(spikes_df,node_set,skip_ms=skip_ms,ax=ax1,case=case)
        if save_plots:
            f_name = 'case' + str(case) + '_analysis.png'
            print("saving " + f_name)
            plt.savefig(f_name, bbox_inches='tight')
        if show_plots:
            print("showing plots...")
            fig.tight_layout()
            plt.show()
         
    else:
        spike_frequency_histogram(spikes_df,node_set,end_ms,skip_ms=skip_ms,case=case)

def save_data():
    full_dict = {}
    full_dict['spikes'] = spikes
    full_dict['spike_hist'] = spike_hist
    full_dict['psds'] = psds

    with open("analysis.json", "w") as fp:
        json.dump(full_dict, fp, indent=2)


def final_plots():
    # SETUP
    analysis_json = 'analysis.json'
    print("loading " + analysis_json)
    with open(analysis_json,'r') as fp:
        analysis = json.load(fp)
    print(analysis_json + " loaded")
    spikes = analysis['spikes']
    spike_hist = analysis['spike_hist']
    psds = analysis['psds']
    psd_power = {}

    fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,4.8))#6.4,4.8 default
    #fig.suptitle('Amygdala Analysis')    
    
    # AX1 - Base raster
    ax1.set_title("Base Raster")
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Cell ID")
    #spikes[case] = spikes_df.to_json()
    case = "1"
    spikes_df = pd.DataFrame({'timestamps':spikes[case]['timestamps'], 'node_ids':spikes[case]['node_ids']})    
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        ax1.scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])

    handles,labels = ax1.get_legend_handles_labels()
    ax1.legend(reversed(handles), reversed(labels))
    ax1.grid(True)

    # AX2 - Theta band
    ax2.set_title("Theta Band PSD by Case")
    ax2.set_xlabel("Hz")
    ax2.set_ylabel("PSD [V^2/Hz]")
    labels = {
        "1": "1. Base",
        "2": "2. VPSI",
        "3": "3. ACH High",
        "4": "4. ACH Low",
        "5": "5. ACH High/No VPSI",
        "6": "6. ACH Low/No VPSI"
    }
    for i in range(6):
        case = str(i + 1)
        psd_case = psds[case]
        f = np.array(psd_case['f'])
        fx = f[np.where((f>=4) & (f<=8))]
        pxx = np.array(psd_case['pxx'])
        theta = pxx[np.where((f>=4) & (f<=8))]
        ax2.plot(fx,theta,linewidth=0.6,label=labels[case])
        psd_power[case] = scipy.integrate.simps(theta)
        print(f"PSD Theta Power for case {case}: {psd_power[case]}")
    ax2.legend()

    # AX3 - Power
    ax3.set_title("Theta Band Power by Case")
    ax3.set_xlabel("Case")
    ax3.set_ylabel("Power")
    for i in range(6):
        case = str(i + 1)
        ax3.bar(i+1,psd_power[case], label=labels[case])
    #ax3.legend()
    plt.tight_layout() 
    plt.savefig('analysis_final.png', bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    show_plots = False
    save_plots = True
    if '--show-plots' in sys.argv:
        show_plots = True
    if '--save-plots' in sys.argv:
        save_plots = True
    if '--final' in sys.argv:
        final_plots()
    else:
        for i in range(6):
            run(i+1, show_plots = show_plots, save_plots = save_plots)
        save_data()
