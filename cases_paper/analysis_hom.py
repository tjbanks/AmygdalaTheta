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

def raster(spikes_df,node_set,skip_ms=0,ax=None,case=None,title=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms] 
    #spikes[case] = spikes_df.to_json()
    if not case in spikes.keys():
        spikes[case] = []

    spikes[case].append({'timestamps':spikes_df['timestamps'].tolist(), 'node_ids':spikes_df['node_ids'].tolist()})
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        ax.scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels))
    ax.grid(True)
    ax.set_xlim(8500, 9000)
    if title:
        ax.set_title(title)

def raw_ecp(lfp):
    pass

def ecp_psd(ecps,skip_n=0,downsample=10,nfft=1024,fs=1000,noverlap=0,ax=None,case=None):
    
    # Average all the runs
    #ecp = np.mean(np.array(ecps),axis=0)
    #skip_n first few
    pxxs = []
    pxxs_raw = []
    for ecp in ecps:
        data = ecp[skip_n:]

        #downsample the data to fit ms (steps used 20=1/.05 step)
        lfp_d = decimate(data,downsample)
        #lfp_d = data
        raw_ecp(lfp_d)
        win = hanning(nfft, True)

        f,pxx = welch(lfp_d,fs,window=win,noverlap=noverlap,nfft=nfft)
        f_raw, pxx_raw = welch(ecp,fs=1000,nfft=1024)
        pxxs.append(pxx)    
        pxxs_raw.append(pxx_raw)

    pxx = np.mean(np.array(pxxs),axis=0)
    pxx_raw = np.mean(np.array(pxxs_raw),axis=0)
    # check ach theta adjustment match up
    #if case == 1:
    #    pxx = pxx - 0.000000054
    #if case == 2:
    #    pxx = pxx - 0.000000053    

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(f, pxx*1000,linewidth=0.6)
    ax.set_ylim([0,0.1])
    
    theta = pxx[np.where((f>=4) & (f<=12))]*1000
    gamma = pxx[np.where((f>=50) & (f<=60))]*1000
    mean_theta = theta.mean()
    peak_theta = theta.max() 
    mean_gamma = gamma.mean()
    peak_gamma = gamma.max()
    print('')
    print("Mean theta (4Hz-12Hz)  : " + str(round(mean_theta,6)))
    #print("Mean gamma (50Hz-60Hz) : " + str(round(mean_gamma,6)))     
    #print('')
    print("Peak theta (4Hz-12Hz)  : " + str(round(peak_theta,6)))
    #print("Peak gamma (50Hz-60Hz) : " + str(round(peak_gamma,6)))
    print('')

    psds[case] = {'f':f.tolist(),'pxx':(pxx*1000).tolist(),
                  'f_raw':f_raw.tolist(), 'pxx_raw':(pxx_raw*1000).tolist()}

def spike_frequency_histogram(spikes_df,node_set,ms,skip_ms=0,ax=None,n_bins=10,case=None,title=None):
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
    if title:
        ax.set_title(title)
        
        
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


def final_plots(num_cases=6):
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

    fig, ax = plt.subplots(3,3,figsize=(15,9.6))#6.4,4.8 default
    #fig.suptitle('Amygdala Analysis')    
    fig2,ax2 = plt.subplots(2,3,figsize=(15,9.6))
    fig3,ax3 = plt.subplots(2,3,figsize=(15,9.6))
    fig4,ax4 = plt.subplots(2,3,figsize=(15,9.6))

    # AX1 - Base raster
    ax[0,0].set_title("Base Raster")
    ax[0,0].set_xlabel("Time")
    ax[0,0].set_ylabel("Cell ID")
    #spikes[case] = spikes_df.to_json()
    case = "1"
    spikes_df = pd.DataFrame({'timestamps':spikes[case][0]['timestamps'], 'node_ids':spikes[case][0]['node_ids']})    
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        ax[0,0].scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])

    handles,labels = ax[0,0].get_legend_handles_labels()
    ax[0,0].legend(reversed(handles), reversed(labels))
    ax[0,0].grid(True)

    # AX2 - Theta band
    ax[0,1].set_title("Theta Band PSD by Case")
    ax[0,1].set_xlabel("Hz")
    ax[0,1].set_ylabel("PSD [V^2/Hz]")
    labels = {
        "1": "1.GABA + no Theta + ACh baseline",
        "2": "2.GABA + Theta + ACh baseline",
        "3": "3.GABA + Theta + ACH High",
        "4": "4.GABA + Theta + ACH Low",
        "5": "5.GABA + no Theta + ACH High",
        "6": "6.GABA + no Theta + ACH Low"
    }
    from fooof import FOOOF
    from fooof.sim.gen import gen_aperiodic
    use_fooof = True
    use_peak = True

    def plot_psd(axis, case_int, max_freq=150, legend=True):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        
    #for i in range(num_cases):
        case = str(case_int)
        psd_case = psds[case]
        if not use_fooof:
            f = np.array(psd_case['f'])
            fx = f[np.where((f>=4) & (f<=12))]
            pxx = np.array(psd_case['pxx'])
            theta = pxx[np.where((f>=4) & (f<=12))]
            axis.plot(fx,theta,linewidth=0.6,label=labels[case],color=colors[case_int-1])
        else:
            freqs,spectrum = np.array(psd_case['f']),np.array(psd_case['pxx'])
            fm = FOOOF(aperiodic_mode='knee')
            fm.fit(freqs, spectrum, [1,max_freq])
            ap_fit = fm._ap_fit
            print(ap_fit)
            residual_spec = spectrum[0:max_freq+2] - 10**ap_fit
            #ax2.plot([i for i in range(4,13)],residual_spec[4:13])
            axis.plot([i for i in range(len(residual_spec))],residual_spec,color=colors[case_int-1], label=labels[case])
            theta = residual_spec[4:13]
            # original
            #freqs,spectrum = np.array(psd_case['f_raw']),np.array(psd_case['pxx_raw'])
            #fm = FOOOF(aperiodic_mode='knee')
            #fm.add_data(freqs,spectrum,[1,150])
            #fm.fit(freqs, spectrum, [1,150])
            #ap_fit = gen_aperiodic(fm.freqs, fm._robust_ap_fit(fm.freqs, fm.power_spectrum))
            #init_flat_spec = fm.power_spectrum - 10**ap_fit
            #theta = (init_flat_spec[4:13])
            #ax2.plot([i for i in range(4,13)],theta)
        axis.set_xscale('log')
        if legend:
            axis.legend()
        return theta

    for case in range(1,num_cases+1):
        theta = plot_psd(ax[0,1],case,legend=False)
        if use_peak:
            psd_power[str(case)] = max(theta)
        else: # integrage
            psd_power[str(case)] = scipy.integrate.simps(theta)
        print(f"PSD Theta Power for case {case}: {psd_power[str(case)]}")

    # AX3 - Power
    if use_peak:
        ax[0,2].set_title("Theta Band Peak by Case")
    else:
        ax[0,2].set_title("Theta Band Power by Case")
    ax[0,2].set_xlabel("Case")
    if use_peak:
        ax[0,2].set_ylabel("[V^2/Hz]")
    else:
        ax[0,2].set_ylabel("Power")
    for i in range(num_cases):
        case = str(i + 1)
        ax[0,2].bar(i+1,psd_power[case], label=labels[case])
    #ax3.legend()

    comparisons = [
            {'axis':ax[1,0], "cases":[1,2]},
            {'axis':ax[1,1], "cases":[1,3]},
            {'axis':ax[1,2], "cases":[1,4]},
            {'axis':ax[2,0], "cases":[1,5]},
            {'axis':ax[2,1], "cases":[1,6]},
            {'axis':ax[2,2], "cases":[3,5]},
        ]
    for compare in comparisons:
        title = "Case "
        for not_first, case in enumerate(compare["cases"]):
            plot_psd(compare["axis"], case)
            if not_first:
                title = title + " vs Case "
            title = title + str(case)
        compare["axis"].set_title(title)


    plt.savefig('analysis_final.svg', bbox_inches='tight', format='svg')
    #plt.show()

    dt = 0.05
    steps_per_ms = 1/dt
    skip_seconds = 5
    skip_ms = skip_seconds*1000
    skip_n = int(skip_ms * steps_per_ms)
    end_ms = 15000

    
    ### Firing rates
    def plot_f_rates(spikes_df, node_set, ax=None, bin_size=100, title=None,skip_ms=0):

        for node in node_set:
            cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
            cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]
            spiketimes = spikes_df['timestamps']
            n_cells = len(cells) #len(spikes_df.node_ids.unique()) # we want to include those that didn't fire
            
            # do histogram (choose bin size)
        
            n,b = np.histogram(spiketimes,bins=np.arange(skip_ms,np.max(spiketimes),bin_size))

            def moving_average(a, n=5):
                ret = np.cumsum(a, dtype=float)
                ret[n:] = ret[n:] - ret[:-n]
                return ret[n - 1:] / n

            spikes_second_cell = n/((bin_size/1000)*n_cells) # spikes per second per cell
            ax.plot(spikes_second_cell, label = node['name']) #'raw bins')
            #ax.plot(moving_average(spikes_second_cell), label = node['name']) #'moving average')

            #ax.plot(something, label = node['name'])

        #ax.set_xticks(ticks = np.arange(5000,15000,bin_size), labels = ['{}'.format(5*i) for i in np.arange(5000,15000,bin_size)])
        ax.set_xlabel('time (ms)')
        ax.set_ylabel('firing rate (Hz)')
        if title:
            ax.set_title(title)
        ax.legend()

    for i, case in enumerate(range(1,num_cases+1)):
        case = str(case)
        spikes_df = pd.DataFrame({'timestamps':spikes[case][0]['timestamps'], 'node_ids':spikes[case][0]['node_ids']})
        col = i % 3
        row = 0 if i < 3 else 1
        # figure 2
        spike_frequency_histogram(spikes_df,node_set,end_ms,skip_ms=skip_ms,case=case,ax=ax2[row,col],title=labels[case])
        # figure 3
        raster(spikes_df,node_set,skip_ms=skip_ms,ax=ax3[row,col],case=case,title=labels[case])
        # figure 4
        plot_f_rates(spikes_df,node_set,ax=ax4[row,col],title=labels[case],skip_ms=skip_ms)


    #plt.tight_layout()
    fig.set_layout_engine('tight')
    fig2.set_layout_engine('tight')
    fig3.set_layout_engine('tight')
    fig4.set_layout_engine('tight')
    plt.show()

if __name__ == '__main__':
    show_plots = False
    save_plots = True
    num_cases = 6
    if '--show-plots' in sys.argv:
        show_plots = True
    if '--save-plots' in sys.argv:
        save_plots = True
    if '--final' in sys.argv:
        final_plots(num_cases)
    else:
        for i in range(num_cases):
            run(i+1, show_plots = show_plots, save_plots = save_plots)
        save_data()
