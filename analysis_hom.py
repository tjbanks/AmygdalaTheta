"""
A python implementation of matlab/analysis.m

TB - 8/4/21
"""
from fooof import FOOOF
from fooof.sim.gen import gen_aperiodic
from scipy.signal.windows import hann as hanning
from scipy.signal import welch,decimate
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

scale = 1

def raster(spikes_df,node_set,skip_ms=0,ax=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms] 
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        ax.scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels))
    ax.grid(True)
    ax.set_xlim(8500, 9000)

def raw_ecp(lfp):
    pass

def ecp_psd(ecp,skip_n=0,downsample=10,nfft=1024,fs=1000,noverlap=0,ax=None,use_fooof=True):
    
    #skip_n first few
    data = ecp[skip_n:]

    #downsample the data to fit ms (steps used 20=1/.05 step)
    lfp_d = decimate(data,downsample)
    #raw_ecp(lfp_d)
    if not use_fooof:
        win = hanning(nfft, True)
        f,pxx = welch(lfp_d,fs,window=win,noverlap=noverlap,nfft=nfft)
        #ax.set_xscale('log')
        ax.set_yscale('log')
        ax.plot(f, pxx*1000,linewidth=0.6)
        ax.set_ylim([0,0.1])

        theta = pxx[np.where((f>=4) & (f<=12))]*1000
        gamma = pxx[np.where((f>=50) & (f<=60))]*1000
        mean_theta = theta.mean()
        peak_theta = theta.max()
        mean_gamma = gamma.mean()
        peak_gamma = gamma.max()
    else:
        f,pxx = welch(lfp_d,fs=1000,nfft=1024)
        freqs,spectrum = np.array(f),np.array(pxx)
        fm = FOOOF(aperiodic_mode='knee')
        fm.fit(freqs, spectrum, [1,500])
        ap_fit = fm._ap_fit
        residual_spec = spectrum[0:500] #- 10**ap_fit

        # Plot
        #plt.plot([i for i in range(len(residual_spec))],residual_spec)
        ax.plot(freqs[:len(residual_spec)], residual_spec)
        ax.grid()
        #ax.set_xlim(1, 25)
        #ax.plot([i for i in range(4,13)],residual_spec[4:13]) # Only theta range
        
        theta = residual_spec[4:13]
        peak_theta = max(theta)
        mean_theta = sum(theta)/len(theta)
        gamma = residual_spec[50:61]
        peak_gamma = max(gamma)
        mean_gamma = sum(gamma)/len(gamma)
    
    text = f"""
    
Mean theta (4Hz-8Hz)  : {str(round(mean_theta,8))}
Mean gamma (50Hz-60Hz) : {str(round(mean_gamma,8))} 
    
Peak theta (4Hz-8Hz)  : {str(round(peak_theta,8))}
Peak gamma (50Hz-60Hz) : {str(round(peak_gamma,8))}
    
    """
    return text

def spike_frequency_histogram(spikes_df,node_set,ms,skip_ms=0,ax=None,n_bins=10):
    return_text = "Type : mean (std)\n"
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
        
        label = "{} : {:.2f} ({:.2f})".format(node['name'],spikes_mean,spikes_std)
        #print(label)
        return_text = return_text + label + '\n'
        c = "tab:" + node['color']
        if ax:
            ax.hist(spike_counts_per_second,n_bins,density=True,histtype='bar',label=label,color=c)
    if ax:
        ax.set_xscale('log')
        ax.legend() 
        
    return return_text 

def run(show_plots=False,save_plots=False,slack=True):
    

    dt = 0.05
    steps_per_ms = 1/dt
    skip_seconds = 5
    skip_ms = skip_seconds*1000
    skip_n = int(skip_ms * steps_per_ms)
    end_ms = 15000

    spikes_location = 'outputECP/spikes.h5'
    
    print("loading " + spikes_location)
    f = h5py.File(spikes_location)
    spikes_df = pd.DataFrame({'node_ids':f['spikes']['BLA']['node_ids'],'timestamps':f['spikes']['BLA']['timestamps']}) 
    print("done")

    if show_plots or save_plots:
        ecp_h5_location = 'outputECP/ecp.h5'
        print("loading " + ecp_h5_location)
        ecp_channel = 0
        f = h5py.File(ecp_h5_location)
        data_raw = np.array(f['ecp']['data'])
        ecp = data_raw.T[ecp_channel] #flip verts and grab channel 0
        print("done")

    node_set = [
        {"name":"PN","start":0*scale,"end":799*scale,"color":"blue"},
        {"name":"PV","start":800*scale,"end":892*scale,"color":"red"},
        {"name":"SOM","start":893*scale,"end":943*scale,"color":"green"},
        {"name":"CR","start":944*scale,"end":999*scale,"color":"purple"}
    ]
    
    if show_plots or save_plots:
        print("plotting...")
        fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,4.8))#6.4,4.8 default
        fig.suptitle('Amygdala Theta Analysis')
        output_text = ecp_psd(ecp, skip_n=skip_n, ax=ax2)
        sfh_text = spike_frequency_histogram(spikes_df,node_set,end_ms,skip_ms=skip_ms,ax=ax3)
        output_text = output_text + sfh_text
        print(output_text)
        raster(spikes_df,node_set,skip_ms=skip_ms,ax=ax1)
        if save_plots:
            f_name = 'analysis.png'
            print("saving " + f_name)
            plt.savefig(f_name, bbox_inches='tight')
        if show_plots:
            print("showing plots...")
            fig.tight_layout()
            plt.show()
        if slack:
            import upload_analysis_slack
            import subprocess
            output_text = output_text + '\n\n' + str(subprocess.check_output(['git', 'diff', 'components_homogenous/synaptic_models/']))
            upload_analysis_slack.upload(output_text)
    else:
        spike_frequency_histogram(spikes_df,node_set,end_ms,skip_ms=skip_ms)


if __name__ == '__main__':
    show_plots = False
    save_plots = False
    if '--show-plots' in sys.argv:
        show_plots = True
    if '--save-plots' in sys.argv:
        save_plots = True
        
    run(show_plots = show_plots, save_plots = save_plots)
