import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter,filtfilt,hilbert

scale = 1
node_set = [
        #{"name":"PN","start":0*scale,"end":799*scale,"color":"blue"},
        {"name":"PNa","start":0*scale,"end":639*scale,"color":"blue"},
        {"name":"PNc","start":640*scale,"end":799*scale,"color":"blue"},
        {"name":"PV","start":800*scale,"end":892*scale,"color":"red"},
        {"name":"SOM","start":893*scale,"end":943*scale,"color":"green"},
        {"name":"CR","start":944*scale,"end":999*scale,"color":"purple"}
    ]

def butter_bandpass(lowcut, highcut, fs, order=5):
    return butter(order, [lowcut, highcut], fs=fs, btype='band')

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def get_band(ecp_h5_location, low, high, fs):
    print("loading " + ecp_h5_location)
    ecp_channel = 0
    f = h5py.File(ecp_h5_location)
    data_raw = np.array(f['ecp']['data'])
    ecp = data_raw.T[ecp_channel] #flip verts and grab channel 0

    # Get theta band for plotting
    theta_band = butter_bandpass_filter(ecp,low,high,fs).tolist()
    return theta_band

def get_spikes(spikes_h5_location,skip_ms=0):
    print("loading " + spikes_h5_location)
    f = h5py.File(spikes_h5_location)
    spikes_df = pd.DataFrame({'node_ids':f['spikes']['BLA']['node_ids'],'timestamps':f['spikes']['BLA']['timestamps']})
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms]
    return spikes_df

def plot_raster(spikes_df,node,skip_ms=0,ax=None,start=None, end=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms]
    if start and end:
        spikes_df = spikes_df[(spikes_df['timestamps']>start) & (spikes_df['timestamps']<end)]

    ax.scatter(spikes_df['timestamps'],spikes_df['node_ids'],
            c='tab:'+node['color'],s=0.25, label=node['name'])

    handles,labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels))
    ax.grid(True)
    return spikes_df

def plot_phase(ecp_h5_location, spikes_h5_location, tstart=5000, tend=15000, low_band=4, high_band=12, fs=1000, dt=0.1, bin_size=0.1, top_percentage=0.2, show=False, title=None,plot_vpsi=True,vpsi_location='./vpsi_inh_spikes.h5'):

    theta_band = get_band(ecp_h5_location, low_band, high_band, fs)
    spikes_df = get_spikes(spikes_h5_location,skip_ms=tstart)
    hilbert_trans = hilbert(theta_band)

    start = int(np.round(tstart / dt))
    end = int(np.round(tend / dt))
    extra_plots = 2 if plot_vpsi else 1
    fig,ax = plt.subplots(len(node_set)+extra_plots,2,figsize=(10,9.6))

    x_ax = np.arange(start,end)
    #ax[0].plot(x_ax,theta_band[start:end])
    hilbert_power = np.abs(hilbert_trans[start:end])
    print(f"Band mean power {hilbert_power.mean()} | std power {hilbert_power.std()}")
    print(f"Band max power {hilbert_power.max()} | min power {hilbert_power.min()} ")
    #ax[0].plot(x_ax,hilbert_power)
    ax[0][0].set_title("Theta")
    #ax[0,1].set_title("Theta")

    hilbert_phase = np.angle(hilbert_trans[start:end])

    nbins = int(2*np.pi//bin_size)+1
    #print(nbins)
    bins = np.linspace(-np.pi,np.pi,nbins+1)

    x = np.linspace(-np.pi,3*np.pi,1000)
    y = np.cos(x)
    ax[0][0].plot(x,y)
    pd.DataFrame(np.array([x,y]).T,columns=['x','y']).to_csv(os.path.join(output_folder,'figure6b_cosin_wave.csv'),index=False)
    #ax[0,1].plot(x,y)

    def plot_phase_inner(cell_spikes, ax, power_threshold=0):
        # scatter plot
        #ax[i+1].scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],c='tab:'+node['color'])
        # spike time histogram
        #n,b = np.histogram(cell_spikes['timestamps'],bins=np.arange(start,end,bin_size))
        #ax[i+1].stairs(n,b,color=node['color'])
        spike_times = np.round((cell_spikes['timestamps']-tstart)/dt).astype(int)
        
        phases = hilbert_phase[spike_times]
        
        if power_threshold: # if set, we only want to count those that occur in high peaks of the filtered signal
            powers = hilbert_power[spike_times]
            powerful_indicies = [i for i, elem in enumerate(powers) if elem > power_threshold]
            phases = phases[powerful_indicies]
        
        n, b = np.histogram(phases, bins=bins)
        # normalize by total number
        n = n / n.sum()
        # going to repeat to extend
        b = np.concatenate([b, b[1:]+2*np.pi])
        n = np.tile(n, 2)
        
        ax2 = ax.twinx()
        ax2.plot(x,y,linestyle='dashed')
        
        ax.stairs(n,b,color=node['color'],fill=True)
        return n, b
    
    n_list = []
    b_list =[]
    spikes_df_list = []
    vpsi_spikes_df = None 

    for i, node in enumerate(node_set):
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]
        top_spiking_cells = cell_spikes.node_ids.value_counts()[:int(top_percentage*len(cells))].index.tolist()
        top_cell_spikes = spikes_df[spikes_df['node_ids'].isin(top_spiking_cells)]

        ax[i+1][0].set_title(node['name'])
        #ax[i+1,1].set_title('top ' + str(int(top_percentage*100)) + '% ' + node['name'])

        n, b = plot_phase_inner(cell_spikes, ax[i+1][0])
        n_list.append(n)
        b_list.append(b)
        #plot_phase_inner(top_cell_spikes, ax[i+1,1])
        spikes_df_ret = plot_raster(cell_spikes, node, ax=ax[i+1][1],start=8000,end=8500)
        spikes_df_list.append(spikes_df_ret)

    if plot_vpsi:
        ax[6][0].set_title("VPSI input")
        vpsi_input = h5py.File(vpsi_location)
        vpsi_spikes = pd.DataFrame({'node_ids':vpsi_input['spikes']['vpsi_inh']['node_ids'],
                                'timestamps':vpsi_input['spikes']['vpsi_inh']['timestamps']})
        vpsi_spikes = vpsi_spikes[vpsi_spikes['timestamps']>5000]
        n_vpsi, b_vpsi = plot_phase_inner(vpsi_spikes, ax[6][0])
        n_list.append(n_vpsi)
        b_list.append(b_vpsi)
        node['name']='VPSI'
        vpsi_spikes_df = plot_raster(vpsi_spikes, node , ax=ax[6][1],start=8000,end=8500)
        spikes_df_list.append(vpsi_spikes_df)

    #plt.tight_layout()
    if title:
        fig.suptitle(title)
    fig.set_layout_engine('tight')
    if show:
        plt.show()

    return n_list, b_list, spikes_df_list 

if __name__ == '__main__':
    output_folder = "figures_data"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if '--all' in sys.argv:
        for case in range(1,7):
            plot_vpsi = True if case in [2,3,4] else False
            (n_list, b_list, spikes_df_list) = plot_phase(f'./case{case}/run1/ecp.h5', f'./case{case}/run1/spikes.h5', show=False, title=f"Case {case}", plot_vpsi=plot_vpsi, vpsi_location='../vpsi_inh_spikes.h5')
            if case in [2]:
                nodes = [n["name"] for n in node_set]
                nodes.append("vpsi")
                # save the bounds
                pd.DataFrame(np.array(b_list[0]).T,columns=['bounds']).to_csv(os.path.join(output_folder, f'figure6b_case_{case}_phase_bounds.csv'), index=False)
                # save the values
                pd.DataFrame(np.array(n_list).T,columns=nodes).to_csv(os.path.join(output_folder, f'figure6b_case_{case}_phase_values.csv'), index=False)
                # save the spikes_list
                for node_name, spikes_df in zip(nodes, spikes_df_list):
                    spikes_df.to_csv(os.path.join(output_folder, f'figure6b_case_{case}_{node_name}_raster.csv'), index=False)
        plt.show()
    else:
        plot_phase('./outputECP/ecp.h5', './outputECP/spikes.h5', show=True)
