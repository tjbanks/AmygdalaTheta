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

def run(ecp_h5_location, spikes_h5_location, tstart=5000, tend=300001, low_band=4, high_band=12, fs=1000, bin_size=0.1, top_percentage=0.2):

    theta_band = get_band(ecp_h5_location, low_band, high_band, fs)
    spikes_df = get_spikes(spikes_h5_location,skip_ms=tstart)
    hilbert_trans = hilbert(theta_band)

    start = int(np.round(tstart / 1000 * fs))
    end = int(np.round(tend / 1000 * fs))
    fig,ax = plt.subplots(len(node_set)+1,2,figsize=(10,9.6))

    x_ax = np.arange(start,end)
    #ax[0].plot(x_ax,theta_band[start:end])
    hilbert_power = np.abs(hilbert_trans[start:end])
    #ax[0].plot(x_ax,hilbert_power)
    ax[0,0].set_title("Theta")
    ax[0,1].set_title("Theta")

    hilbert_phase = np.angle(hilbert_trans[start:end])

    nbins = int(2*np.pi//bin_size)+1
    print(nbins)
    bins = np.linspace(-np.pi,np.pi,nbins+1)

    x = np.linspace(-np.pi,3*np.pi,1000)
    y = np.cos(x)
    ax[0,0].plot(x,y)
    ax[0,1].plot(x,y)

    def plot_phase(cell_spikes, ax):
        # scatter plot
        #ax[i+1].scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],c='tab:'+node['color'])
        # spike time histogram
        #n,b = np.histogram(cell_spikes['timestamps'],bins=np.arange(start,end,bin_size))
        #ax[i+1].stairs(n,b,color=node['color'])
        spike_times = np.round((cell_spikes['timestamps']-tstart) / 1000 * fs).astype(int)

        phases = hilbert_phase[spike_times]
        n, b = np.histogram(phases, bins=bins)
        # normalize by total number
        n = n / n.sum()
        # going to repeat to extend
        b = np.concatenate([b, b[1:]+2*np.pi])
        n = np.tile(n, 2)
        
        ax2 = ax.twinx()
        ax2.plot(x,y,linestyle='dashed')
        
        ax.stairs(n,b,color=node['color'],fill=True)
            
    for i, node in enumerate(node_set):
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]
        top_spiking_cells = cell_spikes.node_ids.value_counts()[:int(top_percentage*len(cells))].index.tolist()
        top_cell_spikes = spikes_df[spikes_df['node_ids'].isin(top_spiking_cells)]
        #import pdb;pdb.set_trace()
        ax[i+1,0].set_title(node['name'])
        ax[i+1,1].set_title('top ' + str(int(top_percentage*100)) + '% ' + node['name'])

        plot_phase(cell_spikes, ax[i+1,0])    
        plot_phase(top_cell_spikes, ax[i+1,1])

    plt.tight_layout()
    plt.show()
if __name__ == '__main__':
    run('./outputECP/ecp.h5', './outputECP/spikes.h5')
