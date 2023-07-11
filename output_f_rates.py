import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# generate spike times (get these from your model)
#spiketimes = np.random.choice(np.arange(0,10000,1),1000)
#n_cells = 100
scale = 1

spikes_location = 'outputECP/spikes.h5'
f = h5py.File(spikes_location)
spikes_df = pd.DataFrame({'node_ids':f['spikes']['BLA']['node_ids'],'timestamps':f['spikes']['BLA']['timestamps']})
spiketimes = spikes_df['timestamps']
n_cells = 1000 * scale #len(spikes_df.node_ids.unique()) # we want to include those that didn't fire

def raster(spikes_df,skip_ms=0,ax=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms]

    plt.figure()
    plt.scatter(spikes_df['timestamps'],spikes_df['node_ids'],
                s=0.25 )

    plt.grid(True)

# do histogram (choose bin size)
bin_size=5
n,b = np.histogram(spiketimes,bins=np.arange(0,np.max(spiketimes),bin_size))

def moving_average(a, n=5):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

spikes_second_cell = n/((bin_size/1000)*n_cells) # spikes per second per cell

plt.figure()
plt.plot(spikes_second_cell, label = 'raw bins')
plt.plot(moving_average(spikes_second_cell), label = 'moving average')
plt.xticks(ticks = np.arange(0,2250,250), labels = ['{}'.format(5*i) for i in np.arange(0,2250,250)])
plt.xlabel('time (ms)')
plt.ylabel('firing rate (Hz)')
plt.legend()

plt.figure()
plt.psd(moving_average(spikes_second_cell), Fs = 200)

raster(spikes_df)
plt.show()
