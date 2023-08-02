import matplotlib.pyplot as plt

# Import time-frequency functions
from neurodsp.timefrequency import freq_by_time

fs = 1000

# Set the frequency range to be used
f_range = (13, 30)

scale = 1

spikes_location = 'vpsi_inh_spikes.h5'
f = h5py.File(spikes_location)
spikes_df = pd.DataFrame({'node_ids':f['spikes']['vpsi_inh']['node_ids'],'timestamps':f['spikes']['vpsi_inh']['timestamps']})
spiketimes = spikes_df['timestamps']
n_cells = 893 * scale

# Compute instantaneous frequency from a signal
i_f = freq_by_time(sig, fs, f_range)
