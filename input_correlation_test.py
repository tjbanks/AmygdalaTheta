import random

import h5py
import pandas as pd
from scipy.stats import pearsonr, spearmanr

input_file = 'thalamus_pyr_spikes.h5'
spikes_path = 'thalamus_pyr'

#number of sample timestamps to pull
sample_size = 300

f = h5py.File(input_file)

spikes_df = pd.DataFrame({'node_ids':f['spikes'][spikes_path]['node_ids'],'timestamps':f['spikes'][spikes_path]['timestamps']}) 
num_cells = max(spikes_df['node_ids'])

cell_1 = int(random.random()*num_cells)
cell_2 = int(random.random()*num_cells)

cell_1_times = spikes_df[spikes_df['node_ids'] == cell_1]['timestamps'].tolist()
cell_2_times = spikes_df[spikes_df['node_ids'] == cell_2]['timestamps'].tolist()

cell_1_times_sample = random.sample(cell_1_times, sample_size)
cell_2_times_sample = random.sample(cell_2_times, sample_size)


corr, _ = pearsonr(cell_1_times_sample, cell_2_times_sample)
print('pearsonr correlation: %.3f' % corr)

corr, _ = spearmanr(cell_1_times_sample, cell_2_times_sample)
print('spearmanr correlation: %.3f' % corr)
