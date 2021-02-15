from bmtk.analyzer.compartment import plot_traces

_ = plot_traces(config_file='simulation_config.json', node_ids=[2], report_name='v_report')

exit()

import h5py
e_s = h5py.File('mthalamus_spikes.h5')
e_nodes = list(e_s['spikes']['mthalamus']['node_ids'])
e_times = list(e_s['spikes']['mthalamus']['timestamps'])

i_nodes = [i+2 for i in i_nodes]

import matplotlib.pyplot as plt
fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(e_times, e_nodes, color='r')
ax.set_xlabel('node_id')
ax.set_ylabel('timestamps')
ax.set_title('input_spike times')
plt.show()

