from bmtk.analyzer.compartment import plot_traces

_ = plot_traces(config_file='simulation_config.json', node_ids=[0,1,2], report_name='v_report')



import h5py
e_s = h5py.File('mthalamus_spikes.h5')
i_s = h5py.File('exc_bg_bask_spikes.h5')
e_nodes = list(e_s['spikes']['mthalamus']['node_ids'])
e_times = list(e_s['spikes']['mthalamus']['timestamps'])
i_nodes = list(i_s['spikes']['exc_bg_bask']['node_ids'])
i_times = list(i_s['spikes']['exc_bg_bask']['timestamps'])

i_nodes = [i+2 for i in i_nodes]

import matplotlib.pyplot as plt
girls_grades = [89, 90, 70, 89, 100, 80, 90, 100, 80, 34]
boys_grades = [30, 29, 49, 48, 100, 48, 38, 45, 20, 30]
grades_range = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(e_times, e_nodes, color='r')
ax.scatter(i_times, i_nodes, color='r')
ax.set_xlabel('node_id')
ax.set_ylabel('timestamps')
ax.set_title('input_spike times')
plt.show()

