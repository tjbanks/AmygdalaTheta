import h5py
import numpy as np
import matplotlib.pyplot as plt

scale=1

def convert(num_pn,inp='./spikesmatrix_op_ryt', out='../vpsi_inh_spikes.h5'):
    matlab_file = open(inp, 'r')
    lines = matlab_file.readlines()
    matlab_file.close()

    vpsi = h5py.File(out, 'w')

    vpsi_spikes = vpsi.create_group("spikes")
    vpsi_spikes_vp = vpsi_spikes.create_group("vpsi_inh")

    spikes_node_ids = []
    spikes_timestamps = []
 
    count = 0
    for line in lines:
        times = line.strip().split('\t')
        if len(times)==1 and times[0]=='':
            continue # There were no spikes for this cell
        
        ids = [count for _ in range(len(times))]

        spikes_timestamps = spikes_timestamps + times
        spikes_node_ids = spikes_node_ids + ids
             
        count = count + 1
 
    nodes=np.array(spikes_node_ids).astype(np.int)
    timestamps=np.array(spikes_timestamps).astype(np.float)

    vpsi_spikes_vp.create_dataset("node_ids", data=nodes)
    vpsi_spikes_vp.create_dataset("timestamps", data=timestamps)

    vpsi.close()

    # plot a scatter plot
    fig, ax = plt.subplots()
    ax.scatter(timestamps,nodes,s=2)
    ax.set_xlabel('Spike Timestamp', fontsize=12)
    ax.set_ylabel('Affrent Neuron ID', fontsize=12)
    ax.set_title('8Hz Input', fontsize=12)

    #ax.grid(True)
    fig.tight_layout()
    plt.show()
    
if __name__ == '__main__':
    convert(num_pn=800*scale)
