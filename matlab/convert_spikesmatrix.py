import h5py
import numpy as np

def convert(num_pn,inp='./spikesmatrix_op_ryt', out_pn='../vpsi_pyr_inh_spikes.h5', out_pv='../vpsi_pv_inh_spikes.h5'):
    file1 = open(inp, 'r')
    lines = file1.readlines()
    file1.close()

    pn = h5py.File(out_pn, 'w')
    pv = h5py.File(out_pv, 'w')    
    pn_key = 'vpsi_pn_inh'
    pv_key = 'vpsi_pv_inh'

    pn_spikes = pn.create_group("spikes")
    pv_spikes = pv.create_group("spikes")
    pn_spikes_vp = pn_spikes.create_group(pn_key)
    pv_spikes_vp = pv_spikes.create_group(pv_key)

    pn_spikes_node_ids = []
    pn_spikes_timestamps = []
    pv_spikes_node_ids = []
    pv_spikes_timestamps = []
 
    count = 0
    for line in lines:
        times = line.strip().split('\t')
        if len(times)==1 and times[0]=='':
            continue
        ids = [count for _ in range(len(times))]

        if count < num_pn:
            pn_spikes_timestamps = pn_spikes_timestamps + times
            pn_spikes_node_ids = pn_spikes_node_ids + ids
        else:
            ids = [count-num_pn for _ in range(len(times))]
            pv_spikes_timestamps = pv_spikes_timestamps + times
            pv_spikes_node_ids = pv_spikes_node_ids + ids
             
       
        count = count + 1
 
    pn_spikes_vp.create_dataset("node_ids", data=np.array(pn_spikes_node_ids).astype(np.int))
    pn_spikes_vp.create_dataset("timestamps", data=np.array(pn_spikes_timestamps).astype(np.float))

    pv_spikes_vp.create_dataset("node_ids", data=np.array(pv_spikes_node_ids).astype(np.int))
    pv_spikes_vp.create_dataset("timestamps", data=np.array(pv_spikes_timestamps).astype(np.float))

    pn.close()
    pv.close()

if __name__ == '__main__':
    convert(num_pn=800)
