import sys

import numpy as np
import pandas as pd
import h5py

from bmtools.cli.plugins.util.util import relation_matrix

f = h5py.File('outputECP/spikes.h5')
spikes = pd.DataFrame({'node_ids':f['spikes']['BLA']['node_ids'],'timestamps':f['spikes']['BLA']['timestamps']})

edges = None

def conn_info(**kwargs):
    global edges
    _edges = kwargs["edges"]
    source_id_type = kwargs["sid"]
    target_id_type = kwargs["tid"]
    source_id = kwargs["source_id"]
    target_id = kwargs["target_id"]
    t_list = kwargs["target_nodes"]
    s_list = kwargs["source_nodes"]
    
    if edges is None: #want to accumulate connections
        edges = _edges
    else:
        edges = edges.append(_edges).drop_duplicates()

    #center_x = kwargs["center_x"]
    #center_y = kwargs["center_y"]
    #center_z = kwargs["center_z"]    
    
    #step = kwargs["step"]
    #iterations = kwargs["iterations"]
    
    center_x = 600
    center_y = 600
    center_z = 600

    iterations = 13
    step = 25
    ms = 15000
    skip_ms = 5000

    cons = edges[(edges[source_id_type] == source_id) & (edges[target_id_type]==target_id)]
    # Disable this next assignment ifor most normal situatipns
    cons = edges[(edges[target_id_type]==target_id)]
    #total_cons = cons.count().source_node_id
    #import pdb;pdb.set_trace() 
    print(source_id + " -> " + target_id + " step " + str(step*2))
    print("from locations => to locations : n number of cells : connections mean (std) : spikes mean (std)")
    
    biggest_jump = 0
    biggest_jump_str = ""
    last = 0
    last_str = ""
    
    dist = np.linalg.norm(edges[['target_pos_x', 'target_pos_y', 'target_pos_z']].values.astype(float)-edges[['source_pos_x', 'source_pos_y', 'source_pos_z']].values.astype(float),axis=1) 
    
    for inner_dist in range(step,iterations*step,step):
        
        #cube
        min_x = center_x - inner_dist
        max_x = center_x + inner_dist
        min_y = center_y - inner_dist
        max_y = center_y + inner_dist
        min_z = center_z - inner_dist
        max_z = center_z + inner_dist

        connections = cons[(cons['target_pos_x'] > min_x) & (cons['target_pos_y'] > min_x ) & (cons['target_pos_z'] > min_z)
                         & (cons['target_pos_x'] < max_x) & (cons['target_pos_y'] < max_y ) & (cons['target_pos_z'] < max_z)]
        
        connection_counts = connections['target_node_id'].value_counts()
        mean_connections = connection_counts.mean()
        std_connections = connection_counts.std()  
        cells = list(connections.drop_duplicates(subset=['target_node_id'])['target_node_id'])
        num_cells = len(cells)
        
        cell_spikes = spikes[spikes['node_ids'].isin(cells)]
        
        #skip the first few ms
        cell_spikes = cell_spikes[cell_spikes['timestamps']>skip_ms]
        spike_counts = cell_spikes.node_ids.value_counts()
        total_seconds = (ms-skip_ms)/1000
        spike_counts_per_second = spike_counts / total_seconds 

        spikes_mean = spike_counts_per_second.mean()
        spikes_std = spike_counts_per_second.std()
        
 
        location_str = str((min_x,min_y,min_z)) + "=>" + str((max_x,max_y,max_z)) 
        
        print(location_str + " : " +  str(num_cells) +" cells : {:.2f}".format(mean_connections) + " ({:.2f})".format(std_connections) +
              ": {:.2f}".format(spikes_mean) + " ({:.2f})".format(spikes_std))

        if last-spikes_mean > biggest_jump:
            biggest_jump = last-spikes_mean
            biggest_jump_str = last_str + " and " + location_str
        last = spikes_mean
        last_str = location_str

    print("Biggest jump (" + "{:.2f}".format(biggest_jump) + ") occurs between " + biggest_jump_str)


    return 0
    
def run(config):

    nodes = None
    edges = None 
    sources = ['BLA','shell']
    targets = ['BLA']
    sids = ['model_type','model_type']
    #sids = ['a_name']
    tids = ['a_name']
    prepend_pop = True
    
    #center_x = 300
    #center_y = 300
    #center_z = 300

    #iterations = 13
    #step = 25
    #print("\ttotal\tuni\tbi") 
    ret = relation_matrix(config,nodes,edges,sources,targets,sids,tids,prepend_pop,relation_func=conn_info)
    
    return

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('simulation_configECP.json')

