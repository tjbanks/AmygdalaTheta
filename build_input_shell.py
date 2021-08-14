import sys

import numpy as np
import pandas as pd

from bmtools.cli.plugins.util.util import relation_matrix
from build_inpu import lognorm_fr_list

psg = None
t_sim = 15000
output_h5 = 'shell_spikes.h5'
mean_std = {
               'PyrA':[1,1],
               'PyrC':[1,1],
               'Bask':[1,1],
               'SOM':[1,1],
               'CR':[1,1]
           }

def add_inputs(**kwargs):
    edges = kwargs["edges"]
    source_id_type = kwargs["sid"]
    target_id_type = kwargs["tid"]
    source_id = kwargs["source_id"]
    target_id = kwargs["target_id"]
    t_list = kwargs["target_nodes"]
    s_list = kwargs["source_nodes"]

    cons = edges[(edges[source_id_type] == source_id) & (edges[target_id_type]==target_id)]
    total_cons = cons.count().source_node_id

    
    #psg.add(node_ids=node_ids,
    #    firing_rate=lognorm_fr_list(len(node_ids),mean,std),
    #    times=(0.0, t_sim/1000.0))

    return 
    

def run(config):

    population = 'shell'    
    
    nodes = None
    edges = None
    sources = [population]
    targets = ['BLA']
    sids = ['pop_name']
    tids = ['pop_name']
    prepend_pop = True

    psg = PoissonSpikeGenerator(population=population)

    relation_matrix(config,nodes,edges,sources,targets,sids,tids,prepend_pop,relation_func=add_inputs)

    psg.to_sonata(output_h5)

    return

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('simulation_configECP.json')
