import sys

import numpy as np
import pandas as pd

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
from bmtools.cli.plugins.util.util import relation_matrix
from build_input import lognorm_fr_list

psg = None
t_sim = 15000
output_h5 = 'shell_spikes.h5'
mean_std = {
               'PyrA':[1,0.8],
               'PyrC':[1,0.8],
               'PV':[30,13],
               'SOM':[2,1],
               'CR':[20,4]
           }

def add_inputs(**kwargs):
    edges = kwargs["edges"]
    source_id_type = kwargs["sid"]
    target_id_type = kwargs["tid"]
    source_id = kwargs["source_id"]
    target_id = kwargs["target_id"]
    t_list = kwargs["target_nodes"]
    s_list = kwargs["source_nodes"]
    global psg
    cons = edges[(edges[source_id_type] == source_id) & (edges[target_id_type]==target_id)]

    node_ids = cons.source_node_id.unique() # Most important line - get all unique source ids 
    mean = mean_std[source_id][0]
    std = mean_std[source_id][1]
     
    psg.add(node_ids=node_ids,
        firing_rate=lognorm_fr_list(len(node_ids),mean,std),
        times=(0.0, t_sim/1000.0))

    return 
    

def build_shell_inputs(config='simulation_configECP.json'):

    population = 'shell'    
    
    nodes = None
    edges = None
    sources = [population]
    targets = ['BLA']
    sids = ['pop_name']
    tids = ['model_type']
    prepend_pop = True
    global psg
    psg = PoissonSpikeGenerator(population=population)

    # Relation matrix is kind of my secret formula for loading up all the edge files and sticking
    # them into a nice and easy to read pandas dataframe. The 'relation_func' is what that dataframe
    # gets sent to and called for each source/dest combination
    # In this case for each shell cell type, and every BLA cell call add_inputs
    relation_matrix(config,nodes,edges,sources,targets,sids,tids,prepend_pop,relation_func=add_inputs)

    psg.to_sonata(output_h5)

    return

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        build_shell_inputs(sys.argv[-1])
    else:
        build_shell_inputs('simulation_configECP.json')
