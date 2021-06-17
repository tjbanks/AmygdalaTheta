import sys

import numpy as np
import pandas as pd

from bmtools.cli.plugins.util.util import relation_matrix

def conn_info(**kwargs):

    edges = kwargs["edges"]
    source_id_type = kwargs["sid"]
    target_id_type = kwargs["tid"]
    source_id = kwargs["source_id"]
    target_id = kwargs["target_id"]
    t_list = kwargs["target_nodes"]
    s_list = kwargs["source_nodes"]
    
    if source_id == 'SOM' and target_id == 'PN':
        cons = edges[(edges[source_id_type] == source_id) & (edges[target_id_type]==target_id)]
        print("valid pn targets for som:")
        print(cons.target_node_id.unique())   

    return 0

def run(config):

    nodes = None
    edges = None 
    sources = ['BLA']
    targets = ['BLA']
    sids = ['a_name']
    tids = ['a_name']
    prepend_pop = True
    
    ret = relation_matrix(config,nodes,edges,sources,targets,sids,tids,prepend_pop,relation_func=conn_info)
    
    return

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('simulation_configECP.json')

