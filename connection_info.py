import sys

import numpy as np
import pandas as pd

from bmtools.cli.plugins.util.util import relation_matrix

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
    
    if edges is None:
        edges = _edges
    else:
        edges = edges.append(_edges).drop_duplicates()

    cons = edges[(edges[source_id_type] == source_id) & (edges[target_id_type]==target_id)]
    total_cons = cons.count().source_node_id
    
    # to determine reciprocal connectivity
    # create a copy and flip source/dest
    cons_flip = edges[(edges[source_id_type] == target_id) & (edges[target_id_type]==source_id)]
    cons_flip = cons_flip.rename(columns={'source_node_id':'target_node_id','target_node_id':'source_node_id'})
    # append to original 
    cons_recip = cons.append(cons_flip)

    # determine dropped duplicates (keep=False)
    cons_recip_dedup = cons_recip.drop_duplicates(subset=['source_node_id','target_node_id'])

    # note counts
    num_bi = (cons_recip.count().source_node_id - cons_recip_dedup.count().source_node_id)
    num_uni = total_cons - num_bi    

    num_sources = s_list.apply(pd.Series.value_counts)[source_id_type].dropna().sort_index().loc[source_id]
    num_targets = t_list.apply(pd.Series.value_counts)[target_id_type].dropna().sort_index().loc[target_id]

    total = round(total_cons / (num_sources*num_targets) * 100,2)
    uni = round(num_uni / (num_sources*num_targets) * 100,2)
    bi = round(num_bi / (num_sources*num_targets) * 100,2)

    print(str(source_id) + '->' + str(target_id) + "\t" + str(total) + "\t" + str(uni) + "\t" + str(bi))
    #if source_id == 'PN' and target_id == 'PN':
        #import pdb;pdb.set_trace()
        #print(cons)   
    return total

def run(config):

    nodes = None
    edges = None 
    sources = ['BLA','shell']
    targets = ['BLA']
    sids = ['a_name','a_name']
    tids = ['a_name']
    prepend_pop = True
    
    print("\ttotal\tuni\tbi") 
    ret = relation_matrix(config,nodes,edges,sources,targets,sids,tids,prepend_pop,relation_func=conn_info)
    
    return

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('simulation_configECP.json')

