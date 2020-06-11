"""
Date: 6/11/2020
"""
import cell_positions as p
import synapses

from bmtk.builder.auxi.edge_connectors import distance_connector, connect_random
from bmtk.builder.networks import NetworkBuilder
import math
import numpy as np
import random

def build(output='network', seed=None, experiment="0"):
    
    seed_default = 2020 

    #if seed is not None:#Not variable anymore to prevent bad networks (Feng recommended)
    #    seed_default = seed

    random.seed(seed_default)
    #np.random.seed(seed_default)

    #initialize the networks
    amygdala_net = NetworkBuilder('amygdala')
    input_net = NetworkBuilder('input')
    
    ###########################################################
    # Build node locations
    ###########################################################
    
    #define cell numbers
    

    #get cell positions
    

    #verify cell numbers match cell position list sizes
    

    ###########################################################
    # Experiment Logic
    ###########################################################

    build_node = {}
    build_edge = {}
    nodes = []
    edges = []

    # NODES TO BE BUILT
    
    
    #EDGES TO BE BUILT
    

    #EXPERIMENTS - Add the nodes and edges you want for each experiment
    

    #MERGE -- nodes and edges specified will be built
    build_node.update(dict.fromkeys(nodes, True))
    build_edge.update(dict.fromkeys(edges, True))

    ###########################################################
    # Build nodes
    ###########################################################

    #Add noteds to the network
    

    ###########################################################
    # Build custom synapses
    ###########################################################
    #See https://github.com/AllenInstitute/bmtk/blob/develop/bmtk/simulator/bionet/default_setters/synapse_models.py

    synapses.load()
    syn = synapses.syn_params_dicts()
    syn_list = []

    ###########################################################
    # Build custom connection rules
    ###########################################################
    #See bmtk.builder.auxi.edge_connectors
    

    ###########################################################
    # Build individual edge properties
    ###########################################################
    #Individual edge properties (See bmtk.docs.tutorials.NetworkBuilder_Intro.ipynb)

    
    ###########################################################
    # Build connections
    ###########################################################
    
    ###########################################################
    # Build recurrent connection rules
    ###########################################################
 

    ########################################################### 
    # Build recurrent connections
    ###########################################################
  
    ###########################################################
    # Build strict connection rules
    ###########################################################
 

    ########################################################### 
    # Build strict connections
    ###########################################################

  
    ###########################################################
    # Create External Networks
    ###########################################################
 
    ###########################################################
    # Theta Network Background Noise
    ###########################################################



    ###########################################################
    # Build networks
    ###########################################################
    amygdala_net.build()
    amygdala_net.save(output_dir=output)    
    
    input_net.build()
    input_net.save(output_dir=output)
    
if __name__ == "__main__":
    build()
