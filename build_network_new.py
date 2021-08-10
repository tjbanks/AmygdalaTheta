from bmtk.builder import NetworkBuilder, network
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list, xiter_random
import synapses
import math
import pdb
import random
import pandas as pd
import os

from .connectors import (one_to_one, one_to_one_offset, syn_dist_delay_feng, syn_dist_delay_feng_section, syn_uniform_delay_section,
                        syn_percent, syn_percent_o2a, recurrent_connector, recurrent_connector_o2a)

np.random.seed(123412)

# Initialize our network
net = NetworkBuilder("BLA")

scale = 1

#Number of cells in each population
numPN_A = 569 * scale #640 * scale #4114#15930
numPN_C = 231 * scale #260 * scale #4115#6210
numBask = 93 * scale #100 * scale #854#4860
numSOM = 51 * scale #42 * scale
numCR = 56 * scale #42 * scale
# add_properties = False
# do_pos = False
num_cells = numPN_A + numPN_C + numBask + numSOM + numCR

dist_constraint = True 
min_conn_dist = 0.0
max_conn_dist = 9999.9

if dist_constraint:
    max_conn_dist = 300.0 #9999.9# Distance constraint for all cells

i2i_gap = False

# Create the possible x,y,z coordinates
x_start, x_end = 0,600
y_start, y_end = 0,600
z_start, z_end = 0,600
pos_list = np.random.rand(num_cells,3)
pos_list[:,0] = pos_list[:,0]*x_end - x_start
pos_list[:,1] = pos_list[:,1]*y_end - y_start
pos_list[:,2] = pos_list[:,2]*z_end - z_start

networks = {} #Place to store NetworkBuilder objects referenced by name
network_definitions = [
    {
        'network_name':'BLA',
        'positions_list':pos_list,
        'cells':[
            {   # Pyramidal Cells - Type A
                'N':numPN_A,
                'pop_name':'PyrA',
                'a_name':'PN',
                'rotation_angle_zaxis':xiter_random(N=numPN_A, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=numPN_A, min_x=0.0, max_x=2*np.pi),
                'model_type':'biophysical',
                'model_template':'hoc:Cell_Af'
            },
            {   # Pyramidal Cells - Type C
                'N':numPN_C,
                'pop_name':'PyrC',
                'a_name':'PN',
                'rotation_angle_zaxis':xiter_random(N=numPN_C, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=numPN_C, min_x=0.0, max_x=2*np.pi),
                'model_type':'biophysical',
                'model_template':'hoc:Cell_Cf'
            },
            {   # Interneuron - fast spiking PV
                'N':numBask,
                'pop_name':'Bask',
                'a_name':'PV',
                'rotation_angle_zaxis':xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
                'model_type':'biophysical',
                'model_template':'hoc:InterneuronCellf'
            },
            {   # Interneuron - SOM Cell
                'N':numSOM,
                'pop_name':'SOM',
                'a_name':'SOM',
                'rotation_angle_zaxis':xiter_random(N=numSOM, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=numSOM, min_x=0.0, max_x=2*np.pi),
                'model_type':'biophysical',
                'model_template':'hoc:SOM_Cell'
            },
            {   # Interneuron - CR Cell
                'N':numCR,
                'pop_name':'CR',
                'a_name':'CR',
                'rotation_angle_zaxis':xiter_random(N=numCR, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=numCR, min_x=0.0, max_x=2*np.pi),
                'model_type':'biophysical',
                'model_template':'hoc:CR_Cell'
            }
        ] # End cells
    }, # End BLA
    {   # VPSI INPUTS PYR
        'network_name':'vpsi_pyr',
        'positions_list':None,
        'cells':[
            {
                'N':numPN_A+numPN_C,
                'pop_name':'pyr_inp',
                'pop_group':'vpsi_pyr',
                'model_type':'virtual'
            }
        ]
    },
    {   # VPSI INPUTS PV
        'network_name':'vpsi_pv',
        'positions_list':None,
        'cells':[
            {
                'N':numBask,
                'pop_name':'pv_inp',
                'pop_group':'vpsi_pv',
                'model_type':'virtual'
            }
        ]
    },
    {   # VPSI INPUTS Inhibitory
        'network_name':'vpsi_inh',
        'positions_list':None,
        'cells':[
            {
                'N':numPN_A+numPN_C+numBask,
                'pop_name':'inh_inp',
                'pop_group':'vpsi_inh',
                'model_type':'virtual'
            }
        ]
    },
    {
        # Thalamic PYR INPUTS
        'network_name':'thalamus_pyr',
        'positions_list':None,
        'cells':[
            {
                'N':numPN_A+numPN_C,
                'pop_name':'pyr_inp',
                'pop_group':'thalamus_pyr',
                'model_type':'virtual'
            }
        ]
    },
    {
        # Thalamic PV INPUTS
        'network_name':'thalamus_pv',
        'positions_list':None,
        'cells':[
            {
                'N':numBask,
                'pop_name':'pv_inp',
                'pop_group':'thalamus_pv',
                'model_type':'virtual'
            }
        ]
    },
    {
        # Thalamic SOM INPUTS
        'network_name':'thalamus_som',
        'positions_list':None,
        'cells':[
            {
                'N':numSOM,
                'pop_name':'som_inp',
                'pop_group':'thalamus_som',
                'model_type':'virtual'
            }
        ]
    },
    {
        # Thalamic CR INPUTS
        'network_name':'thalamus_cr',
        'positions_list':None,
        'cells':[
            {
                'N':numCR,
                'pop_name':'cr_inp',
                'pop_group':'thalamus_cr',
                'model_type':'virtual'
            }
        ]
    }
]

int2int_temp_list = []
uncoupled_bi_track = []

edge_params = {
    'PYR2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.04, 'angle_dist':True, 'max_dist':max_conn_dist},
        'syn_weight':1,
        'dynamics_params':'PN2PN_feng_min.json',
        'distance_range':[0,max_conn_dist],
        'target_sections':['basal']
    },
    'INT2INT': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.19,'no_recip':True,'track_list':int2int_temp_list},#0.19
        'syn_weight':1,
        'dynamics_params':'INT2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'INT2INT_bi_1': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.025, 'track_list':uncoupled_bi_track},#0.03
        'syn_weight':1,
        'dynamics_params':'INT2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'INT2INT_bi_2': {
        'iterator':'one_to_all',
        'connection_rule':recurrent_connector_o2a,
        'connection_params':{'p':1, 'all_edges':uncoupled_bi_track},#p:1
        'syn_weight':1,
        'dynamics_params':'INT2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'INT2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.36},#{'p':0.40},
        'syn_weight':1,
        'dynamics_params':'INT2PN_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    }
} # edges referenced by name

edge_add_properties = {
    'syn_dist_delay_feng_section_default': {
        'names':['delay','sec_id','sec_x'],
        'rule':syn_dist_delay_feng_section,
        'rule_params':{'sec_id':1, 'sec_x':0.9},
        'dtypes':[np.float, np.int32, np.float]
    }
}

edge_definitions = [
    {   # Pyramidal to Pyramidal Connections
        'network':'BLA'
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        }
        'param': 'PYR2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Unidirectional
        'network':'BLA'
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['Bask']}
        }
        'param': 'INT2INT',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Bidirectional Pair
        'network':'BLA'
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['Bask']}
        }
        'param': 'INT2INT_bi_1',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Bidirectional Pair
        'network':'BLA'
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['Bask']}
        }
        'param': 'INT2INT_bi_2',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PYR Unidirectional 
        'network':'BLA'
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        }
        'param': 'INT2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    }

]



##########################################################################
############################### PYR2INT ##################################
########################### UNIDIRECTIONAL ###############################

# Create connections between Pyr --> Bask cells

if connect["PYR2INT"]:
    #dynamics_file = 'PN2INT.json'
    dynamics_file = 'PN2INT_feng_min.json'

    conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'Bask'},
                iterator = 'one_to_all',
                connection_rule=syn_percent_o2a,
                connection_params={'p':0.48, 'angle_dist':True, 'max_dist':max_conn_dist},#'p':0.24 before angle_dist
                syn_weight=1,
                dynamics_params=dynamics_file,
                model_template=syn[dynamics_file]['level_of_detail'],
                distance_range=[min_conn_dist,max_conn_dist],
                target_sections=['basal'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
                rule=syn_dist_delay_feng_section,
                rule_params={'sec_id':1, 'sec_x':0.9},
                dtypes=[np.float, np.int32, np.float])

##########################################################################
############################### PYR2INT ##################################
############################ BIDIRECTIONAL ###############################

pyr_int_bi_list = []

if connect["INT2PYR"]:
    dynamics_file = 'INT2PN_feng_min.json'
    conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['PyrA','PyrC']},
                iterator = 'one_to_all',
                connection_rule=syn_percent_o2a,
                connection_params={'p':0.16,'track_list':pyr_int_bi_list},
                syn_weight=1,
                dynamics_params=dynamics_file,
                model_template=syn[dynamics_file]['level_of_detail'],
                distance_range=[min_conn_dist,max_conn_dist],
                target_sections=['somatic'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
                rule=syn_dist_delay_feng_section,
                rule_params={'sec_id':0, 'sec_x':0.9},
                dtypes=[np.float, np.int32, np.float])


if connect["PYR2INT"]:
    dynamics_file = 'PN2INT_feng_min.json'
    conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'Bask'},
                iterator = 'one_to_all',
                connection_rule=recurrent_connector_o2a,
                connection_params={'p':0.3,'all_edges':pyr_int_bi_list},#was 1
                syn_weight=1,
                dynamics_params=dynamics_file,
                model_template=syn[dynamics_file]['level_of_detail'],
                distance_range=[min_conn_dist,max_conn_dist],
                target_sections=['basal'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
                rule=syn_dist_delay_feng_section,
                rule_params={'sec_id':1, 'sec_x':0.9},
                dtypes=[np.float, np.int32, np.float])


"""
dynamics_file = 'INT2PN_feng_min.json'
conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_all',
              connection_rule=recurrent_connector_o2a,
              connection_params={'p':1,'all_edges':pyr_int_bi_list},
              syn_weight=1,
              delay = 0.1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['somatic'],
              sec_id=0,
              sec_x=0.9)

conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':0, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])

"""
##########################################################################
############################### PYR2SOM ##################################


if connect["PYR2SOM"]:
    dynamics_file = 'PN2SOM_tyler.json'
    conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': ['SOM']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.618, 'angle_dist':True, 'max_dist':max_conn_dist},#0.309 before angle_dist
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['basal'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':1, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])


##########################################################################
############################### SOM2PYR ##################################

if connect["SOM2PYR"]:
    dynamics_file = 'SOM2PN_tyler.json'
    conn = net.add_edges(source={'pop_name': ['SOM']}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.066},#0.066
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['somatic'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':0, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])

##########################################################################
############################### INT2SOM ##################################


if connect["INT2SOM"]:
    dynamics_file = 'INT2SOM_tyler.json'
    conn = net.add_edges(source={'pop_name': ['Bask']}, target={'pop_name': ['SOM']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.55},# Dr Unal suggested .1 -> .55 based on 7/1/21 email
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['soma'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':0, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])



##########################################################################
############################### PYR2CR  ##################################


if connect["PYR2CR"]:
    dynamics_file = 'PN2CR_tyler.json'
    conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': ['CR']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.363, 'angle_dist':True, 'max_dist':max_conn_dist},#0.183 before angle_dist
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['basal'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':1, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])


##########################################################################
############################### CR2PYR  ##################################

if connect["CR2PYR"]:
    dynamics_file = 'CR2PN_tyler.json'
    conn = net.add_edges(source={'pop_name': ['CR']}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.116},#0.116
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['somatic'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':0, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])


##########################################################################
############################### CR2INT  ##################################

if connect["CR2INT"]:
    dynamics_file = 'CR2INT_tyler.json'
    conn = net.add_edges(source={'pop_name': ['CR']}, target={'pop_name': ['Bask']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.297},#.297
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['basal'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':1, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])

##########################################################################
############################### CR2SOM  ##################################

if connect["CR2SOM"]:
    dynamics_file = 'CR2SOM_tyler.json'
    conn = net.add_edges(source={'pop_name': ['CR']}, target={'pop_name': ['SOM']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.764},#.764
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['basal'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':1, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])

##########################################################################
######################### BACKGROUND INPUT ###############################

############################ VPSI INPUT ##################################

if connect["VPSIinh2PYR"]:
    dynamics_file = 'VPSI2PN_inh_tyler_min.json'

    conn = net.add_edges(source=vpsi_inh.nodes(), target=net.nodes(pop_name=['PyrA']),
                   connection_rule=one_to_one,
                   syn_weight=1,
                   target_sections=['basal'],
                   distance_range=[0.0, 9999.9],
                   dynamics_params=dynamics_file,
                   model_template=syn[dynamics_file]['level_of_detail'])


    conn.add_properties(names=['delay','sec_id','sec_x'],
                  rule=syn_uniform_delay_section,
                  rule_params={'sec_id':1, 'sec_x':0.9},
                  dtypes=[np.float, np.int32, np.float])

    conn = net.add_edges(source=vpsi_inh.nodes(), target=net.nodes(pop_name=['PyrC']),
                   connection_rule=one_to_one,
                   syn_weight=1,
                   target_sections=['basal'],
                   distance_range=[0.0, 9999.9],
                   dynamics_params=dynamics_file,
                   model_template=syn[dynamics_file]['level_of_detail'])


    conn.add_properties(names=['delay','sec_id','sec_x'],
                  rule=syn_uniform_delay_section,
                  rule_params={'sec_id':1, 'sec_x':0.9},
                  dtypes=[np.float, np.int32, np.float])

if connect["VPSIinh2INT"]:
    dynamics_file = 'VPSI2PV_inh_tyler_min.json'

    conn = net.add_edges(source=vpsi_inh.nodes(), target=net.nodes(pop_name='Bask'),
                   iterator='one_to_all',
                   connection_rule=syn_percent_o2a,
                   connection_params={'p':0.012}, # We need aprox 10 aff to each PV
                   syn_weight=1,
                   target_sections=['basal'],
                   distance_range=[0.0, 9999.9],
                   dynamics_params=dynamics_file,
                   model_template=syn[dynamics_file]['level_of_detail'])


    conn.add_properties(names=['delay','sec_id','sec_x'],
                  rule=syn_uniform_delay_section,
                  rule_params={'sec_id':1, 'sec_x':0.9},
                  dtypes=[np.float, np.int32, np.float])

######################### THALAMIC INPUT ###############################

if connect["THALAMUS2PYR"]:
        
    dynamics_file='BG2PNe_thalamus_min.json'

    conn = net.add_edges(source=thalamus_pyr.nodes(), target=net.nodes(pop_name='PyrA'),
                    connection_rule=one_to_one,
                    syn_weight=1,
                    target_sections=['basal'],
                    distance_range=[0.0, 9999.9],
                    dynamics_params=dynamics_file,
                    model_template=syn[dynamics_file]['level_of_detail'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
                    rule=syn_uniform_delay_section,
                    rule_params={'sec_id':1, 'sec_x':0.9},
                    dtypes=[np.float, np.int32, np.float])

    conn = net.add_edges(source=thalamus_pyr.nodes(), target=net.nodes(pop_name='PyrC'),
                    connection_rule=one_to_one,
                    syn_weight=1,
                    target_sections=['basal'],
                    distance_range=[0.0, 9999.9],
                    dynamics_params=dynamics_file,
                    model_template=syn[dynamics_file]['level_of_detail'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
                    rule=syn_uniform_delay_section,
                    rule_params={'sec_id':1, 'sec_x':0.9},
                    dtypes=[np.float, np.int32, np.float])


if connect["THALAMUS2SOM"]:

    dynamics_file = 'BG2SOM_thalamus_min.json'

    conn = net.add_edges(source=thalamus_som.nodes(), target=net.nodes(pop_name='SOM'),
                    connection_rule=one_to_one_offset,
                    connection_params={'offset':numPN_A+numPN_C+numBask},
                    syn_weight=1,
                    target_sections=['basal'],
                    distance_range=[0.0, 9999.9],
                    dynamics_params=dynamics_file,
                    model_template=syn[dynamics_file]['level_of_detail'])


    conn.add_properties(names=['delay','sec_id','sec_x'],
                    rule=syn_uniform_delay_section,
                    rule_params={'sec_id':1, 'sec_x':0.9},
                    dtypes=[np.float, np.int32, np.float])

if connect["THALAMUS2CR"]:

    dynamics_file = 'BG2CR_thalamus_min.json'

    conn = net.add_edges(source=thalamus_cr.nodes(), target=net.nodes(pop_name='CR'),
                    connection_rule=one_to_one_offset,
                    connection_params={'offset':numPN_A+numPN_C+numBask+numSOM},
                    syn_weight=1,
                    target_sections=['basal'],
                    distance_range=[0.0, 9999.9],
                    dynamics_params=dynamics_file,
                    model_template=syn[dynamics_file]['level_of_detail'])


    conn.add_properties(names=['delay','sec_id','sec_x'],
                    rule=syn_uniform_delay_section,
                    rule_params={'sec_id':1, 'sec_x':0.9},
                    dtypes=[np.float, np.int32, np.float])

##########################################################################
###############################  BUILD  ##################################


# Build the networks
for net_def in network_definitions:
    network_name = net_def['network_name']
    networks[network_name] = NetworkBuilder(network_name)
    pos_list = net_def['positions_list']

    # Add cells to the network
    for cell in net_def['cells']:
        num_cells = cell['N']
        extra_kwargs = {}
        if pos_list:
            inds = np.random.choice(np.arange(0,np.size(pos_list,0)),num_cells,replace=False)
            pos = pos_list[inds,:]
            # Get rid of coordinates already used
            pos_list = np.delete(pos_list,inds,0)
            extra_kwargs['positions'] = positions_list(positions=pos)

        networks[network_name].add_nodes(**cell,**extra_kwargs)

# Load synapse dictionaries
synapses.load()
syn = synapses.syn_params_dicts()

# Build the edges
for edge in edge_definitions:
    network_name = edge['network']
    edge_src_trg = edge['edge']
    edge_params  = edge_params[edge['param']]
    dynamics_file = edge_params['dynamics_file']
    model_template = syn[dynamics_file]['level_of_detail']

    model_template_kwarg = {'model_template':model_template}

    net = networks[network_name]

    conn = net.add_edges(**edge_src_trg,**edge_params,**model_template_kwarg)
    
    if edge.get('add_properties'):
        edge_add_properties = edge['add_properties']
        conn.add_properties(**edge_add_properties)


network_dir = 'network'
for f in os.listdir(network_dir):
    os.remove(os.path.join(network_dir, f))

net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

print("Internal nodes and edges built")

# Create connections between "vpsi_pyr" and Pyramidals
# First define the connection rule

# Build and save our network

vpsi_pyr.build()
vpsi_pyr.save_nodes(output_dir='network')

vpsi_pv.build()
vpsi_pv.save_nodes(output_dir='network')

vpsi_inh.build()
vpsi_inh.save_nodes(output_dir='network')

thalamus_pyr.build()
thalamus_pyr.save_nodes(output_dir='network')

thalamus_pv.build()
thalamus_pv.save_nodes(output_dir='network')

thalamus_som.build()
thalamus_som.save_nodes(output_dir='network')

thalamus_cr.build()
thalamus_cr.save_nodes(output_dir='network')
#
#print("External nodes and edges built")
t_sim = 15000.0

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
		network_dir='./network',
		tstop=t_sim, dt = 0.05,
		report_vars = ['v'],
                v_init = -70.0,
                celsius = 31.0,
		spikes_inputs=[('vpsi_pyr','vpsi_pyr_spikes.h5'),  # Name of population which spikes will be generated for, file
                       ('vpsi_pv','vpsi_pv_spikes.h5'),
                       ('vpsi_inh','vpsi_inh_spikes.h5'),
                       ('thalamus_pyr','thalamus_pyr_spikes.h5'),
                       ('thalamus_pv', 'thalamus_pv_spikes.h5'),
                       ('thalamus_som','thalamus_som_spikes.h5'),
                       ('thalamus_cr','thalamus_cr_spikes.h5'),
                       ],
		components_dir='components',
                config_file='simulation_config.json',
		compile_mechanisms=True)


