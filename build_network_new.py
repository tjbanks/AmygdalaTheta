from bmtk.builder import NetworkBuilder, network
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list, xiter_random
from bmtk.utils.sim_setup import build_env_bionet
import synapses
import math
import random
import os

from .connectors import (one_to_one, one_to_one_offset, syn_dist_delay_feng, syn_dist_delay_feng_section, syn_uniform_delay_section,
                        syn_percent, syn_percent_o2a, recurrent_connector, recurrent_connector_o2a)

np.random.seed(123412)

network_dir = 'network'
t_sim = 15000.0
dt = 0.05
scale = 1

edge_effects = False

#Number of cells in each population
numPN_A = 569 * scale #640 * scale #4114#15930
numPN_C = 231 * scale #260 * scale #4115#6210
numBask = 93 * scale #100 * scale #854#4860
numSOM = 51 * scale #42 * scale
numCR = 56 * scale #42 * scale
num_cells = numPN_A + numPN_C + numBask + numSOM + numCR

min_conn_dist = 0.0
max_conn_dist = 300.0 #9999.9# Distance constraint for all cells

# Create the possible x,y,z coordinates
x_start, x_end = 0,600
y_start, y_end = 0,600
z_start, z_end = 0,600
pos_list = np.random.rand(num_cells,3)
pos_list[:,0] = pos_list[:,0]*x_end - x_start
pos_list[:,1] = pos_list[:,1]*y_end - y_start
pos_list[:,2] = pos_list[:,2]*z_end - z_start


def build_networks(network_definitions):
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
    
    return networks

def build_edges(edge_definitions,edge_params,edge_add_properties,syn=None):
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

def save_networks(networks,network_dir):
    for f in os.listdir(network_dir):
        os.remove(os.path.join(network_dir, f))

    for network_name, network in enumerate(networks):
        print('Building ' + network_name)
        network.build()
        network.save_nodes(output_dir=network_dir)
        network.save_edges(output_dir=network_dir)

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

##########################################################################
############################  EDGE EFFECTS  ##############################

if edge_effects:
    core_x,core_y,core_z = (x_end-x_start),(y_end-y_start),(z_end-z_start)
    core_volume =  core_x * core_y * core_z
    shell_x_start,shell_y_start,shell_start = x_start - max_conn_dist, x_start - max_conn_dist, z_start - max_conn_dist
    shell_x_end,shell_y_end,shell_z_end = x_end + max_conn_dist, y_end + max_conn_dist, z_end + max_conn_dist
    shell_x,shell_y,shell_z = (shell_x_end-shell_x_start),(shell_y_end-shell_y_start),(shell_z_end-shell_z_start)
    shell_volume = shell_x * shell_y * shell_z
    shell_multiplier = (shell_volume/core_volume - 1) # Determine the size difference between core and shell

    virt_numPN_A = numPN_A * shell_multiplier
    virt_numPN_C = numPN_C * shell_multiplier
    virt_numBask = numBask * shell_multiplier
    virt_numSOM  = numSOM * shell_multiplier
    virt_numCR   = numCR * shell_multiplier
    virt_numCells = virt_numPN_A + virt_numPN_C + virt_numBask + virt_numSOM + virt_numCR

    # TODO - EXCLUDE POSITIONS IN THE CORE
    virt_pos_list = np.random.rand(virt_num_cells,3)
    virt_pos_list[:,0] = pos_list[:,0]*shell_x_end - shell_x_start
    virt_pos_list[:,1] = pos_list[:,1]*shell_y_end - shell_y_start
    virt_pos_list[:,2] = pos_list[:,2]*shell_z_end - shell_z_start


    edge_network = {
        'network_name':'shell',
        'positions_list':virt_pos_list,
        'cells':[
            {   # Pyramidal Cells - Type A
                'N':virt_numPN_A,
                'pop_name':'virt_PyrA',
                'a_name':'virt_PN',
                'rotation_angle_zaxis':xiter_random(N=virt_numPN_A, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=virt_numPN_A, min_x=0.0, max_x=2*np.pi),
                'model_type':'virtual'
            },
            {   # Pyramidal Cells - Type C
                'N':virt_numPN_C,
                'pop_name':'virt_PyrC',
                'a_name':'virt_PN',
                'rotation_angle_zaxis':xiter_random(N=virt_numPN_C, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=virt_numPN_C, min_x=0.0, max_x=2*np.pi),
                'model_type':'virtual'
            },
            {   # Interneuron - fast spiking PV
                'N':virt_numBask,
                'pop_name':'virt_Bask',
                'a_name':'virt_PV',
                'rotation_angle_zaxis':xiter_random(N=virt_numBask, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=virt_numBask, min_x=0.0, max_x=2*np.pi),
                'model_type':'virtual'
            },
            {   # Interneuron - SOM Cell
                'N':virt_numSOM,
                'pop_name':'virt_SOM',
                'a_name':'virt_SOM',
                'rotation_angle_zaxis':xiter_random(N=virt_numSOM, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=virt_numSOM, min_x=0.0, max_x=2*np.pi),
                'model_type':'virtual'
            },
            {   # Interneuron - CR Cell
                'N':virt_numCR,
                'pop_name':'virt_CR',
                'a_name':'virt_CR',
                'rotation_angle_zaxis':xiter_random(N=virt_numCR, min_x=0.0, max_x=2*np.pi),
                'rotation_angle_yaxis':xiter_random(N=virt_numCR, min_x=0.0, max_x=2*np.pi),
                'model_type':'virtual'
            }
        ]
    }

    networks.append(edge_network)

##########################################################################
##########################################################################


networks = build_networks(network_definitions)

int2int_temp_list = []
uncoupled_bi_track = []
pyr_int_bi_list = []


edge_definitions = [
    {   # Pyramidal to Pyramidal Connections
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        },
        'param': 'PYR2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Unidirectional
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['Bask']}
        },
        'param': 'INT2INT',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Bidirectional Pair
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['Bask']}
        },
        'param': 'INT2INT_bi_1',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Bidirectional Pair
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['Bask']}
        },
        'param': 'INT2INT_bi_2',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PYR Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        },
        'param': 'INT2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PYR to PV Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['Bask']}
        },
        'param': 'PYR2INT',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PYR Bidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        },
        'param': 'INT2PYR_bi',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PYR to PV Bidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['Bask']}
        },
        'param': 'PYR2INT_bi',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PYR to SOM Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['SOM']}
        },
        'param': 'PYR2SOM',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # SOM to PYR Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['SOM']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        },
        'param': 'SOM2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # INT to SOM Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['Bask']}, 
            'target':{'pop_name': ['SOM']}
        },
        'param': 'INT2SOM',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PYR to CR Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['CR']}
        },
        'param': 'PYR2CR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # CR to PYR Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['CR']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        },
        'param': 'CR2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # CR to PV Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['CR']}, 
            'target':{'pop_name': ['BASK']}
        },
        'param': 'CR2INT',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # CR to SOM Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['CR']}, 
            'target':{'pop_name': ['SOM']}
        },
        'param': 'CR2SOM',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },  

        ##################### VPSI INPUT #####################

    {   # VPSI Inhibition to Pyramidal
        'network':'BLA',
        'edge': {
            'source':networks['vpsi_inh'].nodes(),
            'target':networks['BLA'].nodes(pop_name=['PyrA','PyrC'])
        },
        'param': 'VPSIinh2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # VPSI Inhibition to PV
        'network':'BLA',
        'edge': {
            'source':networks['vpsi_inh'].nodes(),
            'target':networks['BLA'].nodes(pop_name=['Bask'])
        },
        'param': 'VPSIinh2INT',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },

        ##################### THALAMIC INPUT #####################

    {   # VPSI Inhibition to Pyramidal
        'network':'BLA',
        'edge': {
            'source':networks['thalamus_pyr'].nodes(),
            'target':networks['BLA'].nodes(pop_name=['PyrA','PyrC'])
        },
        'param': 'THALAMUS2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'        
    },
    {   # VPSI Inhibition to SOM
        'network':'BLA',
        'edge': {
            'source':networks['thalamus_som'].nodes(),
            'target':networks['BLA'].nodes(pop_name='SOM')
        },
        'param': 'THALAMUS2SOM',
        'add_properties': 'syn_dist_delay_feng_section_default'        
    },
    {   # VPSI Inhibition to CR
        'network':'BLA',
        'edge': {
            'source':networks['thalamus_cr'].nodes(),
            'target':networks['BLA'].nodes(pop_name='CR')
        },
        'param': 'THALAMUS2CR',
        'add_properties': 'syn_dist_delay_feng_section_default'        
    }
]

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
    },
    'PYR2INT': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.48, 'angle_dist':True, 'max_dist':max_conn_dist},#'p':0.24 before angle_dist
        'syn_weight':1,
        'dynamics_params':'PN2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'INT2PYR_bi': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.16,'track_list':pyr_int_bi_list},
        'syn_weight':1,
        'dynamics_params':'INT2PN_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PYR2INT_bi': {
        'iterator':'one_to_all',
        'connection_rule':recurrent_connector_o2a,
        'connection_params':{'p':0.3,'all_edges':pyr_int_bi_list},#was 1
        'syn_weight':1,
        'dynamics_params':'PN2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'PYR2SOM': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.618, 'angle_dist':True, 'max_dist':max_conn_dist},#0.309 before angle_dist
        'syn_weight':1,
        'dynamics_params':'PN2SOM_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'SOM2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.066},#0.066
        'syn_weight':1,
        'dynamics_params':'SOM2PN_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'INT2SOM': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.55},# Dr Unal suggested .1 -> .55 based on 7/1/21 email
        'syn_weight':1,
        'dynamics_params':'INT2SOM_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['soma']
    },
    'PYR2CR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.363, 'angle_dist':True, 'max_dist':max_conn_dist},#0.183 before angle_dist
        'syn_weight':1,
        'dynamics_params':'PN2CR_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'CR2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.116},#0.116
        'syn_weight':1,
        'dynamics_params':'CR2PN_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'CR2INT': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.297},#.297
        'syn_weight':1,
        'dynamics_params':'CR2INT_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'CR2SOM': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.764},#.764
        'syn_weight':1,
        'dynamics_params':'CR2SOM_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'VPSIinh2PYR': {
        'connection_rule'=one_to_one,
        'syn_weight'=1,
        'dynamics_params'='VPSI2PN_inh_tyler_min.json'
        'distance_range'=[0.0, 9999.9],
        'target_sections'=['basal'],
    },
    'VPSIinh2INT': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.012}, # We need aprox 10 aff to each PV
        'syn_weight':1,
        'dynamics_params'='VPSI2PV_inh_tyler_min.json',
        'distance_range'=[0.0, 9999.9],
        'target_sections':['basal']
    },
    'THALAMUS2PYR': {
        'connection_rule'=one_to_one,
        'syn_weight'=1,
        'dynamics_params'='BG2PNe_thalamus_min.json',
        'distance_range'=[0.0, 9999.9],
        'target_sections'=['basal']
    },
    'THALAMUS2SOM': {
        'connection_rule':one_to_one_offset,
        'connection_params':{'offset':numPN_A+numPN_C+numBask},
        'syn_weight':1,
        'target_sections':['basal'],
        'distance_range':[0.0, 9999.9],
        'dynamics_params':'BG2SOM_thalamus_min.json'
    },
    'THALAMUS2CR': {
        'connection_rule'=one_to_one_offset,
        'connection_params'={'offset':numPN_A+numPN_C+numBask+numSOM},
        'syn_weight'=1,
        'target_sections'=['basal'],
        'distance_range'=[0.0, 9999.9],
        'dynamics_params'='BG2CR_thalamus_min.json'
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

##########################################################################
############################  EDGE EFFECTS  ##############################

if edge_effects:
    
    virt_edges = [
        {   # Pyramidal to Pyramidal Connections
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_PyrA','virt_PyrC']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'PYR2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # PV to PV Uncoupled Unidirectional
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_Bask']}), 
                'target':{'pop_name': ['Bask']}
            },
            'param': 'INT2INT',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
            # PV to PV Uncoupled Bidirectional Pair
            # PV to PV Uncoupled Bidirectional Pair
        {   # PV to PYR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_Bask']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'INT2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # PYR to PV Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_PyrA','virt_PyrC']}), 
                'target':{'pop_name': ['Bask']}
            },
            'param': 'PYR2INT',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
            # PV to PYR Bidirectional 
            # PYR to PV Bidirectional    
        {   # PYR to SOM Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_PyrA','virt_PyrC']}), 
                'target':{'pop_name': ['SOM']}
            },
            'param': 'PYR2SOM',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # SOM to PYR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_SOM']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'SOM2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # INT to SOM Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_Bask']}), 
                'target':{'pop_name': ['SOM']}
            },
            'param': 'INT2SOM',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # PYR to CR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_PyrA','virt_PyrC']}), 
                'target':{'pop_name': ['CR']}
            },
            'param': 'PYR2CR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # CR to PYR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_CR']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'CR2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # CR to PV Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_CR']}), 
                'target':{'pop_name': ['BASK']}
            },
            'param': 'CR2INT',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # CR to SOM Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['virt_CR']}), 
                'target':{'pop_name': ['SOM']}
            },
            'param': 'CR2SOM',
            'add_properties': 'syn_dist_delay_feng_section_default'
        }
    ]

    edge_definitions = edge_definitions + virt_edges
##########################################################################
########################## END EDGE EFFECTS ##############################


##########################################################################
###############################  BUILD  ##################################

# Load synapse dictionaries
synapses.load()
syn = synapses.syn_params_dicts()

build_edges(edge_definitions,edge_params,edge_add_properties,syn)

save_networks(networks,network_dir)

build_env_bionet(base_dir='./',
		network_dir=network_dir,
		tstop=t_sim, dt = dt,
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