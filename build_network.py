import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import cProfile
import pstats
import numpy as np
import synapses
import math
import random
import os, sys

from bmtk.builder.auxi.node_params import positions_cuboid, positions_list, xiter_random
from bmtk.utils.sim_setup import build_env_bionet
from bmtk.builder import NetworkBuilder

from connectors import (simple_rand, one_to_one, one_to_one_offset, syn_dist_delay_feng_section, syn_uniform_delay_section,
                        syn_percent_o2a, recurrent_connector_o2a)

from connectors import init_connectors

profiler = cProfile.Profile()
profiler.enable()

np.random.seed(123412)

network_dir = 'network'
components_dir = 'components'
if os.environ.get("TSTOP"):
    t_sim = float(os.environ.get("TSTOP"))
else:
    t_sim = 15000.0
print(f'TSTOP set to {t_sim}')
dt = 0.1
scale = 1

min_conn_dist = 0.0
max_conn_dist = 9999.0 #300.0 #300.0 #9999.9# Distance constraint for all cells

# When enabled, a shell of virtual cells will be created around the core network.
edge_effects = False 
net_size = 1400#um]=[]

# Use gap junctions?
use_gap_junctions = False

#Number of cells in each population
numPN_A = 569
numPN_C = 231
numPV = 93
numSOM = 51
numCR = 56

numVPSI = 100

if __name__ == '__main__':
    if 'feng_homogenous' in sys.argv:
        numPN_A = 640
        numPN_C = 260
        numPV = 100
        numSOM = 0
        numCR = 0

        network_dir = network_dir + '_feng_homogenous'
        components_dir = components_dir + '_feng_homogenous'
        scale = 1
        max_conn_dist = 9999.9
        edge_effects = False
        net_size = 600

        print("Building Feng's original homogenous network")
        print("Ommiting SOM and CR cells")

    elif 'homogenous' in sys.argv:
        network_dir = network_dir + '_homogenous'
        components_dir = components_dir + '_homogenous'
        scale = 1
        max_conn_dist = 9999.9
        edge_effects = False
        net_size = 600

        print('Building homogenous network')

    else:
        print('Building full network')

#Scale the number of cells in each population
numPN_A = int(numPN_A * scale) #640 * scale #4114#15930
numPN_C = int(numPN_C * scale) #260 * scale #4115#6210
numPV = int(numPV * scale) #100 * scale #854#4860
numSOM = int(numSOM * scale) #42 * scale
numCR = int(numCR * scale) #42 * scale
num_cells = numPN_A + numPN_C + numPV + numSOM + numCR #Only used to populate an overall position list
numVPSI = int(numVPSI * scale)
num_shell_cells = 0

# Create the possible x,y,z coordinates
x_start, x_end = 0+max_conn_dist,net_size+max_conn_dist
y_start, y_end = 0+max_conn_dist,net_size+max_conn_dist
z_start, z_end = 0+max_conn_dist,net_size+max_conn_dist

if not edge_effects:
    x_start, x_end = 0,net_size
    y_start, y_end = 0,net_size
    z_start, z_end = 0,net_size    

pos_list = np.random.rand(num_cells,3)
pos_list[:,0] = pos_list[:,0]*(x_end - x_start) + x_start
pos_list[:,1] = pos_list[:,1]*(y_end - y_start) + y_start
pos_list[:,2] = pos_list[:,2]*(z_end - z_start) + z_start


def build_networks(network_definitions: list) -> dict:
    # network_definitions should be a list of dictionaries eg:[{}]
    # Keys should include an arbitrary 'network_name', a positions_list (if any),
    # And 'cells'. 'cells' should contain a list of dictionaries, and the dictionary 
    # should corrospond with any valid input for BMTK's NetworkBuilder.add_nodes method 
    # A dictionary of NetworkBuilder BMTK objects will be returned, reference by individual network_name
    for net_def in network_definitions:
        network_name = net_def['network_name']
        networks[network_name] = NetworkBuilder(network_name)
        pos_list = net_def.get('positions_list',None)

        # Add cells to the network
        for cell in net_def['cells']:
            num_cells = cell['N']
            extra_kwargs = {}
            if pos_list is not None:
                inds = np.random.choice(np.arange(0,np.size(pos_list,0)),num_cells,replace=False)
                pos = pos_list[inds,:]
                # Get rid of coordinates already used
                pos_list = np.delete(pos_list,inds,0)
                extra_kwargs['positions'] = positions_list(positions=pos)

            networks[network_name].add_nodes(**cell,**extra_kwargs)
    
    return networks

def build_edges(networks,edge_definitions,edge_params,edge_add_properties,syn=None,use_gap_junctions=True):
    # Builds the edges for each network given a set of 'edge_definitions'
    # edge_definitions examples shown later in the code
    for edge in edge_definitions:
        network_name = edge['network']
        edge_src_trg = edge['edge']
        edge_params_val  = edge_params[edge['param']]

        net = networks[network_name]

        if edge_params_val.get('is_gap_junction'):
            if use_gap_junctions: # we to add gap
                del edge_params_val['is_gap_junction']
                conn = net.add_gap_junctions(**edge_src_trg,**edge_params_val)
            continue # we're not using gaps or we added, continue

        model_template_kwarg = {}
        if edge_params_val.get('dynamics_params'):
            dynamics_file = edge_params_val['dynamics_params']
            model_template = syn[dynamics_file]['level_of_detail']
            model_template_kwarg = {'model_template':model_template}

        conn = net.add_edges(**edge_src_trg,**edge_params_val,**model_template_kwarg)
        
        if edge.get('add_properties'):
            edge_add_properties_val = edge_add_properties[edge['add_properties']]
            conn.add_properties(**edge_add_properties_val)

def save_networks(networks,network_dir):
    # Remove the existing network_dir directory
    for f in os.listdir(network_dir):
        os.remove(os.path.join(network_dir, f))

    # Run through each network and save their nodes/edges
    for i, (network_name, network) in enumerate(networks.items()):
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
                'rotation_angle_zaxis':xiter_random(N=numPN_A, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=numPN_A, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'biophysical',
                'model_template':'hoc:Cell_Af'
            },
            {   # Pyramidal Cells - Type C
                'N':numPN_C,
                'pop_name':'PyrC',
                'a_name':'PN',
                'rotation_angle_zaxis':xiter_random(N=numPN_C, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=numPN_C, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'biophysical',
                'model_template':'hoc:Cell_Cf'
            },
            {   # Interneuron - fast spiking PV
                'N':numPV,
                'pop_name':'PV',
                'a_name':'PV',
                'rotation_angle_zaxis':xiter_random(N=numPV, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=numPV, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'biophysical',
                'model_template':'hoc:InterneuronCellf'
            },
            {   # Interneuron - SOM Cell
                'N':numSOM,
                'pop_name':'SOM',
                'a_name':'SOM',
                'rotation_angle_zaxis':xiter_random(N=numSOM, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=numSOM, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'biophysical',
                'model_template':'hoc:SOM_Cell'
            },
            {   # Interneuron - CR Cell
                'N':numCR,
                'pop_name':'CR',
                'a_name':'CR',
                'rotation_angle_zaxis':xiter_random(N=numCR, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=numCR, min_x=np.pi/2, max_x=np.pi/2),
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
                'N':numVPSI,
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
                'N':numPV,
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

if edge_effects: # When enabled, a shell of virtual cells will be created around the core network.
    
    # compute the core volume
    core_x,core_y,core_z = (x_end-x_start),(y_end-y_start),(z_end-z_start)
    core_volume =  core_x * core_y * core_z

    # compute the outer shell volume. The absolute max_conn_dist will extend each dimension of the core by 2*max_conn_dist
    shell_x_start,shell_y_start,shell_z_start = x_start - max_conn_dist, x_start - max_conn_dist, z_start - max_conn_dist
    shell_x_end,shell_y_end,shell_z_end = x_end + max_conn_dist, y_end + max_conn_dist, z_end + max_conn_dist
    shell_x,shell_y,shell_z = (shell_x_end-shell_x_start),(shell_y_end-shell_y_start),(shell_z_end-shell_z_start)
    shell_volume = shell_x * shell_y * shell_z

    # Determine the size difference between core and shell
    shell_multiplier = (shell_volume/core_volume) 

    # Increase the number of original cells based on the shell_multiplier
    virt_numPN_A = int(numPN_A * shell_multiplier)
    virt_numPN_C = int(numPN_C * shell_multiplier)
    virt_numPV = int(numPV * shell_multiplier)
    virt_numSOM  = int(numSOM * shell_multiplier)
    virt_numCR   = int(numCR * shell_multiplier)
    virt_num_cells = virt_numPN_A + virt_numPN_C + virt_numPV + virt_numSOM + virt_numCR
    
    # Create a positions list for each cell in the shell, this includes positions in the core
    virt_pos_list = np.random.rand(virt_num_cells,3)
    virt_pos_list[:,0] = virt_pos_list[:,0]*(shell_x_end - shell_x_start) + shell_x_start
    virt_pos_list[:,1] = virt_pos_list[:,1]*(shell_y_end - shell_y_start) + shell_y_start
    virt_pos_list[:,2] = virt_pos_list[:,2]*(shell_z_end - shell_z_start) + shell_z_start
    
    # EXCLUDE POSITIONS IN THE CORE - We remove all virtual cells located in the core (accounting for no -1 on shell_multiplier)
    in_core = np.where(((virt_pos_list[:,0] > x_start) & (virt_pos_list[:,0] < x_end)) & 
                       ((virt_pos_list[:,1] > y_start) & (virt_pos_list[:,1] < y_end)) & 
                       ((virt_pos_list[:,2] > z_start) & (virt_pos_list[:,2] < z_end)))
    virt_pos_list = np.delete(virt_pos_list,in_core,0)

    # Bring down the number of shell cells to create by scaling
    # This ensures we have enough positions in virt_pos_list for all of our cells
    # Old density multiplied by new number of cells
    new_virt_num_cells = len(virt_pos_list)
    virt_numPN_A = int(virt_numPN_A/virt_num_cells*new_virt_num_cells)
    virt_numPN_C = int(virt_numPN_C/virt_num_cells*new_virt_num_cells)
    virt_numPV = int(virt_numPV/virt_num_cells*new_virt_num_cells)
    virt_numSOM = int(virt_numSOM/virt_num_cells*new_virt_num_cells)
    virt_numCR = int(virt_numCR/virt_num_cells*new_virt_num_cells)
    virt_num_cells = virt_numPN_A + virt_numPN_C + virt_numPV + virt_numSOM + virt_numCR

    # This should always be true, virt_num_cells is now equal to a scaled down number
    # While new_virt_num_cells is the length of the available cells
    assert(virt_num_cells <= new_virt_num_cells)    

    # This network should contain all the same properties as the original network, except
    # the cell should be virtual. For connectivity, you should name the cells the same as
    # the original network because connection rules defined later will require it
    shell_network = {
        'network_name':'shell',
        'positions_list':virt_pos_list,
        'cells':[
            {   # Pyramidal Cells - Type A
                'N':virt_numPN_A,
                'pop_name':'PyrA',
                'a_name':'PN',
                'rotation_angle_zaxis':xiter_random(N=virt_numPN_A, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=virt_numPN_A, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'virtual'
            },
            {   # Pyramidal Cells - Type C
                'N':virt_numPN_C,
                'pop_name':'PyrC',
                'a_name':'PN',
                'rotation_angle_zaxis':xiter_random(N=virt_numPN_C, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=virt_numPN_C, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'virtual'
            },
            {   # Interneuron - fast spiking PV
                'N':virt_numPV,
                'pop_name':'PV',
                'a_name':'PV',
                'rotation_angle_zaxis':xiter_random(N=virt_numPV, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=virt_numPV, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'virtual'
            },
            {   # Interneuron - SOM Cell
                'N':virt_numSOM,
                'pop_name':'SOM',
                'a_name':'SOM',
                'rotation_angle_zaxis':xiter_random(N=virt_numSOM, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=virt_numSOM, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'virtual'
            },
            {   # Interneuron - CR Cell
                'N':virt_numCR,
                'pop_name':'CR',
                'a_name':'CR',
                'rotation_angle_zaxis':xiter_random(N=virt_numCR, min_x=np.pi/2, max_x=np.pi/2),
                'rotation_angle_yaxis':xiter_random(N=virt_numCR, min_x=np.pi/2, max_x=np.pi/2),
                'model_type':'virtual'
            }
        ]
    }
    for cell in shell_network['cells']:
        print(f"Number of shell {cell['pop_name']} cells to be built: {cell['N']}")
        num_shell_cells += cell['N']
    print(f"Total shell cells: {sum([cell['N'] for cell in shell_network['cells']])}")
    # Add the shell to our network definitions
    network_definitions.append(shell_network)

##########################################################################
##########################################################################

init_connectors(num_cells+num_shell_cells)

# Build and save our NetworkBuilder dictionary
networks = build_networks(network_definitions)

# A few connectors require a list for tracking synapses that are recurrent, declare them here 
int2int_temp_list = []
uncoupled_bi_track = []
pyr_int_bi_list = []


# Whole reason for restructuring network building lies here, by separating out the
# source and target params from the remaining parameters in NetworkBuilder's 
# add_edges function we can reuse connectivity rules for the virtual shell
# or elsewhere
# [
#    {
#       'network':'network_name', # => The name of the network that these edges should be added to (networks['network_name'])
#       'edge': {
#                    'source': {},
#                    'target': {}
#               }, # should contain source and target only, any valid add_edges param works
#       'param': 'name_of_edge_parameter' # to be coupled with when add_edges is called
#       'add_properties': 'prop_name' # name of edge_add_properties for adding additional connection props, like delay
#    }
# ]

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
            'source':{'pop_name': ['PV']}, 
            'target':{'pop_name': ['PV']}
        },
        'param': 'PV2PV',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Bidirectional Pair
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PV']}, 
            'target':{'pop_name': ['PV']}
        },
        'param': 'PV2PV_bi_1',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PV Uncoupled Bidirectional Pair
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PV']}, 
            'target':{'pop_name': ['PV']}
        },
        'param': 'PV2PV_bi_2',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PYR Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PV']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        },
        'param': 'PV2PYR',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PYR to PV Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['PV']}
        },
        'param': 'PYR2PV',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PV to PYR Bidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PV']}, 
            'target':{'pop_name': ['PyrA','PyrC']}
        },
        'param': 'PV2PYR_bi',
        'add_properties': 'syn_dist_delay_feng_section_default'
    },
    {   # PYR to PV Bidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PyrA','PyrC']}, 
            'target':{'pop_name': ['PV']}
        },
        'param': 'PYR2PV_bi',
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
    {   # PV to SOM Unidirectional 
        'network':'BLA',
        'edge': {
            'source':{'pop_name': ['PV']}, 
            'target':{'pop_name': ['SOM']}
        },
        'param': 'PV2SOM',
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
            'target':{'pop_name': ['PV']}
        },
        'param': 'CR2PV',
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
        ################## GAP JUNCTIONS #####################

    {   # SOM to SOM GAP JUNCTION
        'network':'BLA',
        'edge': {
            'source':networks['BLA'].nodes(pop_name=['SOM']), 
            'target':networks['BLA'].nodes(pop_name=['SOM'])
        },
        'param': 'SOM2SOM_GAP',
    },

        ##################### VPSI INPUT #####################

    {   # VPSI Inhibition to Pyramidal
        'network':'BLA',
        'edge': {
            'source':networks['vpsi_inh'].nodes(),
            'target':networks['BLA'].nodes(pop_name=['PyrA','PyrC'])
        },
        'param': 'VPSIinh2PYR',
        'add_properties': 'syn_uniform_delay_section_default'
    },
    {   # VPSI Inhibition to PV
        'network':'BLA',
        'edge': {
            'source':networks['vpsi_inh'].nodes(),
            'target':networks['BLA'].nodes(pop_name=['PV'])
        },
        'param': 'VPSIinh2PV',
        'add_properties': 'syn_uniform_delay_section_default'
    },

        ##################### THALAMIC INPUT #####################

    {   # Thalamus to Pyramidal
        'network':'BLA',
        'edge': {
            'source':networks['thalamus_pyr'].nodes(),
            'target':networks['BLA'].nodes(pop_name=['PyrA','PyrC'])
        },
        'param': 'THALAMUS2PYR',
        'add_properties': 'syn_uniform_delay_section_default'        
    },
    {   # Thalamus to Pyramidal
        'network':'BLA',
        'edge': {
            'source':networks['thalamus_pv'].nodes(),
            'target':networks['BLA'].nodes(pop_name=['PV'])
        },
        'param': 'THALAMUS2PV',
        'add_properties': 'syn_uniform_delay_section_default'
    },
    {   # Thalamus to SOM
        'network':'BLA',
        'edge': {
            'source':networks['thalamus_som'].nodes(),
            'target':networks['BLA'].nodes(pop_name='SOM')
        },
        'param': 'THALAMUS2SOM',
        'add_properties': 'syn_uniform_delay_section_default'        
    },
    {   # Thalamus  to CR
        'network':'BLA',
        'edge': {
            'source':networks['thalamus_cr'].nodes(),
            'target':networks['BLA'].nodes(pop_name='CR')
        },
        'param': 'THALAMUS2CR',
        'add_properties': 'syn_uniform_delay_section_default'        
    }
]

# edge_params should contain additional parameters to be added to add_edges calls
edge_params = {
    'PYR2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.02/scale, 'angle_dist':False, 'min_dist':0, 'max_dist':max_conn_dist, 'angle_dist_radius': 200},
        'syn_weight':1,
        'dynamics_params':'PN2PN_feng_min.json',
        'distance_range':[0,max_conn_dist],
        'target_sections':['basal']
    },
    'PV2PV': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.17/scale,'no_recip':True,'track_list':int2int_temp_list, 'max_dist':max_conn_dist},#0.19
        'syn_weight':1,
        'dynamics_params':'INT2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PV2PV_bi_1': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.0275/scale, 'track_list':uncoupled_bi_track, 'max_dist':max_conn_dist},#0.03
        'syn_weight':1,
        'dynamics_params':'INT2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PV2PV_bi_2': {
        'iterator':'one_to_all',
        'connection_rule':recurrent_connector_o2a,
        'connection_params':{'p':1, 'all_edges':uncoupled_bi_track},#p:1
        'syn_weight':1,
        'dynamics_params':'INT2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PV2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.40/scale, 'max_dist':max_conn_dist},#{'p':0.40},
        'syn_weight':1,
        'dynamics_params':'INT2PN_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PYR2PV': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.22/scale, 'angle_dist':False, 'max_dist':max_conn_dist, 'angle_dist_radius': 100},#'p':0.24
        'syn_weight':1,
        'dynamics_params':'PN2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'PV2PYR_bi': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.09/scale,'track_list':pyr_int_bi_list, 'max_dist':max_conn_dist},
        'syn_weight':1,
        'dynamics_params':'INT2PN_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PYR2PV_bi': {
        'iterator':'one_to_all',
        'connection_rule':recurrent_connector_o2a,
        'connection_params':{'p':1,'all_edges':pyr_int_bi_list},#was 1
        'syn_weight':1,
        'dynamics_params':'PN2INT_feng_min.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'PYR2SOM': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.31/scale, 'angle_dist':False, 'max_dist':max_conn_dist, 'angle_dist_radius': 100},#0.309
        'syn_weight':1,
        'dynamics_params':'PN2SOM_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'SOM2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.066/scale, 'max_dist':max_conn_dist},#0.066
        'syn_weight':1,
        'dynamics_params':'SOM2PN_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PV2SOM': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.55/scale, 'max_dist':max_conn_dist},# Dr Unal suggested .1 -> .55 based on 7/1/21 email
        'syn_weight':1,
        'dynamics_params':'INT2SOM_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'PYR2CR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.185/scale, 'angle_dist':False, 'max_dist':max_conn_dist, 'angle_dist_radius': 100},#0.183
        'syn_weight':1,
        'dynamics_params':'PN2CR_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'CR2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.116/scale, 'max_dist':max_conn_dist},#0.116
        'syn_weight':1,
        'dynamics_params':'CR2PN_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['somatic']
    },
    'CR2PV': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.297/scale, 'max_dist':max_conn_dist},#.297
        'syn_weight':1,
        'dynamics_params':'CR2INT_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'CR2SOM': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.764/scale, 'max_dist':max_conn_dist},#.764
        'syn_weight':1,
        'dynamics_params':'CR2SOM_tyler.json',
        'distance_range':[min_conn_dist,max_conn_dist],
        'target_sections':['basal']
    },
    'VPSIinh2PYR': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.10/scale}, # connect a PN cell to a VPSI cell 10% of the time
        'syn_weight':1,
        'dynamics_params':'VPSI2PN_inh_tyler_min.json',
        'distance_range':[0.0, 9999.9],
        'target_sections':['basal'],
    },
    'VPSIinh2PV': {
        'iterator':'one_to_all',
        'connection_rule':syn_percent_o2a,
        'connection_params':{'p':0.90/scale}, # connect a PV cell to a VPSI cell 90% of the time #{'p':0.012/scale}, 
        'syn_weight':1,
        'dynamics_params':'VPSI2PV_inh_tyler_min.json',
        'distance_range':[0.0, 9999.9],
        'target_sections':['basal']
    },
    'THALAMUS2PYR': {
        'connection_rule':one_to_one,
        'syn_weight':1,
        'dynamics_params':'BG2PNe_thalamus_min.json',
        'distance_range':[0.0, 9999.9],
        'target_sections':['basal']
    },
    'THALAMUS2PV': {
        'connection_rule':one_to_one_offset,
        'connection_params':{'offset':numPN_A+numPN_C},
        'syn_weight':1,
        'dynamics_params':'BG2PNi_feng_min.json',
        'distance_range':[0.0, 9999.9],
        'target_sections':['basal']
    },
    'THALAMUS2SOM': {
        'connection_rule':one_to_one_offset,
        'connection_params':{'offset':numPN_A+numPN_C+numPV},
        'syn_weight':1,
        'target_sections':['basal'],
        'distance_range':[0.0, 9999.9],
        'dynamics_params':'BG2SOM_thalamus_min.json'
    },
    'THALAMUS2CR': {
        'connection_rule':one_to_one_offset,
        'connection_params':{'offset':numPN_A+numPN_C+numPV+numSOM},
        'syn_weight':1,
        'target_sections':['basal'],
        'distance_range':[0.0, 9999.9],
        'dynamics_params':'BG2CR_thalamus_min.json'
    },
    ### GAP JUNCTIONS ###
    'SOM2SOM_GAP': {
        'connection_rule':simple_rand,
        'connection_params':{'p':0.08},
        'resistance':0.0005, # Conductance
        'is_gap_junction':True,
    },
} # edges referenced by name

# Will be called by conn.add_properties for the associated connection
edge_add_properties = {
    'syn_dist_delay_feng_section_default': {
        'names':['delay','sec_id','sec_x'],
        'rule':syn_dist_delay_feng_section,
        'rule_params':{'sec_x':0.9},
        'dtypes':[float, np.int32, float]
    },
    'syn_uniform_delay_section_default': {
        'names':['delay','sec_id','sec_x'],
        'rule':syn_uniform_delay_section,
        'rule_params':{'sec_x':0.9},
        'dtypes':[float, np.int32, float]
    },
}

##########################################################################
############################  EDGE EFFECTS  ##############################

if edge_effects:
    # These rules are for edge effect edges. They should directly mimic the connections
    # created previously, re-use the params set above. This keeps our code DRY
    virt_edges = [
        {   # Pyramidal to Pyramidal Connections
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['PyrA','PyrC']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'PYR2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # PV to PV Uncoupled Unidirectional
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['PV']}), 
                'target':{'pop_name': ['PV']}
            },
            'param': 'PV2PV',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
            # PV to PV Uncoupled Bidirectional Pair N/A
            # PV to PV Uncoupled Bidirectional Pair N/A
        {   # PV to PYR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['PV']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'PV2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # PYR to PV Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['PyrA','PyrC']}), 
                'target':{'pop_name': ['PV']}
            },
            'param': 'PYR2PV',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
            # PV to PYR Bidirectional N/A
            # PYR to PV Bidirectional N/A
        {   # PYR to SOM Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['PyrA','PyrC']}), 
                'target':{'pop_name': ['SOM']}
            },
            'param': 'PYR2SOM',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # SOM to PYR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['SOM']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'SOM2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # PV to SOM Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['PV']}), 
                'target':{'pop_name': ['SOM']}
            },
            'param': 'PV2SOM',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # PYR to CR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['PyrA','PyrC']}), 
                'target':{'pop_name': ['CR']}
            },
            'param': 'PYR2CR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # CR to PYR Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['CR']}), 
                'target':{'pop_name': ['PyrA','PyrC']}
            },
            'param': 'CR2PYR',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # CR to PV Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['CR']}), 
                'target':{'pop_name': ['PV']}
            },
            'param': 'CR2PV',
            'add_properties': 'syn_dist_delay_feng_section_default'
        },
        {   # CR to SOM Unidirectional 
            'network':'BLA',
            'edge': {
                'source':networks['shell'].nodes(**{'pop_name': ['CR']}), 
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
# see synapses.py - loads each json's in components/synaptic_models into a 
# dictionary so the properties can be referenced in the files eg: syn['file.json'].get('property')
synapses.load()
syn = synapses.syn_params_dicts()

# Build your edges into the networks
build_edges(networks, edge_definitions,edge_params,edge_add_properties,syn,use_gap_junctions=use_gap_junctions)

# Save the network into the appropriate network dir
save_networks(networks,network_dir)

# Usually not necessary if you've already built your simulation config
build_env_bionet(base_dir='./',
	network_dir=network_dir,
	tstop=t_sim, dt = dt,
	report_vars = ['v'],
        v_init = -70.0,
        celsius = 31.0,
	spikes_inputs=[('vpsi_inh','vpsi_inh_spikes.h5'),# Name of population which spikes will be generated for, file
                       ('thalamus_pyr','thalamus_pyr_spikes.h5'),
                       ('thalamus_som','thalamus_som_spikes.h5'),
                       ('thalamus_cr','thalamus_cr_spikes.h5'),
                       ('shell','shell_spikes.h5')],
	components_dir=components_dir,
        config_file='simulation_config.json',
	compile_mechanisms=False,
        overwrite_config=True)

if edge_effects:
    from build_input_shell import build_shell_inputs

    # This needs to be called after building and saving your networks
    # There's an optimization in there that determines which shells in the
    # shell are connected to bio cells, and only delivers a spike train to
    # those, excluding others and speeding up the simulation.
    # Since edges are semi-random, we want to deliver to the correct cells
    # each time.
    build_shell_inputs()

from build_input import build_input

build_input(t_sim,
            numPN_A = numPN_A,
            numPN_C = numPN_C,
            numPV = numPV,
            numSOM = numSOM,
            numCR = numCR,
            scale = scale)

#Generate ACH input nodes
from generate_node_sets import generate_node_sets
generate_node_sets(scale)

#Build VSPI input spikes
from build_input_vpsi_inh_spikes import build_vpsi_input, build_vpsi_input_jitter
#build_vpsi_input(t_sim=t_sim, n_cells=numVPSI, depth_of_mod=1, output='vpsi_inh_spikes.h5')
build_vpsi_input_jitter(t_sim=t_sim, n_cells=numVPSI, output='vpsi_inh_spikes.h5')
build_vpsi_input(t_sim=t_sim,n_cells=numVPSI, depth_of_mod=0, output='vpsi_inh_spikes_0_depth.h5')

profile_stats = pstats.Stats(profiler).sort_stats('tottime')
profile_stats.print_stats(100)
