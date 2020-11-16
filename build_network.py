from bmtk.builder import NetworkBuilder
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list, xiter_random
import synapses
import pdb
import random

np.random.seed(123412)

# Initialize our network
net = NetworkBuilder("BLA")

#Number of cells in each population
numPN_A = 640#4114#15930
numPN_C = 260#4115#6210
numBask = 100#854#4860
# add_properties = False
# do_pos = False
num_cells = numPN_A + numPN_C + numBask

# Create the possible x,y,z coordinates

x_start, x_end = 0,600
y_start, y_end = 0,600
z_start, z_end = 0,600

pos_list = np.random.rand(num_cells,3)
pos_list[:,0] = pos_list[:,0]*x_end - x_start
pos_list[:,1] = pos_list[:,1]*y_end - y_start
pos_list[:,2] = pos_list[:,2]*z_end - z_start

# Load synapse dictionaries
synapses.load()
syn = synapses.syn_params_dicts()
###################################################################################
####################################Pyr Type A#####################################

# Pick coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numPN_A,replace=False)
pos = pos_list[inds,:]

pyra_pos = pos.copy()

# Add a population of numPN_A nodes (all of which share model_type, dynamics_params, etc.)
net.add_nodes(N=numPN_A, pop_name='PyrA',
              positions=positions_list(positions=pos),
	      rotation_angle_zaxis=xiter_random(N=numPN_A, min_x=0.0, max_x=2*np.pi),
	      rotation_angle_yaxis=xiter_random(N=numPN_A, min_x=0.0, max_x=2*np.pi),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:feng_typeA',
              morphology=None)

##################################################################################
###################################Pyr Type C#####################################

# Get rid of coordinates already used
pos_list = np.delete(pos_list,inds,0)

# Pick new coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numPN_C,replace=False)
pos = pos_list[inds,:]

pyrc_pos = pos.copy()

# Add a population of numPN_A nodes (all of which share model_type, dynamics_params, etc.)
net.add_nodes(N=numPN_C, pop_name='PyrC',
              positions=positions_list(positions=pos),
	      rotation_angle_zaxis=xiter_random(N=numPN_C, min_x=0.0, max_x=2*np.pi),
	      rotation_angle_yaxis=xiter_random(N=numPN_C, min_x=0.0, max_x=2*np.pi),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:feng_typeC',
              morphology=None)
#################################################################################
############################# Chandelier ########################################
if False:
    # Get rid of coordinates already used
    pos_list = np.delete(pos_list,inds,0)

    # Pick new coordinates
    inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numAAC,replace=False)
    pos = pos_list[inds,:]

    aac_pos = pos.copy()

    # Add a population of numAAC nodes
    net.add_nodes(N=numAAC, pop_name='AAC',
                positions=positions_list(positions=pos),
            rotation_angle_zaxis=xiter_random(N=numAAC, min_x=0.0, max_x=2*np.pi),
            rotation_angle_yaxis=xiter_random(N=numAAC, min_x=0.0, max_x=2*np.pi),
                mem_potential='e',
                model_type='biophysical',
                model_template='hoc:chandelier',
                morphology=None)
    # Get rid of coordinates already used
    pos_list = np.delete(pos_list,inds,0)
    
#################################################################################
###########################Fast - spiking PV ints################################

# Pick new coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numBask,replace=False)
pos = pos_list[inds,:]

bask_pos = pos.copy()
#nid_pos = np.concatenate([pyra_pos, pyrc_pos, aac_pos, bask_pos])
nid_pos = np.concatenate([pyra_pos, pyrc_pos, bask_pos])
#import pdb; pdb.set_trace()

# Add a population of numBask nodes
net.add_nodes(N=numBask, pop_name='Bask',
              positions=positions_list(positions=pos),
	      rotation_angle_zaxis=xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
	      rotation_angle_yaxis=xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:basket',
              morphology=None)
################################################################################
############################# BACKGROUND INPUTS ################################

# External inputs
thalamus = NetworkBuilder('mthalamus')
thalamus.add_nodes(N=numPN_A+numPN_C,
                   pop_name='tON',
                   pop_group='mthalamus',
                   potential='exc',
                   model_type='virtual')

# External inputs
exc_bg_bask = NetworkBuilder('exc_bg_bask')
exc_bg_bask.add_nodes(N=numBask,
                   pop_name='tON',
                   pop_group='exc_bg_bask',
                   potential='exc',
                   model_type='virtual')
##############################################################################
############################## CONNECT CELLS #################################

def dist_conn_perc_angle(src, trg, min_dist=0.0, max_dist=300.0, min_syns=1, max_syns=2, A=0.2, B=0.2):
    
    sid = src.node_id
    tid = trg.node_id
    # No autapses
    if sid==tid:
        return None
    else:
        src_pos = src['positions']
        trg_pos = trg['positions']
    #dist =np.sqrt((src_pos[0]-trg_pos[0])**2+(src_pos[1]-trg_pos[1])**2+(src_pos[2]-trg_pos[2])**2)
        src_angle_x = src['rotation_angle_zaxis'] 
        src_angle_y = src['rotation_angle_yaxis']

        # We must calculate vec_pos, the coordinates of the point on the tip
        # of the unit vector pointing from src_pos
        vec_pos = [np.cos(src_angle_x), np.sin(src_angle_y), np.sin(src_angle_x)] 
        perp_dist = np.linalg.norm(np.cross((trg_pos-src_pos),(trg_pos-vec_pos)))/np.linalg.norm((vec_pos - src_pos))

        #print("src_pos: {} trg_pos: {} dist: {}".format(src_pos,trg_pos,dist))        
    prob = A

    if perp_dist <= max_dist and np.random.uniform() < prob:
        tmp_nsyn = np.random.randint(min_syns, max_syns)
        #print("creating {} synapse(s) between cell {} and {}".format(tmp_nsyn,sid,tid))
    else:
        tmp_nsyn = 0

    return tmp_nsyn

def dist_conn_perc(src, trg, min_dist=0.0, max_dist=300.0, min_syns=1, max_syns=2, A=0.2, B=0.2):
    
    sid = src.node_id
    tid = trg.node_id
    # No autapses
    if sid==tid:
        return None
    else:
        src_pos = src['positions']
        trg_pos = trg['positions']
    dist =np.sqrt((src_pos[0]-trg_pos[0])**2+(src_pos[1]-trg_pos[1])**2+(src_pos[2]-trg_pos[2])**2)
        #print("src_pos: {} trg_pos: {} dist: {}".format(src_pos,trg_pos,dist))        
    prob = A*np.exp(-B*dist)

    if dist <= max_dist and np.random.uniform() < prob:
        tmp_nsyn = np.random.randint(min_syns, max_syns)
        #print("creating {} synapse(s) between cell {} and {}".format(tmp_nsyn,sid,tid))
    else:
        tmp_nsyn = 0

    return tmp_nsyn

def one_to_one(source, target):
    
    sid = source.node_id
    tid = target.node_id
    if sid == tid:
    #print("connecting cell {} to {}".format(sid,tid))
        tmp_nsyn = 1
    else:
        return None

    return tmp_nsyn


def syn_dist_delay(source, target, min_delay):#, min_weight, max_weight):

	x_ind,y_ind,z_ind = 0,1,2

	dx = target['positions'][x_ind] - source['positions'][x_ind]
	dy = target['positions'][y_ind] - source['positions'][y_ind]
	dz = target['positions'][z_ind] - source['positions'][z_ind]

	dist = np.sqrt(dx**2 + dy**2 + dz**2)
	distDelay = dist/1000
	#print("delay = {}".format(distDelay)) 
	return float(min_delay) + distDelay

def syn_dist_delay_section(source, target, min_delay, sec_id=0, sec_x=0.9):
    return syn_dist_delay(source, target, min_delay), sec_id, sec_x

def syn_delay(source, target, min_delay):
    return [syn_dist_delay(source, target, min_delay)]
    #return [0.1]

def syn_percent(source,target,p):
    return 1 if random.random() < p else 0

# Create connections between Pyr --> Pyr cells
add_delays = []#Says whether the next add_edges should have delays added by distance.
min_delays = []#Stores min_delay for each synapse type to be used later.

dynamics_file = 'PN2PN.json'

add_delays.append(True)
min_delays.append(syn[dynamics_file]['delay'])

conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_one',
              #connection_rule=dist_conn_perc,
              #connection_params={'min_dist':0.0,'max_dist':50.0,
			  #       'min_syns':1,'max_syns':2,'A':0.01366,'B':0.008618},
              connection_rule=syn_percent,
              connection_params={'p':0.03},
              syn_weight=1,
	      delay=0.1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              #distance_range=[0.0, 300.0],
              distance_range=[0.0, 9999.9],
              target_sections=['basal'],
              sec_id=0,
              sec_x=0.5)

# if add_properties:
#     if do_pos:
#         conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                         rule=syn_dist_delay_section,
#                         rule_params={'min_delay':syn[dynamics_file]['delay'], 
#                         'sec_id':0, 'sec_x':0.9},
#                         dtypes=[np.float, np.int32, np.float])
#     else:
#         conn.add_properties(names=['delay'], rule=syn_delay, 
#             rule_params={'min_delay':syn[dynamics_file]['delay']}, dtypes=[np.float])

# Create connections between Pyr --> Bask cells
dynamics_file = 'PN2INT.json'

add_delays.append(True)
min_delays.append(syn[dynamics_file]['delay'])

conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'Bask'},
              iterator = 'one_to_one',
	          #connection_rule=dist_conn_perc,
              #connection_params={'min_dist':0.0,'max_dist':50.0,
			  #       'min_syns':1,'max_syns':2,'A':0.3217,'B':0.005002},
              connection_rule=syn_percent,
              connection_params={'p':0.12},
              syn_weight=1,
	          delay = 0.1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              #distance_range=[0.0, 300.0],
              distance_range=[0.0, 9999.9],
              target_sections=['somatic'],
              sec_id=0,
              sec_x=0.5)

# if add_properties:
#     if do_pos:
#         conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                         rule=syn_dist_delay_section,
#                         rule_params={'min_delay':syn[dynamics_file]['delay'], 
#                         'sec_id':0, 'sec_x':0.9},
#                         dtypes=[np.float, np.int32, np.float])
#     else:
#         conn.add_properties(names=['delay'], rule=syn_delay, 
#             rule_params={'min_delay':syn[dynamics_file]['delay']}, dtypes=[np.float])



# Create connections between Pyr --> AAC cells
if False:
    dynamics_file = 'PN2INT.json'

    add_delays.append(True)
    min_delays.append(syn[dynamics_file]['delay'])

    conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'AAC'},
                iterator = 'one_to_one',
            connection_rule=dist_conn_perc,
                connection_params={'min_dist':0.0,'max_dist':50.0,
                        'min_syns':1,'max_syns':2,'A':0.3217,'B':0.005002},
                syn_weight=1,
            delay = 0.1,
                dynamics_params=dynamics_file,
                model_template=syn[dynamics_file]['level_of_detail'],
                distance_range=[0.0, 300.0],
                target_sections=['somatic'],
                sec_id=0,
                sec_x=0.5)

# Create connections between Bask --> Pyr cells
dynamics_file = 'INT2PN.json'

add_delays.append(True)
min_delays.append(syn[dynamics_file]['delay'])

conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_one',
              #connection_rule=dist_conn_perc,
              #connection_params={'min_dist':0.0,'max_dist':100000.0,
			  #   'min_syns':1,'max_syns':2,'A':0.313,'B':0.004029},
              connection_rule=syn_percent,
              connection_params={'p':0.34},
              syn_weight=1,
              delay=0.1,
              dynamics_params='INT2PN.json',
              model_template=syn['INT2PN.json']['level_of_detail'],
              #distance_range=[0.0, 300.0],
              distance_range=[0.0, 9999.9],
              target_sections=['somatic'],
              sec_id=0,
              sec_x=0.5)

# if add_properties:
#     if do_pos:
#         conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                         rule=syn_dist_delay_section,
#                         rule_params={'min_delay':syn[dynamics_file]['delay'], 
#                         'sec_id':0, 'sec_x':0.9},
#                         dtypes=[np.float, np.int32, np.float])
#     else:
#         conn.add_properties(names=['delay'], rule=syn_delay, 
#             rule_params={'min_delay':syn[dynamics_file]['delay']}, dtypes=[np.float])




# Create connections between AAC --> Pyr cells
if False:
    dynamics_file = 'INT2PN.json'

    add_delays.append(True)
    min_delays.append(syn[dynamics_file]['delay'])

    conn = net.add_edges(source={'pop_name': 'AAC'}, target={'pop_name': ['PyrA','PyrC']},
                iterator = 'one_to_one',
                connection_rule=dist_conn_perc,
                connection_params={'min_dist':0.0,'max_dist':100000.0,
                    'min_syns':1,'max_syns':2,'A':0.313,'B':0.004029},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[0.0, 300.0],
                target_sections=['somatic'],
                sec_id=0,
                sec_x=0.5)

# Create connections between Bask --> Bask cells
dynamics_file = 'INT2INT.json'

add_delays.append(True)
min_delays.append(syn[dynamics_file]['delay'])

conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['Bask']},
              iterator = 'one_to_one',
              #connection_rule=dist_conn_perc,
              #connection_params={'min_dist':0.0,'max_dist':300.0,
			  #   'min_syns':1,'max_syns':2,'A':0.24,'B':0.0},
              connection_rule=syn_percent,
              connection_params={'p':0.26},
              syn_weight=5.0e-03,
              weight_function='lognormal',
              weight_sigma=1.0e-03,
	          delay=0.1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              #distance_range=[0.0, 300.0],
              distance_range=[0.0, 9999.9],
              target_sections=['somatic'],
              sec_id=0,
              sec_x=0.5)

# if add_properties:
#     if do_pos:
#         conn.add_properties(names=['delay', 'sec_id', 'sec_x'],
#                         rule=syn_dist_delay_section,
#                         rule_params={'min_delay':syn[dynamics_file]['delay'], 
#                         'sec_id':0, 'sec_x':0.9},
#                         dtypes=[np.float, np.int32, np.float])
#     else:
#         conn.add_properties(names=['delay'], rule=syn_delay, 
#             rule_params={'min_delay':syn[dynamics_file]['delay']}, dtypes=[np.float])

add_delays.append(False)
min_delays.append(-1)#Want to append -1 if not adding delays.

#conn = net.add_gap_junctions(source={'pop_name': ['Bask']}, 
#		      target={'pop_name': ['Bask']},
# 		      resistance = 0.0001, target_sections=['somatic'], 
#		      connection_rule=dist_conn_perc,
#		      connection_params={'min_dist':0.0,
#					'max_dist':300.0,'min_syns':1,
#					'max_syns':2,'A':0.08,'B':0.0})
#conn._edge_type_properties['sec_id'] = 0
#conn._edge_type_properties['sec_x'] = 0.9

"""

net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrA'),
                   connection_rule=one_to_one,
                   syn_weight=7.0e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='exp2syn')

net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrC'),
                   connection_rule=one_to_one,
                   syn_weight=7.0e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='exp2syn')

net.add_edges(source=exc_bg_bask.nodes(), target=net.nodes(pop_name='Bask'),
                   connection_rule=one_to_one,
                   syn_weight=7.0e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='exp2syn')

"""

net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrA'),
                   connection_rule=one_to_one,
                   syn_weight=8.5e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 9999.9],
                   dynamics_params='BG2PN_feng.json',
                   model_template='bg2pyr',
                   sec_x=0.9)

net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrC'),
                   connection_rule=one_to_one,
                   syn_weight=8.5e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-03,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 9999.9],
                   dynamics_params='BG2PN_feng.json',
                   model_template='bg2pyr',
                   sec_x=0.9)

net.add_edges(source=exc_bg_bask.nodes(), target=net.nodes(pop_name='Bask'),
                   connection_rule=one_to_one,
                   syn_weight=1.0e-03,
                   weight_function='lognormal',
                   weight_sigma=2.0e-04,
                   target_sections=['basal'],
                   delay=0.1,
                   distance_range=[0.0, 9999.9],
                   dynamics_params='BG2PN_feng.json',
                   model_template='bg2pyr',
                   sec_x=0.9)

net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

print("Internal nodes and edges built")

# Create connections between "thalamus" and Pyramidals
# First define the connection rule

# Build and save our network

thalamus.build()
thalamus.save_nodes(output_dir='network')

exc_bg_bask.build()
exc_bg_bask.save_nodes(output_dir='network')
#
#print("External nodes and edges built")
t_sim = 5000

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
		network_dir='./network',
		tstop=t_sim, dt = 0.1,
		report_vars = ['v'],
		spikes_inputs=[('mthalamus',   # Name of population which spikes will be generated for
                                'mthalamus_spikes.h5'),('exc_bg_bask','exc_bg_bask_spikes.h5')],
		#current_clamp={     
                #     'gids': [0],
                #     'amp': [0.5], 
                #     'delay': 100.0, 
                #     'duration': 50.0 
                # },
		components_dir='components',
		compile_mechanisms=True)


from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

psg = PoissonSpikeGenerator(population='mthalamus')
psg.add(node_ids=range(numPN_A+numPN_C),  # Have nodes to match mthalamus
        firing_rate=0.002,    # 15 Hz, we can also pass in a nonhomoegenous function/array
        times=(0.0, t_sim))    # Firing starts at 0 s up to 3 s
psg.to_sonata('mthalamus_spikes.h5')

psg = PoissonSpikeGenerator(population='exc_bg_bask')
psg.add(node_ids=range(numBask),  # Have nodes to match mthalamus
        firing_rate=0.002,    # 15 Hz, we can also pass in a nonhomoegenous function/array
        times=(0.0, t_sim))    # Firing starts at 0 s up to 3 s
psg.to_sonata('exc_bg_bask_spikes.h5')


exit()

def syn_dist_delay(source, target, min_delay, pos):
    x_ind,y_ind,z_ind = 0,1,2
    dx = pos[target][x_ind] - pos[source][x_ind]
    dy = pos[target][y_ind] - pos[source][y_ind]
    dz = pos[target][z_ind] - pos[source][z_ind]

    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    distDelay = dist/500
    #print("delay = {}".format(distDelay)) 
    return float(min_delay) + distDelay


import h5py
import numpy as np
f = h5py.File('network/BLA_BLA_edges.h5', 'r+')
#import pdb; pdb.set_trace()
edges = f['edges']['BLA_to_BLA']
del edges['0']
#f.move(edges.name+'/0', edges.name+'/1')
edges.create_group('0')
edges.create_group('1')
del_group = edges['0']

types = edges['edge_type_id'].value - 100
sources = edges['source_node_id'].value
targets = edges['target_node_id'].value
group_id = []
group_index = []
delays = []
num_no_d = 0
num_d = 0 

for i in range(len(sources)):
    if add_delays[types[i]]:
        delays.append(syn_dist_delay(sources[i], targets[i], min_delays[types[i]], nid_pos))
        group_id.append(0)
        group_index.append(num_d)
        num_d += 1
    else:
        group_id.append(1)
        group_index.append(num_no_d)
        num_no_d += 1

#delays = np.array([syn_dist_delay(sources[i], targets[i], types[i], dyns, nid_pos) for i in range(len(sources)) if types[i] != 104])
# np.save("delays_test", np.array(delays))
# np.save("gids_test", np.array(group_id))
# np.save("index_test", np.array(group_index))
edges['1'].create_dataset('nsyns', data=np.full(num_no_d, 1), dtype='uint16')

del_group.create_dataset('delay', data=np.array(delays))
#del_group.create_dataset('sec_x', data=np.full(len(delays), 0.9))
#del_group.create_dataset('sec_id', data=np.full(len(delays), 0))
edges['edge_group_id'][...] = np.array(group_id)
edges['edge_group_index'][...] = np.array(group_index)
f.close()
