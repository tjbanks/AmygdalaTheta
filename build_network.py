from bmtk.builder import NetworkBuilder
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list, xiter_random
import synapses
import math
import pdb
import random
import pandas as pd

np.random.seed(123412)

# Initialize our network
net = NetworkBuilder("BLA")

scale = 1

all_synapses = pd.DataFrame([],columns=['source_gid','target_gid'])

#Number of cells in each population
numPN_A = 640 * scale #4114#15930
numPN_C = 260 * scale #4115#6210
numBask = 100 * scale #854#4860
numSOM = 42 * scale
numCR = 42 * scale
# add_properties = False
# do_pos = False
num_cells = numPN_A + numPN_C + numBask

dist_constraint = True 
min_conn_dist = 0.0
max_conn_dist = 9999.9

if dist_constraint:
    max_conn_dist = 300.0 #9999.9# Distance constraint for all cells

i2i_gap = False
connect_som = False

if connect_som:
    num_cells = num_cells + numSOM
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
              model_template='hoc:Cell_Af',
              morphology=None)

# Get rid of coordinates already used
pos_list = np.delete(pos_list,inds,0)

##################################################################################
###################################Pyr Type C#####################################

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
              model_template='hoc:Cell_Cf',
              morphology=None)

# Get rid of coordinates already used
pos_list = np.delete(pos_list,inds,0)

    
#################################################################################
###########################Fast - spiking PV ints################################

# Pick new coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numBask,replace=False)
pos = pos_list[inds,:]

bask_pos = pos.copy()
nid_pos = np.concatenate([pyra_pos, pyrc_pos, bask_pos])

# Add a population of numBask nodes
net.add_nodes(N=numBask, pop_name='Bask',
              positions=positions_list(positions=pos),
	      rotation_angle_zaxis=xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
	      rotation_angle_yaxis=xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:InterneuronCellf',
              morphology=None)

pos_list = np.delete(pos_list,inds,0)

#################################################################################
################################# SOM Cells #####################################

if connect_som:
    # Pick new coordinates
    inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numSOM,replace=False)
    pos = pos_list[inds,:]

    som_pos = pos.copy()
    nid_pos = np.concatenate([pyra_pos, pyrc_pos, bask_pos,som_pos])

    # Add a population of numBask nodes
    net.add_nodes(N=numSOM, pop_name='SOM',
              positions=positions_list(positions=pos),
              rotation_angle_zaxis=xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
              rotation_angle_yaxis=xiter_random(N=numBask, min_x=0.0, max_x=2*np.pi),
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:SOM_Cell',
              morphology=None)


    pos_list = np.delete(pos_list,inds,0)
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

def one_to_one(source, target):
    
    sid = source.node_id
    tid = target.node_id
    if sid == tid:
    #print("connecting cell {} to {}".format(sid,tid))
        tmp_nsyn = 1
    else:
        return None

    return tmp_nsyn

def one_to_one_offset(source, target, offset=0):

    sid = source.node_id
    tid = target.node_id - offset
    if sid == tid:
        #print("connecting cell {} to {}".format(sid,tid))
        tmp_nsyn = 1
    else:
        #print("NOT connecting cell {} to {}".format(sid,tid))
        return None

    return tmp_nsyn


def syn_dist_delay_feng(source, target):
    #if not dist_constraint:
    #    return 0.1

    dt = 0.05
    min_delay=0.8   #////define minmum delay,ms
    #maxdis=2.425   #/// mm sqrt((1.4)^2+(1.4)^2+(1.4)^2)
    x = float(x_end - x_start)#/1000
    y = float(y_end - y_start)#/1000
    z = float(z_end - z_start)#/1000
    max_dist = 2.425#np.sqrt(x**2 + y**2 + z**2)
    max_delay=2.425 #////define maximum delay,ms

    x_ind,y_ind,z_ind = 0,1,2

    dx = target['positions'][x_ind] - source['positions'][x_ind]
    dy = target['positions'][y_ind] - source['positions'][y_ind]
    dz = target['positions'][z_ind] - source['positions'][z_ind]

    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    
    del_fluc = np.random.uniform(-0.1,0.1)

    delay=(dist/max_dist)*max_delay+min_delay+del_fluc+dt

    return delay

def syn_dist_delay_feng_section(source, target, sec_id=0, sec_x=0.9):
    return syn_dist_delay_feng(source, target), sec_id, sec_x

def syn_uniform_delay_section(source, target, sec_id=0, sec_x=0.9, mean=0.5,std=1):
    return np.random.uniform(mean,std), sec_id, sec_x

def syn_percent(source,target,p,track_list=None):
    """
    track_list: supply a list to append and track the synapses with
    """
    global all_synapses

    sid = source.node_id
    tid = target.node_id
    # No autapses
    if sid==tid:
        return None

    # Do not add synapses if they already exist, we don't want duplicates
    if ((all_synapses.source_gid == sid) & (all_synapses.target_gid == tid)).any():
        return None

    if random.random() < p:
        all_synapses = all_synapses.append({'source_gid':sid,'target_gid':tid},ignore_index=True)
        if track_list is not None:#we only want to track synapses that may have a recurrent connection, will speed up build time considerably
            track_list.append({'source_gid':source['node_id'],'target_gid':target['node_id']})        
        return 1
    else:
        return 0

def syn_percent_o2a(source,targets,p,track_list=None,no_recip=False):
    """
    track_list: supply a list to append and track the synapses with
    one to all connector for increased speed.
    returns a list of len(targets) where the value is the number of synapses at the index=cellid
    if no_recip is set to true, any connection that has been made previously and is reciprical won't be considered
    """

    global all_synapses
    sid = source.node_id 
    tids = np.array([target.node_id for target in targets])

    #List of synapses that will be returned, initialized to 0 synapses
    syns = np.array([0 for _ in range(len(targets))])

    #Get all existing targets for that source that can't be targetted here
    existing = all_synapses[all_synapses.source_gid == sid]
    existing_list = existing[existing.target_gid.isin(tids)].target_gid.tolist()

    #remove existing synapses from available list
    available = tids.copy()
    available = available[~np.isin(available,existing_list)]
    
    if no_recip:    
        recur = [i['source_gid'] for i in track_list if i['target_gid'] == sid]
        available = available[~np.isin(available,recur)]
    
    # of those remaining we want p% chosen
    n = int(len(tids)*p)
    extra = 1 if random.random() < (p*100 - math.floor(p*100)) else 0
    n = n + extra #This will account for half percentages
    chosen = np.random.choice(available,size=n,replace=False) 

    syns[np.isin(tids,chosen)] = 1
    
    #Add to lists
    new_syns = pd.DataFrame(chosen,columns=['target_gid'])
    new_syns['source_gid'] = sid
    all_synapses = all_synapses.append(new_syns,ignore_index=True)

    if track_list is not None:
        #track_list = track_list.append(new_syns,ignore_index=True)
        for target in chosen:
            track_list.append({'source_gid':sid,'target_gid':target})
    #any index selected will be set to 1 and returned
    return syns 

def recurrent_connector(source,target,p,all_edges=[],min_syn=1, max_syn=1):
    """
    General logic:
    1. Given a *potential* source and target
    2. Look through all edges currently made
    3. If any of the current edges contains 
        a. the current source as a previous target of 
        b. the current target as a prevous source
    4. Return number of synapses per this connection, 0 otherwise (no connection)
    """
    global all_synapses

    sid = source.node_id
    tid = target.node_id
    
    # Do not add synapses if they already exist, we don't want duplicates
    if ((all_synapses.source_gid == sid) & (all_synapses.target_gid == tid)).any():
        return None
    
    for e in all_edges: #should probably make this a pandas df to speed up building... and use .any() to search
        if sid == e['target_gid'] and tid == e['source_gid']:
            #print('found recurrent')

            if random.random() < p:
                #print('--------------connecting')
                all_synapses = all_synapses.append({'source_gid':sid,'target_gid':tid},ignore_index=True)
                return random.randint(min_syn,max_syn)
            else:
                return 0
    return 0

def recurrent_connector_o2a(source,targets,p,all_edges=[],min_syn=1,max_syn=1):

    global all_synapses
    sid = source.node_id
    tids = np.array([target.node_id for target in targets])

    #List of synapses that will be returned, initialized to 0 synapses
    syns = np.array([0 for _ in range(len(targets))])

    existing = all_synapses[all_synapses.source_gid == sid]
    existing_list = existing[existing.target_gid.isin(tids)].target_gid.tolist()

    #remove existing synapses from available list
    available = tids.copy()
    available = available[~np.isin(available,existing_list)]

    #remove any connection that is not in the all_edges list from 'available' list
    recur = [i['source_gid'] for i in all_edges if i['target_gid'] == sid]
    available = available[np.isin(available,recur)]
    #import pdb;pdb.set_trace()

    # of those remaining we want p% chosen
    n = int(len(available)*p)
    extra = 1 if random.random() < (p*100 - math.floor(p*100)) else 0
    n = n + extra #This will account for half percentages
    chosen = np.random.choice(available,size=n,replace=False) 

    syns[np.isin(tids,chosen)] = 1
    
    #Add to lists
    new_syns = pd.DataFrame(chosen,columns=['target_gid'])
    new_syns['source_gid'] = sid
    all_synapses = all_synapses.append(new_syns,ignore_index=True)

    #any index selected will be set to 1 and returned
    return syns

# Create connections between Pyr --> Pyr cells
add_delays = []#Says whether the next add_edges should have delays added by distance.
min_delays = []#Stores min_delay for each synapse type to be used later.


##########################################################################
############################### PYR2PYR ##################################

p2p_props = [
    {'min_dist': 0, 'max_dist': 50, 'syn_prob': 0.02},    #0.03
    {'min_dist': 51, 'max_dist': 100, 'syn_prob': 0.02},  #0.02
    {'min_dist': 101, 'max_dist': 200, 'syn_prob': 0.02}, #0.01
    {'min_dist': 201, 'max_dist': 300, 'syn_prob': 0.02}, #0.005
    {'min_dist': 301, 'max_dist': 400, 'syn_prob': 0.02}, #0.005
    {'min_dist': 401, 'max_dist': 500, 'syn_prob': 0.02}, #0.005
    {'min_dist': 501, 'max_dist': 600, 'syn_prob': 0.02}, #0.005
]

p2p_props = [
    {'min_dist':0, 'max_dist':max_conn_dist, 'syn_prob': 0.02}
]

# What we're doing here is looping through the different connection
# probabilities based on distance apart instead of re-writing this
# large block of code several times

for p2p_prop in p2p_props:
    #dynamics_file = 'PN2PN.json'
    dynamics_file = 'PN2PN_feng_min.json'

    conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': ['PyrA','PyrC']},
                iterator = 'one_to_all',
                connection_rule=syn_percent_o2a,
                connection_params={'p':p2p_prop['syn_prob']},
                syn_weight=1,
                dynamics_params=dynamics_file,
                model_template=syn[dynamics_file]['level_of_detail'],
                distance_range=[p2p_prop['min_dist'],p2p_prop['max_dist']],
                target_sections=['basal'])

    conn.add_properties(names=['delay','sec_id','sec_x'],
                rule=syn_dist_delay_feng_section,
                rule_params={'sec_id':1, 'sec_x':0.9},
                dtypes=[np.float, np.int32, np.float])

##########################################################################
############################### GAPJUNC ##################################

gap_list = []

if i2i_gap:
    conn = net.add_gap_junctions(source={'pop_name': ['Bask']}, 
		      target={'pop_name': ['Bask']},
 		      resistance = 0.0001, target_sections=['somatic'], 
		      connection_rule=syn_percent,
		      connection_params={'p':0.08, 'track_list':gap_list})
    conn._edge_type_properties['sec_id'] = 0
    conn._edge_type_properties['sec_x'] = 0.9


##########################################################################
############################### INT2INT ##################################
####################### Coupled Unidirectional Pair ######################

##########################################################################
############################### INT2INT ##################################
####################### Coupled Bidirectional Pair #######################


##########################################################################
############################### INT2INT ##################################
##################### Uncoupled Unidirectional Pair ######################

# Create connections between Bask --> Bask cells
dynamics_file = 'INT2INT.json'
dynamics_file = 'INT2INT_feng_min.json'

temp_list = []
conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['Bask']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.19,'no_recip':True,'track_list':temp_list},#0.19
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
############################### INT2INT ##################################
###################### Uncoupled Bidirectional Pair ######################

uncoupled_bi_track = []

dynamics_file = 'INT2INT_feng_min.json'

conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['Bask']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.025, 'track_list':uncoupled_bi_track},#0.03
              syn_weight=1,
              dynamics_params=dynamics_file,
              model_template=syn[dynamics_file]['level_of_detail'],
              distance_range=[min_conn_dist,max_conn_dist],
              target_sections=['somatic'])

conn.add_properties(names=['delay','sec_id','sec_x'],
              rule=syn_dist_delay_feng_section,
              rule_params={'sec_id':0, 'sec_x':0.9},
              dtypes=[np.float, np.int32, np.float])


conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['Bask']},
              iterator = 'one_to_all',
              connection_rule=recurrent_connector_o2a,
              connection_params={'p':1, 'all_edges':uncoupled_bi_track},#p:1
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
############################### INT2PYR ##################################
########################### UNIDIRECTIONAL ###############################

# Create connections between Bask --> Pyr cells
#dynamics_file = 'INT2PN.json'
dynamics_file = 'INT2PN_feng_min.json'


conn = net.add_edges(source={'pop_name': 'Bask'}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.36},#{'p':0.40},
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
############################### PYR2INT ##################################
########################### UNIDIRECTIONAL ###############################

# Create connections between Pyr --> Bask cells
#dynamics_file = 'PN2INT.json'
dynamics_file = 'PN2INT_feng_min.json'

conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'Bask'},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.24},#{'p':0.12},
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


if connect_som:
    dynamics_file = 'PN2SOM_tyler.json'
    conn = net.add_edges(source={'pop_name': ['PyrA','PyrC']}, target={'pop_name': 'SOM'},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.309},
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

    dynamics_file = 'SOM2PN_tyler.json'
    conn = net.add_edges(source={'pop_name': 'SOM'}, target={'pop_name': ['PyrA','PyrC']},
              iterator = 'one_to_all',
              connection_rule=syn_percent_o2a,
              connection_params={'p':0.066},
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
######################### BACKGROUND INPUT ###############################

dynamics_file='BG2PNe_feng_min.json'
#dynamics_file = 'BG2PN_feng_min.json'

conn = net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrA'),
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

conn = net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='PyrC'),
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


dynamics_file = 'BG2PNi_feng_min.json'

conn = net.add_edges(source=exc_bg_bask.nodes(), target=net.nodes(pop_name='Bask'),
                   connection_rule=one_to_one_offset,
                   connection_params={'offset':numPN_A+numPN_C},
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
t_sim = 10000.0

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
		network_dir='./network',
		tstop=t_sim, dt = 0.05,
		report_vars = ['v'],
                v_init = -70.0,
                celsius = 31.0,
		spikes_inputs=[('mthalamus',   # Name of population which spikes will be generated for
                                'mthalamus_spikes.h5'),('exc_bg_bask','exc_bg_bask_spikes.h5')],
		components_dir='components',
		compile_mechanisms=True)


