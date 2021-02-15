from bmtk.builder import NetworkBuilder
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list, xiter_random
import pdb
import random
import pandas as pd

#A bunch of junk for importing the synapses loader
import sys, os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))
sys.path.insert(0, parentdir)
import synapses
import build_input

np.random.seed(123412)

# Initialize our network
net = NetworkBuilder("BLA")

scale = 1

all_synapses = pd.DataFrame([],columns=['source_gid','target_gid'])

#Number of cells in each population
numPN_A = 1
numPN_C = 1
numBask = 1
num_cells = numPN_A + numPN_C + numBask

dist_constraint = False
min_conn_dist = 0.0
max_conn_dist = 9999.9

if dist_constraint:
    max_conn_dist = 300.0 #9999.9# Distance constraint for all cells

# Load synapse dictionaries
synapses.load()
syn = synapses.syn_params_dicts('../../components/synaptic_models')
###################################################################################
####################################Pyr Type A#####################################

# Add a population of numPN_A nodes (all of which share model_type, dynamics_params, etc.)
net.add_nodes(N=numPN_A, pop_name='PyrA',
              mem_potential='e',
              model_type='biophysical',
              #model_template='hoc:feng_typeA',#Ben's model
              model_template='hoc:Cell_Af',
              morphology=None)


##################################################################################
###################################Pyr Type C#####################################

# Add a population of numPN_A nodes (all of which share model_type, dynamics_params, etc.)
net.add_nodes(N=numPN_C, pop_name='PyrC',
              mem_potential='e',
              model_type='biophysical',
              #model_template='hoc:feng_typeC',#Ben's model
              model_template='hoc:Cell_Cf',
              morphology=None)

#################################################################################
###########################Fast - spiking PV ints################################

# Add a population of numBask nodes
net.add_nodes(N=numBask, pop_name='Bask',
              mem_potential='e',
              model_type='biophysical',
              #model_template='hoc:basket',#Ben's model
              model_template='hoc:InterneuronCellf',
              morphology=None)


################################################################################
############################# BACKGROUND INPUTS ################################

# External inputs
thalamus = NetworkBuilder('mthalamus')
thalamus.add_nodes(N=numPN_A+numPN_C+numBask,
                   pop_name='tON',
                   pop_group='mthalamus',
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

def syn_uniform_delay_section(source, target, sec_id=0, sec_x=0.9, mean=0.5,std=1):
    return np.random.uniform(mean,std), sec_id, sec_x


##########################################################################
######################### BACKGROUND INPUT ###############################


dynamics_file='BG2PNe_feng_min_initw.json'

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


dynamics_file = 'BG2PNi_feng_min_initw.json'

conn = net.add_edges(source=thalamus.nodes(), target=net.nodes(pop_name='Bask'),
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

#
#print("External nodes and edges built")
t_sim = 1000.0

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
		network_dir='./network',
		tstop=t_sim, dt = 0.05,
		report_vars = ['v'],
                v_init = -70.0,
                celsius = 31.0,
		spikes_inputs=[('mthalamus',   # Name of population which spikes will be generated for
                                'mthalamus_spikes.h5')],
		components_dir='../../components',
		compile_mechanisms=False)


#build_input.build_input(t_sim,numPN_A,numPN_C,2)

def build_h5(filename='mthalamus_spikes.h5'):
    """
    >>> f.keys()
    <KeysViewHDF5 ['spikes']>
    >>> f['spikes'].keys()
    <KeysViewHDF5 ['mthalamus']>
    >>> f['spikes']['mthalamus'].keys()
    <KeysViewHDF5 ['node_ids', 'timestamps']>
    >>> f['spikes']['mthalamus'].keys()

    """

    import h5py,numpy as np
    hf = h5py.File(filename,'w')

    gbla = hf.create_group('spikes/mthalamus')

    #timestamps = np.array([250, 250, 250, 500, 500, 500, 550, 550, 500, 600, 600, 600, 900, 900, 900])
    timestamps = np.array([450, 450, 450, 500, 500, 500, 550, 550, 550, 600, 600, 600, 650, 650, 650])


    nodes = np.array([0,1,2,0,1,2,0,1,2,0,1,2,0,1,2])

    gbla.create_dataset('timestamps',data=timestamps)
    gbla.create_dataset('node_ids',data=nodes)

    hf.close()

build_h5()
