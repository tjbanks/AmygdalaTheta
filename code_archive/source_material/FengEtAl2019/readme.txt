This is the readme for the model associated with the journal paper:

Feng F, Headley DB , Amir A, Kanta V, Chen Z, Pare D, Nair S (2019)
Gamma oscillations in the basolateral amygdala: biophysical mechanisms
and computational consequences eNeuro.

This model was contributed by Feng Feng.

This is the first large-scale biologically realistic model of the
basolateral amygdala (BL) with parameters constrained by published
literature. The model matches closely with numerous in vivo testing
results, including average neuron firing rate, gamma oscillation
detected on LFPs, and firing phase entrainments to gamma oscillation
with spatial gradients.

This model has been developed using the NEURON simulator (Hines and
Carnavale, 2008).

Summary of files: 

DATA INPUT FILES:
spikesmatrix_op:   Stores spikes for each external connection (row); it needs to be used together with spikes_ind file
E2P:             Stores P cells gids receiving each external spike train (row)
E2I:             Stores FSI cells gids receiving each external spike train (row)
Cell_type.txt:      This file indicates the type of cell for the PNs and FSIs. The types are defined in the main file. 1 for A PNs, 10 for C PNs and 100 for FSIs.
active_syn_op:     Stores pre-cell gids for each post-cell (row), for internal connections; it needs to be used together with active_syn_ind.
active_syn_GAP_op: Stores connected FSI gids via gap couplings; it needs to be used together with active_syn_ind active_GAP_ind.
GAPid:            Indexes only used for gap junction transfer. Indexes here has nothing to do with neurons’ gid.
sim_length:        Simulation time length.
location.txt:        Stores 3D coordinates (um) for each cell
oritation.txt:        Stores dendrite orientations for each cell
elec_coords.txt:     List electrode coordinates (um) to calculate LFPs. If interesting in multiple electrode recording, need to generate multiple electrode coordinates to write into this file. 
NM.txt:           This file indicates whether the cell has any neuromodulator receptors present. It can be of three types, DA, NE and DANE. (has not been used in this study).
NEURON FILES:
- BL_main.hoc:       Main file to run the simulations
- function_ConnectInputs_invivo_op.hoc:  Establish external connections and feed predefined external spikes to BL network with adjustable synapse weight parameters.
- function_ConnectInternal_simplify_online_op.hoc: Establish internal connections with adjustable synapse weight parameters.
- function_ConnectInternal_gj_simplify: Establish gap coupling between FSIs with adjustable coupling strength parameters.
- function_ConnectInternal_simplify_online_op.hoc: Establish internal connections with adjustable synapse weight parameters.
- function_ConnectTwoCells.hoc: Procedure used to connect two internal neurons.
- function_calcconduc: Function used to calculate conductances for being used to calculate LFPs. 

OUTPUT FILES:

data saves spiking time history of all cells
\LFPs\LFP_elec_0: saves calculated LFP from center electrode.

BEFORE YOU RUN:

Make sure to compile the mod files and also make an empty folder named
LFPs in the directory. NEURON will save the calculated LFP files into
this directory.
