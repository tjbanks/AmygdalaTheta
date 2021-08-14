from bmtk.analyzer import nodes_table
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pdb
from mpl_toolkits.mplot3d import Axes3D
import h5py
import networkx as nx

num_PN = 164
num_PV = 37
num_Chn = 4
# Load files
# Node information
nodes_info = h5py.File('network/SPWR_biophysical_nodes.h5')
# Connection information
edges_info = h5py.File('updated_conns/SPWR_biophysical_SPWR_biophysical_edges.h5')

# Make source and target IDs into numpy arrays
src_ids = np.array(edges_info['edges']['SPWR_biophysical_SPWR_biophysical']['source_node_id'])
tgt_ids = np.array(edges_info['edges']['SPWR_biophysical_SPWR_biophysical']['target_node_id'])

# Make connectivity matrix
conns = np.array([src_ids,tgt_ids]).T

# Get positions
#dset = np.array(nodes_info['nodes']['SPWR_biophysical']['positions'].values)

# Convert positions to a dataframe
#pos_df = pd.DataFrame(np.array(list(dset))[:,0:3])

# Weights
wgts = edges_info['edges']['SPWR_biophysical_SPWR_biophysical']['0']['syn_weight'][:]
wgts = np.reshape(wgts, (-1, 1)) 

# Connectivity/weight histograms
Pyr2PVconns = conns[(conns[:,0]<num_PN) & (conns[:,1]<num_PN+num_PV) & (conns[:,1]>num_PN),:]
Pyr2PVwgts = wgts[(conns[:,0]<num_PN) & (conns[:,1]<num_PN+num_PV) & (conns[:,1]>num_PN),:]

plt.hist(Pyr2PVwgts,bins=20)


# Divergence (Number of cells contacted by the presynaptic cell)
unqvals1, Pyr2PVdiverg = np.unique(Pyr2PVconns[:,0], return_counts=True)

# Convergence (Number of cells the postsynaptic cell receives)
unqvals2, Pyr2PVconverg = np.unique(Pyr2PVconns[:,1], return_counts=True)

convhist = plt.figure(figsize=(7,5))
plt.hist(Pyr2PVdiverg,bins=129)
plt.title('Pyr2PVdivergence')

divhist = plt.figure(figsize=(7,5))
plt.hist(Pyr2PVconverg,bins=500)
plt.title('Pyr2PVconvergence')


# Plot positions
#fig = plt.figure(figsize=(7,5))
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#ax.scatter(pos_df.iloc[num_PN+num_PV:num_PN+num_PV+num_Chn-1,0],
#        pos_df.iloc[num_PN+num_PV:num_PN+num_PV+num_Chn-1,1],
#        pos_df.iloc[num_PN+num_PV:num_PN+num_PV+num_Chn-1,2],color='blue')
#ax.scatter(pos_df.iloc[0:num_PN-1,0],
#        pos_df.iloc[0:num_PN-1,1],
#        pos_df.iloc[0:num_PN-1,2],color=[1,0,0])
#ax.scatter(pos_df.iloc[num_PN:num_PN+num_PV-1,0],
#        pos_df.iloc[num_PN:num_PN+num_PV-1,1],
#        pos_df.iloc[num_PN:num_PN+num_PV-1,2],color=[0,1,0])
#plt.title('entire network')
#
#fig2 = plt.figure(figsize=(7,5))
#ax2 = fig2.add_subplot(1, 1, 1, projection='3d')
#ax2.scatter(pos_df.iloc[num_PN+num_PV:num_PN+num_PV+num_Chn-1,0],
#        pos_df.iloc[num_PN+num_PV:num_PN+num_PV+num_Chn-1,1],
#        pos_df.iloc[num_PN+num_PV:num_PN+num_PV+num_Chn-1,2],color='blue')
#plt.title('chandeliers')
#
#
#fig3 = plt.figure(figsize=(7,5))
#ax3 = fig3.add_subplot(1, 1, 1, projection='3d')
#ax3.scatter(pos_df.iloc[0:num_PN-1,0],
#        pos_df.iloc[0:num_PN-1,1],
#        pos_df.iloc[0:num_PN-1,2],color=[1,0,0])
#plt.title('PNs')
#
#
#fig4 = plt.figure(figsize=(7,5))
#ax4 = fig4.add_subplot(1, 1, 1, projection='3d')
#ax4.scatter(pos_df.iloc[0:num_PN-1,0],
#        pos_df.iloc[0:num_PN-1,1],
#        pos_df.iloc[0:num_PN-1,2],color=[1,0,0])
#ax4.set_xlim3d(600,700)
#ax4.set_ylim3d(600,700)
#ax4.set_zlim3d(200,300)
#plt.title('PNs')
#
#
#fig5 = plt.figure(figsize=(7,5))
#ax5 = fig5.add_subplot(1, 1, 1, projection='3d')
#ax5.scatter(pos_df.iloc[num_PN:num_PN+num_PV-1,0],
#        pos_df.iloc[num_PN:num_PN+num_PV-1,1],
#        pos_df.iloc[num_PN:num_PN+num_PV-1,2],color=[0,1,0])
#plt.title('PVs')

plt.show()
