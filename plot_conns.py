import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pdb
from mpl_toolkits.mplot3d import Axes3D

f = h5py.File('./network/BLA_BLA_edges.h5','r')
g = h5py.File('./network/BLA_nodes.h5','r')

source = f['edges']['BLA_to_BLA']['source_node_id']
target = f['edges']['BLA_to_BLA']['target_node_id']

n_cells = 1000
scale = 1
pn_ids = np.arange(0*scale,799*scale)
itn_ids = np.arange(800*scale,893*scale)

positions = g['nodes']['BLA']['0']['positions']
y_angles = g['nodes']['BLA']['0']['rotation_angle_yaxis']
z_angles = g['nodes']['BLA']['0']['rotation_angle_zaxis']

source_pos = pd.DataFrame(data=np.concatenate((np.arange(0, n_cells).reshape(-1,1), 
						positions[:,:],
						np.where(np.arange(0,n_cells)<pn_ids[-1],'PN','ITN').reshape(-1,1)),axis=1),
			  columns=['source_ID','source_x','source_y','source_z','source_type'])

target_pos = pd.DataFrame(data=np.concatenate((np.arange(0,n_cells).reshape(-1,1), 
						positions[:,:],
						np.where(np.arange(0,n_cells)<pn_ids[-1],'PN','ITN').reshape(-1,1)),axis=1),
			  columns=['target_ID','target_x','target_y','target_z','target_type'])


conns = pd.DataFrame(data=np.concatenate((source[:].reshape(-1,1),target[:].reshape(-1,1)),axis=1),
                     columns=['source_ID','target_ID'])

conns = conns.join(source_pos,on='source_ID',rsuffix='_x')
conns = conns.join(target_pos,on='target_ID',rsuffix='_x')


#print(conns[(conns.source_x.astype(float)<400)&(conns.source_x.astype(float)>200)&
#      (conns.source_y.astype(float)<400)&(conns.source_y.astype(float)>200)&
#      (conns.source_z.astype(float)<400)&(conns.source_z.astype(float)>200)&
#      (conns.source_type=='PN')]['source_ID'])


d = np.sqrt(np.sum((conns[['source_x', 'source_y', 'source_z']].values.astype(float)
	 - conns[['target_x','target_y','target_z']].values.astype(float))**2,axis=1))

conns.loc[:,'distance'] = d

conns.loc[:,'source_x'] = conns.loc[:,'source_x'].astype(float)
conns.loc[:,'source_y'] = conns.loc[:,'source_y'].astype(float)
conns.loc[:,'source_z'] = conns.loc[:,'source_z'].astype(float)

conns.loc[:,'target_x'] = conns.loc[:,'target_x'].astype(float)
conns.loc[:,'target_y'] = conns.loc[:,'target_y'].astype(float)
conns.loc[:,'target_z'] = conns.loc[:,'target_z'].astype(float)

PN2PN = conns[(conns.source_type=='PN') & (conns.target_type=='PN')]
PN2ITN = conns[(conns.source_type=='PN') & (conns.target_type=='ITN')]
ITN2PN = conns[(conns.source_type=='ITN') & (conns.target_type=='PN')]
ITN2ITN = conns[(conns.source_type=='ITN') & (conns.target_type=='ITN')]



def plot_scatter(cell_to_plot):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    src_pos_x = float(conns[conns.source_ID==cell_to_plot]['source_x'].values[0])
    src_pos_y = float(conns[conns.source_ID==cell_to_plot]['source_y'].values[0])
    src_pos_z = float(conns[conns.source_ID==cell_to_plot]['source_z'].values[0])

    ax.scatter(src_pos_x, src_pos_y, src_pos_z, marker='^',color='blue')
    
    src_angle_y = y_angles[cell_to_plot]
    src_angle_x = z_angles[cell_to_plot]
    vec_line_pos = [np.cos(src_angle_x), np.sin(src_angle_y), np.sin(src_angle_x)] 
    vec_line_neg = [np.cos(src_angle_x), np.sin(-src_angle_y), np.sin(src_angle_x)] 

    max_dist = 300
    #ax.plot([vec_line_neg[0],5*src_pos_x], [vec_line_neg[1], 5*src_pos_y], [vec_line_neg[2], 5*src_pos_z],color='m')
    ax.quiver(src_pos_x,src_pos_y,src_pos_z,vec_line_pos[0],vec_line_pos[1],vec_line_pos[2],length=max_dist)
    src_pos = np.array([src_pos_x,src_pos_y,src_pos_z])
    
    pt2 = src_pos + np.array(vec_line_pos) * max_dist
    ax.scatter(pt2[0],pt2[1],pt2[2], marker='^',color='black')
 
    pn2itn_x = PN2ITN[PN2ITN.source_ID==cell_to_plot]['target_x'].values
    pn2itn_y = PN2ITN[PN2ITN.source_ID==cell_to_plot]['target_y'].values
    pn2itn_z = PN2ITN[PN2ITN.source_ID==cell_to_plot]['target_z'].values

    ax.scatter(pn2itn_x,pn2itn_y,pn2itn_z,marker='o',color='red')

    pn2pn_x = PN2PN[PN2PN.source_ID==cell_to_plot]['target_x'].values
    pn2pn_y = PN2PN[PN2PN.source_ID==cell_to_plot]['target_y'].values
    pn2pn_z = PN2PN[PN2PN.source_ID==cell_to_plot]['target_z'].values

    ax.scatter(pn2pn_x,pn2pn_y,pn2pn_z,marker='o',color='blue')

    ax.set_xlim(300,900)
    ax.set_ylim(300,900)
    ax.set_zlim(300,900)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

plot_scatter(99)
#plot_scatter(100)
plt.show()

# PN2PN figure #
plt.figure()
plt.subplot(1,2,1)
density, bins = np.histogram(PN2PN['distance'], normed=True, density=True)
unity_density = density / density.sum()
widths = bins[:-1] - bins[1:]
plt.bar(bins[1:],unity_density,width=widths)
plt.title('PN2PN - pdf')
plt.xlabel('distance (um)')

plt.subplot(1,2,2)
plt.bar(bins[1:],unity_density.cumsum(),width=widths)
plt.title('PN2PN - cdf')
plt.xlabel('distance (um)')

# PN2ITN figure #
plt.figure()
plt.subplot(1,2,1)
density, bins = np.histogram(PN2ITN['distance'], normed=True, density=True)
unity_density = density / density.sum()
widths = bins[:-1] - bins[1:]
plt.bar(bins[1:],unity_density,width=widths)
plt.title('PN2ITN - pdf')
plt.xlabel('distance (um)')

plt.subplot(1,2,2)
plt.bar(bins[1:],unity_density.cumsum(),width=widths)
plt.title('PN2ITN - cdf')
plt.xlabel('distance (um)')


plt.show()

