import numpy as np
import pandas as pd
import random
import math

all_synapses = pd.DataFrame([],columns=['source_gid','target_gid'])

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
    #x = float(x_end - x_start)#/1000
    #y = float(y_end - y_start)#/1000
    #z = float(z_end - z_start)#/1000
    max_dist = 2.425#np.sqrt(x**2 + y**2 + z**2)
    max_delay=2.425 #////define maximum delay,ms

    x_ind,y_ind,z_ind = 0,1,2
    
    dx = target['positions'][x_ind] - source['positions'][x_ind]
    dy = target['positions'][y_ind] - source['positions'][y_ind]
    dz = target['positions'][z_ind] - source['positions'][z_ind]

    dist = np.sqrt(dx**2 + dy**2 + dz**2)/1000
    
    del_fluc = np.random.uniform(-0.1,0.1)

    delay=(dist/max_dist)*max_delay+min_delay+del_fluc+dt

    return delay

def get_target_sec_id(source,target):
    # this simplified target selector could be replaced with something more 
    # robust using a lookup table for more complicated cell types
    
    if source['pop_name'] == 'PyrA' or source['pop_name'] == 'PyrC':
        return 1 # Target Dendrites
    elif source['pop_name'] == 'PV':
        return 0 # Target Soma
    elif source['pop_name'] == 'SOM':
        return 0 # Target Soma
    elif source['pop_name'] == 'CR':
        if target['pop_name'] == 'PyrA' or target['pop_name'] == 'PyrC':
            return 0 # Target Soma
        else:
            return 1 # Target Denditrites
    elif source['model_type'] == 'virtual':
        return 1 # Target Dendrites
    else: # We really don't want a default case so we can catch errors
        #return 0 # Target Soma by default
        import pdb;pdb.set_trace()

def syn_dist_delay_feng_section(source, target, sec_id=None, sec_x=0.9):

    if sec_id is None: # allows for overriding
        sec_id = get_target_sec_id(source, target)

    return syn_dist_delay_feng(source, target), sec_id, sec_x

def syn_uniform_delay_section(source, target, sec_id=None, sec_x=0.9, mean=0.5,std=1):
    
    if sec_id is None: # allows for overriding
        sec_id = get_target_sec_id(source, target)

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

def points_in_cylinder(pt1, pt2, r, q):
    #https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
    vec = pt2 - pt1
    const = r * np.linalg.norm(vec)
    c1 = np.dot(q - pt1, vec) >= 0
    c2 = np.dot(q - pt2, vec) <= 0
    c3 = np.linalg.norm(np.cross(q - pt1, vec),axis=1) <= const 
    return c1 & c2 & c3

def syn_percent_o2a(source,targets,p,track_list=None,no_recip=False, autapses_allowed=False, angle_dist=False,min_dist=0, max_dist=300, angle_dist_radius=100):
    """
    track_list: supply a list to append and track the synapses with
    one to all connector for increased speed.
    returns a list of len(targets) where the value is the number of synapses at the index=cellid
    if no_recip is set to true, any connection that has been made previously and is reciprical won't be considered
    """
    original_p = p

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
    
    if not autapses_allowed:
        available = available[~np.isin(available,sid)] #remove self from possible targets
    
    if no_recip:
        recur = [i['source_gid'] for i in track_list if i['target_gid'] == sid]
        available = available[~np.isin(available,recur)]

    if source.get('positions') is not None:
        src_pos = np.array(source['positions'])
        trg_pos = np.array([target['positions'] for target in targets])

        euclid_dist = np.linalg.norm(trg_pos - src_pos,axis=1)
        euclid_mask_dist = np.array((euclid_dist < max_dist) & (euclid_dist > min_dist))
        euclid_dist_available = np.argwhere(euclid_mask_dist).T[0] # Convert to indicies        
        default_dist_available = euclid_dist_available

        if angle_dist:
            # Checks if points are inside a cylinder
            src_angle_x = np.array(source['rotation_angle_zaxis'])
            src_angle_y = np.array(source['rotation_angle_yaxis'])

            vec_pos = np.array([np.cos(src_angle_x), np.sin(src_angle_y), np.sin(src_angle_x)])
            pt1 = src_pos + vec_pos * min_dist
            pt2 = src_pos + vec_pos * max_dist # Furthest point (max dist away from position of cell)

            mask_dist = points_in_cylinder(pt1, pt2, angle_dist_radius, trg_pos)
            angle_dist_available = np.argwhere(mask_dist).T[0] # Convert to indicies
            default_dist_available = angle_dist_available
        
        # We now want to reduce the number available by distance 
        available = available[np.isin(available,default_dist_available)]

    # of those remaining we want p% chosen from all within the distance around in a sphere, no matter angle
    n = int(len(euclid_dist_available)*p)
    extra = 1 if random.random() < (p*100 - math.floor(p*100)) else 0
    n = n + extra #This will account for half percentages

    if n > len(available):
        n = len(available)

    chosen = np.random.choice(available,size=n,replace=False)
    mask = np.isin(tids,chosen)
    syns[mask] = 1

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

def syn_percent_o2a_old(source,targets,p,track_list=None,no_recip=False, angle_dist=False,min_dist=0, max_dist=300, angle_dist_radius=100, warn=False):
    """
    track_list: supply a list to append and track the synapses with
    one to all connector for increased speed.
    returns a list of len(targets) where the value is the number of synapses at the index=cellid
    if no_recip is set to true, any connection that has been made previously and is reciprical won't be considered
    """
    original_p = p

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

    dist = None    
    mask_dist = None
 
    if source.get('positions') is not None:
        src_pos = np.array(source['positions'])
        trg_pos = np.array([target['positions'] for target in targets])
    
        #if angle_dist: #https://github.com/latimerb/SPWR_BMTK2/blob/master/build_network.py#L148-L176
        #    """
        #    'finding the perpendicular distance from a three dimensional vector ... the goal was simply 
        #     to calculate the perpendicular distance of the target cell from the source cell’s direction vector... 
        #     the Euclidean distance would be the hypotenuse of that right triangle so the 
        #     perpendicular distance should be the opposite side.
        #     the way I was thinking about it was to imagine a cylinder with its center around the [directional] vector
        #     ... and only cells that fall in the cylinder are eligible for connection' - Per Ben
        #    """
        #    vec_pos = np.array([np.cos(src_angle_x), np.sin(src_angle_y), np.sin(src_angle_x)])
        #    dist = np.linalg.norm(np.cross((trg_pos - src_pos), (trg_pos - vec_pos)),axis=1) / np.linalg.norm((vec_pos - src_pos))
        if angle_dist:
            # Checks if points are inside a cylinder
            src_angle_x = np.array(source['rotation_angle_zaxis'])
            src_angle_y = np.array(source['rotation_angle_yaxis'])
            
            vec_pos = np.array([np.cos(src_angle_x), np.sin(src_angle_y), np.sin(src_angle_x)])
            pt1 = src_pos + vec_pos * min_dist
            pt2 = src_pos + vec_pos * max_dist # Furthest point (max dist away from position of cell)
            
            mask_dist = points_in_cylinder(pt1, pt2, angle_dist_radius, trg_pos)            
            
        else: 
            dist = np.linalg.norm(trg_pos - src_pos,axis=1)
            mask_dist = np.array((dist < max_dist) & (dist > min_dist))

        # Since we cut down on the number of available cells due to distance constraints we need to scale up the p
        avg_connected = p*len(targets)
        # new p
        num_available = len(np.where(mask_dist==True)[0])
        if num_available:
            p = avg_connected/num_available
        else:
            p = 1.1

        if p > 1:
            p = 1
            if not angle_dist and warn:
                sorted_dist = np.sort(dist)
                minimum_max_dist = sorted_dist[int(avg_connected)]
                print("Warning: distance constraint (max_dist:" + str(max_dist) + ") between gid " + str(sid) + " and target gids " +
                    str(min(tids)) + "-" + str(max(tids)) + " prevented " + str(original_p*100) + "% overall connectivity. " +
                    "To achive this connectivity, max_dist would have needed to be greater than " + str(minimum_max_dist))

    # of those remaining we want p% chosen
    n = int(len(tids)*p)
    extra = 1 if random.random() < (p*100 - math.floor(p*100)) else 0
    n = n + extra #This will account for half percentages

    if len(available) < n:
        n = len(available)

    chosen = np.random.choice(available,size=n,replace=False)
    mask = np.isin(tids,chosen)
    
    if mask_dist is not None:
        mask = mask & mask_dist
        chosen = np.where(mask==True)[0]

    syns[mask] = 1
 
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
