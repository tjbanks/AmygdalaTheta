import json
import glob
import os

def update_config(simulation_config, new_simulation_config=None, new_network_config=None, output=None, network=None):
    """ 
    Assume we're using manifest $BASE_DIR 
    TODO: Fix later to generalize
    """
    
    edited = False
    
    with open(simulation_config, "r") as json_file:
        data = json.load(json_file)
        
    sim_manifest = data["manifest"]
    
    if output:
        data["manifest"]["$OUTPUT_DIR"] = "$BASE_DIR/" + output
        edited = True
            
    if network:
        net = dereference(sim_manifest,data["network"])
        
        with open(net, "r") as json_file:
            net_data = json.load(json_file)
        
        net_manifest = net_data["manifest"]
        
        net_data["manifest"]["$NETWORK_DIR"] = "$BASE_DIR/" + network
        
        net_dir = dereference(net_manifest, net_data["manifest"]["$NETWORK_DIR"])
        
        node_files = glob.glob(net_dir + "/*_nodes.h5")
        node_types_files = glob.glob(net_dir + "/*_node_types.csv")
        edge_files = glob.glob(net_dir + "/*_edges.h5")
        edge_types_files = glob.glob(net_dir + "/*_edge_types.csv")
        
        nodes = []
        edges = []
        
        for node in node_files:
            #A little unsafe but we can assume build network did its job
            new_nodes = {}
            _,node = os.path.split(node)
            
            new_nodes["nodes_file"] = "$NETWORK_DIR/"+node
            new_nodes["node_types_file"] = "$NETWORK_DIR/"+node.replace("_nodes.h5","")+"_node_types.csv"
            nodes.append(new_nodes)
            
        for edge in edge_files:
            new_edges = {}
            _,edge = os.path.split(edge)
            
            new_edges["edges_file"] = "$NETWORK_DIR/"+edge
            new_edges["edge_types_file"] = "$NETWORK_DIR/"+edge.replace("_edges.h5","")+"_edge_types.csv"
            edges.append(new_edges)
        
        net_data["networks"]["nodes"] = nodes
        net_data["networks"]["edges"] = edges
        
        if new_network_config:
            net = new_network_config
            
        with open(net,"w") as json_file:
            json_file.write(json.dumps(net_data,indent=2, separators=(',', ': ')))
            
    if new_network_config:
        data["network"] = new_network_config
    
    if new_simulation_config:
        simulation_config = new_simulation_config
        
    if edited:
        with open(simulation_config,"w") as json_file:
            json_file.write(json.dumps(data,indent=2, separators=(',', ': ')))
            
            
def dereference(manifest, str_):
    new_str = str_
    for prop in manifest:
        new_str = new_str.replace(prop,manifest[prop])
    return new_str