import os, sys
from bmtk.simulator import bionet
import numpy as np
import synapses
import warnings

import argparse

try:
    import corebmtk
except:
    pass

def run(config_file, coreneuron=False, gpu=False):

    warnings.simplefilter(action='ignore', category=FutureWarning)
    synapses.load()
    from bmtk.simulator.bionet.pyfunction_cache import add_weight_function

    def gaussianBL(edge_props, source, target):
        w0 = edge_props["syn_weight"]
        sigma = edge_props["weight_sigma"]
        return np.random.normal(w0, sigma, 1)
	
    def lognormal(edge_props, source, target):
        m = edge_props["syn_weight"]
        s = edge_props["weight_sigma"]
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return np.random.lognormal(mean, std, 1)

    add_weight_function(lognormal)
    add_weight_function(gaussianBL)


    if coreneuron:
        conf = corebmtk.Config.from_json(config_file, validate=True)        
    else:
        conf = bionet.Config.from_json(config_file, validate=True)

    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)

    # This fixes the morphology error in LFP calculation
    pop = graph._node_populations['BLA']
    for node in pop.get_nodes():
        node._node._node_type_props['morphology'] = node.model_template[1]

    if coreneuron:
        sim = corebmtk.CoreBioSimulator.from_config(conf, network=graph, gpu=gpu)
    else:
        sim = bionet.BioSimulator.from_config(conf, network=graph)
    
    # This calls insert_mechs() on each cell to use its gid as a seed
    # to the random number generator, so that each cell gets a different
    # random seed for the point-conductance noise
    cells = graph.get_local_cells()
    for cell in cells:
        cells[cell].hobj.insert_mechs(cells[cell].gid)
        pass
    
    # Setup ACh
    # add a line to config 
    #   "ACH_level": 0, 1, 2
    # to start executing this
    if conf.get('ACH_level'):
        # 1. Get the list of cells that have been pre-selected as ACH receptive
        # from the node_sets.json file
        ach_receptive_cells = conf['node_sets'].get('ach_cells')
        print(f"{len(ach_receptive_cells)} cells are set to be receptive to ACH")

        # 2. Iterate through each cell in the network and if it's receptive, alter the synapse
        total_modified_synapses = 0
        num_modified_synapses = 0
        ach_receptive_property = "bACH"
        ach_recpetive_property_on_value = conf.get("ACH_level", 0) #default no ACH delivered
        for cell_id, cell in graph._rank_node_ids.items():
            if cell_id in ach_receptive_cells:
                cell_connections = cell.connections()
                for connection in cell_connections:
                    syn = connection._syn
                    if hasattr(syn, ach_receptive_property):
                        setattr(syn, ach_receptive_property, ach_recpetive_property_on_value)
                        num_modified_synapses += 1
                total_modified_synapses += len(cell_connections)

        print(f"ACH turned on for {len(ach_receptive_cells)} muscarinic (synaptic) receptors")
        print(f"{num_modified_synapses}/{total_modified_synapses} ACH synapses modified")
        print(f"{ach_receptive_property} set to value {ach_recpetive_property_on_value}")
        
        # 3. Turn on nicotinic receptors per cell
        nic_cells_turned_on = 0
        for cell in ach_receptive_cells:
            hobj = cells[cell].hobj
            if hasattr(hobj, 'activate_ach'):
                hobj.activate_ach()
                nic_cells_turned_on += 1
        print(f"{nic_cells_turned_on}/{len(ach_receptive_cells)} nicotinic (intrinsic) receptors turned on")


    sim.run()
    bionet.nrn.quit_execution()


if __name__ == '__main__':

    #parser = argparse.ArgumentParser()
    #parser.add_argument("config")
    #parser.add_argument("coreneuron", action="store_true", help="Use CoreNeuron for simulation")
    #parser.add_argument("gpu", action="store_true", help="Use GPU for simulation")
    #args = parser.parse_args()
    
    #run(args.config, coreneuron=False, gpu=False)
    run(sys.argv[-1])
