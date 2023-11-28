import os,pathlib,sys
from bmtk.simulator import bionet
import json
import numpy as np
import synapses
import warnings

import argparse

from neuron import h
pc = h.ParallelContext()    # object to access MPI methods

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

    with open(config_file, 'r') as json_file:
        conf_dict = json.load(json_file)
        # modify the output dir if we specify it OUTPUT_DIR=./cases_paper/case1/run2 sbatch h1_.....sh
        if os.environ.get("OUTPUT_DIR"):
            output_dir = os.path.abspath(os.environ.get('OUTPUT_DIR'))
            pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
            print(f"Output directory updated to {output_dir}")
            conf_dict['manifest']['$OUTPUT_DIR'] = output_dir
        if os.environ.get("TSTOP"):
            tstop = float(os.environ.get("TSTOP"))
            print(f"tstop modified to {tstop}")
            conf_dict['run']['tstop'] = tstop

        if coreneuron:
            conf = corebmtk.Config.from_dict(conf_dict, validate=True)        
        else:
            conf = bionet.Config.from_dict(conf_dict, validate=True)
    
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
        ach_receptive_property = "ACH"
        ach_recpetive_property_on_value = conf.get("ACH_level", 1) #default no ACH delivered
        local_cells = graph.get_local_cells()
        print(f"Rank {pc.id()} modifying {len(local_cells)} local cells")
        for cell_id, cell in local_cells.items():
            if cell_id in ach_receptive_cells:
                cell_connections = cell.connections()
                for connection in cell_connections:
                    syn = connection._syn
                    if hasattr(syn, ach_receptive_property):
                        setattr(syn, ach_receptive_property, ach_recpetive_property_on_value)
                        num_modified_synapses += 1
                    elif 'bg2' in syn.__repr__() or 'bginh' in syn.__repr__():
                        total_modified_synapses -= 1 # not changing bg2pyr
                    else:
                        print(f"{syn} not modified for ACH")
                        total_modified_synapses -= 1
                        #import pdb;pdb.set_trace() # we missed one debug
                total_modified_synapses += len(cell_connections)

        print(f"ACH turned on for {len(ach_receptive_cells)} muscarinic (synaptic) receptors")
        print(f"{num_modified_synapses}/{total_modified_synapses} ACH synapses modified")
        print(f"{ach_receptive_property} set to value {ach_recpetive_property_on_value}")
        
        # 3. Turn on nicotinic receptors per cell
        nic_cells_turned_on = 0
        for cell_id, cell in local_cells.items():
            if cell_id in ach_receptive_cells:
                hobj = cell.hobj
                if hasattr(hobj, 'activate_ach'):
                    hobj.activate_ach()
                    nic_cells_turned_on += 1
        print(f"Rank {pc.id()}: {nic_cells_turned_on}/{len(ach_receptive_cells)} nicotinic (intrinsic) receptors turned on")

        no_pv = False
        no_som = False 
        no_cr = False
        if os.environ.get("NO_PV"):
            print(f"Ablation set for PV")
            no_pv = True
        if os.environ.get("NO_SOM"):
            print(f"Ablation set for SOM")
            no_som = True
        if os.environ.get("NO_CR"):
            print(f"Ablation set for CR")
            no_cr = True
        if no_pv or no_som or no_cr:
            for cell_id, cell in local_cells.items():
                cell_type = ''
                if hasattr(cell.hobj, 'type'):
                    cell_type = cell.hobj.type
                if (cell_type == 'InterneuronCellf' and no_pv) or (cell_type == 'SOM_Cell' and no_pv) or (cell_type == 'CR_Cell' and no_cr):
                    cell_connections = cell.connections()
                    for connection in cell_connections:
                        syn = connection._syn
                        if hasattr(syn, 'initW'):
                            syn.initW = 0.000001
                            print(f"Setting cell {cell_id} syn.initW to 0.000001")
        pc.barrier()

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
