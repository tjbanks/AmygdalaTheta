# -*- coding: utf-8 -*-
import synapses

import os, sys, random
from bmtk.simulator import bionet
from bmtk.simulator.bionet.default_setters.cell_models import loadHOC
import numpy as np
from neuron import h

def run(config_file, seed=None, output=None, experiment=0):
    seed_default = 43 
    
    if seed:
        seed_default = seed
    
    random.seed(seed_default)

    #This is not needed in newer versions of BMTK
    #bionet.pyfunction_cache.add_cell_model(loadHOC, directive='hoc', model_type='biophysical')
    synapses.load()

    #BE SURE dL is appropriate (https://github.com/AllenInstitute/bmtk/blob/develop/bmtk/simulator/bionet/biocell.py#L80)
    # This will divide your segments up incorrectly and cause problems for hoc code
    conf = bionet.Config.from_json(config_file, validate=True)

    if output:
        conf.output['output_dir'] = output

    conf.build_env()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
   
    
    cvode = sim.h.CVode()
    def Commands():
        def _events():
            print('CVode Event Activated')
            return

            effectTime = 1
            cvode.event(effectTime , _events)

    #https://www.neuron.yale.edu/phpbb/viewtopic.php?f=2&t=2236
    fih = []
    fih.append(sim.h.FInitializeHandler(0, Commands))
    
    sim.run()
    bionet.nrn.quit_execution()

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('simulation_config.json')
