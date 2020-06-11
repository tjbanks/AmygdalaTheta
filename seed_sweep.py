import build_network
import update_config
import run_bionet

import os
import sys

import numpy as np
import random

import logging
logging.getLogger().setLevel(logging.INFO)

def run_seed(seed=2020,experiment="0",simulation_config="simulation_config.json",skip_build=False):

    logging.info("Running seed " + str(seed))
    network_dir = str(seed)+"_network"
    output_dir = str(seed)+"_output"
    sim_config = str(seed)+"_simulation_config.json"
    net_config = str(seed)+"_network_config.json"

    logging.info("Building Network")
    if not skip_build:
        build_network.build(output=network_dir,seed=seed,experiment=experiment)

    logging.info("Creating config files")
    update_config.update_config(simulation_config,
            new_simulation_config=sim_config,
            new_network_config=net_config,
            output=output_dir,network=network_dir)

    logging.info("Running simulation")
    try:
        run_bionet.run(sim_config,seed=seed,output=output_dir,experiment=experiment)
    except SystemExit:
        pass

if __name__ == '__main__':
    if len(sys.argv)==5:#Really terrible, fix sometime
        run_seed(int(sys.argv[1]),sys.argv[2],sys.argv[3],skip_build = True)
    elif __file__ != sys.argv[-1]:
        run_seed(int(sys.argv[1]),sys.argv[2],sys.argv[3])
    else:
        run_seed()
