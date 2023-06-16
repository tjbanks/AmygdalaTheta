import glob
import json
import os

from bmtk.simulator.bionet.pyfunction_cache import add_synapse_model
from neuron import h
import numpy as np

def lognorm(mean,std):
    mean = float(mean)
    std = float(std)
    mean_ = np.log(mean) - 0.5 * np.log((std/mean)**2+1)
    std_ = np.sqrt(np.log((std/mean)**2 + 1))
    return float(np.random.lognormal(mean_,std_))

def Bg2Pyr(syn_params, sec_x, sec_id):
    """Create a bg2pyr synapse
    :param syn_params: parameters of a synapse
    :param sec_x: normalized distance along the section
    :param sec_id: target section
    :return: NEURON synapse object
    """

    lsyn = h.bg2pyr(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])
    return lsyn

def bginh(syn_params, sec_x, sec_id):

    lsyn = h.bginh(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])
    
    if syn_params.get("AlphaTmax_gaba"):
        lsyn.AlphaTmax_gaba = float(syn_params['AlphaTmax_gaba'])
    if syn_params.get("Beta_gaba"):
        lsyn.Beta_gaba = float(syn_params['Beta_gaba'])

    return lsyn

def interD2interD_STFD(syn_params, sec_x, sec_id):

    lsyn = h.interD2interD_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn


def interD2pyrD_STFD(syn_params, sec_x, sec_id):

    lsyn = h.interD2pyrD_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn

def pyrD2interD_STFD(syn_params, sec_x, sec_id):
    
    lsyn = h.pyrD2interD_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn

def pyrD2pyrD_STFD(syn_params, sec_x, sec_id):
    
    lsyn = h.pyrD2pyrD_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])
        
    return lsyn

def pyrD2interD_P2SOM_STFD(syn_params, sec_x, sec_id):
    
    lsyn = h.pyrD2interD_P2SOM_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn

def interD2pyrD_SOM2P_STFD(syn_params, sec_x, sec_id):

    lsyn = h.interD2pyrD_SOM2P_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn

def pyrD2interD_P2CR_STFD(syn_params, sec_x, sec_id):

    lsyn = h.pyrD2interD_P2CR_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn

def interD2pyrD_CR2P_STFD(syn_params, sec_x, sec_id):

    lsyn = h.interD2pyrD_CR2P_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn

def interD2interD_SOMPV_STFD(syn_params, sec_x, sec_id):

    lsyn = h.interD2interD_SOMPV_STFD(sec_x, sec=sec_id)
    if syn_params.get('bACH'):
        if hasattr(lsyn, 'bACH'):
            lsyn.bACH = float(syn_params['bACH'])

    if syn_params.get('initW'):
        lsyn.initW = float(syn_params['initW'])
    elif syn_params.get('initW_lognormal_mean') and syn_params.get('initW_lognormal_std'):
        lsyn.initW = lognorm(syn_params['initW_lognormal_mean'],syn_params['initW_lognormal_std'])

    return lsyn


def load():
    add_synapse_model(Bg2Pyr, 'bg2pyr', overwrite=False)
    add_synapse_model(Bg2Pyr, overwrite=False)
    add_synapse_model(bginh, 'bginh', overwrite=False)
    add_synapse_model(bginh, overwrite=False)
    add_synapse_model(interD2interD_STFD, 'interD2interD_STFD', overwrite=False)
    add_synapse_model(interD2interD_STFD, overwrite=False)
    add_synapse_model(interD2pyrD_STFD, 'interD2pyrD_STFD', overwrite=False)
    add_synapse_model(interD2pyrD_STFD, overwrite=False)
    add_synapse_model(pyrD2interD_STFD, 'pyrD2interD_STFD', overwrite=False)
    add_synapse_model(pyrD2interD_STFD, overwrite=False)
    add_synapse_model(pyrD2pyrD_STFD, 'pyrD2pyrD_STFD', overwrite=False)
    add_synapse_model(pyrD2pyrD_STFD, overwrite=False)

    #SOM
    add_synapse_model(pyrD2interD_P2SOM_STFD, 'pyrD2interD_P2SOM_STFD', overwrite=False)
    add_synapse_model(pyrD2interD_P2SOM_STFD, overwrite=False)
    add_synapse_model(interD2pyrD_SOM2P_STFD, 'interD2pyrD_SOM2P_STFD', overwrite=False)
    add_synapse_model(interD2pyrD_SOM2P_STFD, overwrite=False)
    add_synapse_model(interD2interD_SOMPV_STFD, 'interD2interD_SOMPV_STFD', overwrite=False)
    add_synapse_model(interD2interD_SOMPV_STFD, overwrite=False)

    #CR
    add_synapse_model(pyrD2interD_P2CR_STFD, 'pyrD2interD_P2CR_STFD', overwrite=False)
    add_synapse_model(pyrD2interD_P2CR_STFD, overwrite=False)
    add_synapse_model(interD2pyrD_CR2P_STFD, 'interD2pyrD_CR2P_STFD', overwrite=False)
    add_synapse_model(interD2pyrD_CR2P_STFD, overwrite=False)


    return

def syn_params_dicts(syn_dir='components/synaptic_models'):
    """
    returns: A dictionary of dictionaries containing all
    properties in the synapse json files
    """
    files = glob.glob(os.path.join(syn_dir,'*.json'))
    data = {}
    for fh in files:
        with open(fh) as f:
            data[os.path.basename(fh)] = json.load(f) #data["filename.json"] = {"prop1":"val1",...}
    return data
