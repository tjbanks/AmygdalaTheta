import numpy as np
import sys

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

scale=1

def lognorm_fr_list(n,m,s):
    mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
    std = np.sqrt(np.log((s/m)**2 + 1))
    return [np.random.lognormal(mean,std) for i in range(n)]

def build_poisson_input(population,node_ids,mean,std,output_h5,t_sim=15000):
    print('Building input for ' + population + "[" + str(len(node_ids)) + " cells at " + str(mean) + "(" + str(std) + ") Hz]")
    psg = PoissonSpikeGenerator(population=population)
    psg.add(node_ids=node_ids,  
    firing_rate=lognorm_fr_list(len(node_ids),mean,std),
    times=(0.0, t_sim/1000.0))  
    psg.to_sonata(output_h5)

def build_input(t_sim, numPN_A = 569, numPN_C=231, numPV = 93, numSOM=51, numCR=56,scale=1):
    
    # VPSI
    build_poisson_input(population='vpsi_inh',
                        node_ids=range((numPN_A+numPN_C+numPV)*scale),
                        mean=3,std=1,
                        output_h5='vpsi_inh_spikes_nonrhythmic.h5',
                        t_sim=t_sim)

    # THALAMUS
    build_poisson_input(population='thalamus_pyr',
                        node_ids=range((numPN_A+numPN_C)*scale),
                        mean=4,std=1,
                        output_h5='thalamus_pyr_spikes.h5',
                        t_sim=t_sim)

    build_poisson_input(population='thalamus_som',
                        node_ids=range((numSOM)*scale),
                        mean=2,std=1,
                        output_h5='thalamus_som_spikes.h5',
                        t_sim=t_sim) 
 
    build_poisson_input(population='thalamus_cr',
                        node_ids=range((numCR)*scale),
                        mean=2,std=1,
                        output_h5='thalamus_cr_spikes.h5',
                        t_sim=t_sim)

    print("Done")

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        build_input(int(sys.argv[-1]),scale=scale)
    else:
        build_input(100000,scale=scale)
