import numpy as np
import sys

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator


def lognorm_fr_list(n,m,s):
    mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
    std = np.sqrt(np.log((s/m)**2 + 1))
    return [np.random.lognormal(mean,std) for i in range(n)]

#def build_input(t_sim, numPN_A = 640, numPN_C=260, numBask = 100):
def build_input(t_sim, numPN_A = 569, numPN_C=231, numBask = 93, numSOM=51, numCR=56):
    print("Building input for " + str(t_sim) + " (ms)")
    psg = PoissonSpikeGenerator(population='vpsi_pyr')
    psg.add(node_ids=range(numPN_A+numPN_C),  # Have nodes to match mthalamus
    #    firing_rate=2.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
    firing_rate=lognorm_fr_list(numPN_A+numPN_C,2,1),    
    times=(0.0, t_sim/1000.0))    # Firing starts at 0 s up to 3 s
    psg.to_sonata('vpsi_pyr_spikes.h5')

    psg = PoissonSpikeGenerator(population='vpsi_pv')
    psg.add(node_ids=range(numBask),  # Have nodes to match mthalamus
    #    firing_rate=2.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
    firing_rate=lognorm_fr_list(numBask,2,1),
    times=(0.0, t_sim/1000.0))    # Firing starts at 0 s up to 3 s
    psg.to_sonata('vpsi_pv_spikes.h5')

    # THALAMUS

    psg = PoissonSpikeGenerator(population='thalamus_pyr')
    psg.add(node_ids=range(numPN_A+numPN_C),  # Have nodes to match 
    #    firing_rate=2.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
    firing_rate=lognorm_fr_list(numPN_A+numPN_C,2,1),    
    times=(0.0, t_sim/1000.0))    # Firing starts at 0 s up to 3 s
    psg.to_sonata('thalamus_pyr_spikes.h5')

    psg = PoissonSpikeGenerator(population='thalamus_som')
    psg.add(node_ids=range(numSOM),  # Have nodes to match 
    #    firing_rate=2.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
    firing_rate=lognorm_fr_list(numSOM,2,1),
    times=(0.0, t_sim/1000.0))    # Firing starts at 0 s up to 3 s
    psg.to_sonata('thalamus_som_spikes.h5')
    
    psg = PoissonSpikeGenerator(population='thalamus_cr')
    psg.add(node_ids=range(numCR),  # Have nodes to match 
    #    firing_rate=2.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
    firing_rate=lognorm_fr_list(numCR,2,1),
    times=(0.0, t_sim/1000.0))    # Firing starts at 0 s up to 3 s
    psg.to_sonata('thalamus_cr_spikes.h5')

    print("Done")



if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        build_input(int(sys.argv[-1]))
    else:
        build_input(100000)
