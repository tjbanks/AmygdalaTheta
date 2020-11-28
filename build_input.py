import sys

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

def build_input(t_sim, numPN_A = 640, numPN_C=260, numBask = 100):
    print("Building input for " + str(t_sim) + " (ms)")
    psg = PoissonSpikeGenerator(population='mthalamus')
    psg.add(node_ids=range(numPN_A+numPN_C),  # Have nodes to match mthalamus
        firing_rate=2.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
        times=(0.0, t_sim/1000.0))    # Firing starts at 0 s up to 3 s
    psg.to_sonata('mthalamus_spikes.h5')

    psg = PoissonSpikeGenerator(population='exc_bg_bask')
    psg.add(node_ids=range(numBask),  # Have nodes to match mthalamus
        firing_rate=2.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
        times=(0.0, t_sim/1000.0))    # Firing starts at 0 s up to 3 s
    psg.to_sonata('exc_bg_bask_spikes.h5')
    print("Done")

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        build_input(int(sys.argv[-1]))
    else:
        build_input(100000)
