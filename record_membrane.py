import neuron
from neuron import h
import numpy as np

def new_cell(cell_name):
    # Get the cell by name
    invoke_cell = getattr(h, cell_name)
    cell = invoke_cell()
    return cell

def main(hoc_file, mechanisms_dir, cell_name, currents=[], i_clamp_dur=300, i_clamp_delay=10, v_init=-70, tstop=400):
        # Load mechs
        neuron.load_mechanisms(mechanisms_dir)
        # Load the standard hoc file
        h.load_file('stdrun.hoc')
        # Load the .hoc file
        h.load_file(hoc_file)

        cells = []
        mem_potentials = []
        i_clamps = []

        for current in currents:
            cell = new_cell(cell_name)
            recording_section = cell.soma[0](0.5)

            # Initialize vectors to track membrane voltage
            mem_potential = h.Vector()
            mem_potential.record(recording_section._ref_v)

            i_clamp = h.IClamp(recording_section)
            i_clamp.dur = i_clamp_dur
            i_clamp.amp = current
            i_clamp.delay = i_clamp_delay

            cells.append(cell)
            mem_potentials.append(mem_potential)
            i_clamps.append(i_clamp)
        
        steps_per_ms = 5000 / 600
        h.steps_per_ms = steps_per_ms
        h.dt = 1 / steps_per_ms
        h.tstop = tstop
        h.v_init = v_init
        h.run()

        potentials = [list(v) for v in mem_potentials]
        potentials = np.array(potentials).T


if __name__ == '__main__':
        currents = [0.1,0.2,0.3,0.4,0.5]
        main('./components_homogenous/templates/SOM.hoc',
             './components_homogenous/mechanisms/modfiles',
             'SOM_Cell',
             currents)