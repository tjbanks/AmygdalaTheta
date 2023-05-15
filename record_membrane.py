import neuron
from neuron import h
import numpy as np
import pandas as pd

def new_cell(cell_name):
    # Get the cell by name
    invoke_cell = getattr(h, cell_name)
    cell = invoke_cell()
    return cell

def main(hoc_files, mechanisms_dir, cell_names, currents=[], i_clamp_dur=300, i_clamp_delay=10, v_init=-70, tstop=400):
        # Load mechs
        neuron.load_mechanisms(mechanisms_dir)
        # Load the standard hoc file
        h.load_file('stdrun.hoc')
        for hoc_file in hoc_files:
            # Load the .hoc file
            h.load_file(hoc_file)

        cells = []
        mem_potentials = []
        i_clamps = []

        if type(cell_names) == str:
              cell_names = [cell_names]*len(currents)

        for cell_name, current in zip(cell_names, currents):
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
        
        steps_per_ms = 10
        h.steps_per_ms = steps_per_ms
        h.dt = 1 / steps_per_ms
        h.tstop = tstop
        h.v_init = v_init
        h.run()

        potentials = [list(v) for v in mem_potentials]
        potentials = np.array(potentials).T
        potentials_df = pd.DataFrame(potentials,columns=cell_names)
        potentials_df.to_csv('potentials.csv', index=False)


if __name__ == '__main__':
        currents = [0.16,0.16,0.16,0.16,0.16]
        main(['./components_homogenous/templates/feng.hoc','./components_homogenous/templates/SOM.hoc'],
             './components_homogenous/mechanisms/modfiles',
             ['Cell_Af','Cell_Cf','InterneuronCellf','SOM_Cell','CR_Cell'],
             currents)