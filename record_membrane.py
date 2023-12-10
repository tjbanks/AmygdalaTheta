import neuron
from neuron import h
import numpy as np
import pandas as pd

def new_cell(cell_name):
    # Get the cell by name
    invoke_cell = getattr(h, cell_name)
    cell = invoke_cell()
    return cell

def main(cell_names, current, i_clamp_dur=400, i_clamp_delay=200, v_init=-70, tstop=500):

        cells = []
        mem_potentials = []
        i_clamps = []

        if type(cell_names) == str:
              cell_names = [cell_names]*len(currents)

        for cell_name in cell_names:
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

        print(f"Running with {int(current*1000)}pA...")
        h.run()

        potentials = [list(v) for v in mem_potentials]
        potentials = np.array(potentials).T
        potentials_df = pd.DataFrame(potentials,columns=cell_names)
        potentials_df.to_csv('potentials.csv', index=False)
        
        potentials_df.to_csv(f'cases_paper/figures_data/figure4b_{int(current*1000)}pA_potentials.csv', index=False)

if __name__ == '__main__':
        currents = [0.1,0.2,0.3,0.4,0.5,1]
        hoc_files = ['./components_homogenous/templates/feng.hoc','./components_homogenous/templates/SOM.hoc'] 
        # Load mechs
        #neuron.load_mechanisms(mechanisms_dir)
        h.nrn_load_dll('./components_homogenous/mechanisms/x86_64/libnrnmech.so')
        # Load the standard hoc file
        h.load_file('stdrun.hoc')
        for hoc_file in hoc_files:
            # Load the .hoc file
            h.load_file(hoc_file)
        
        for current in currents:
            main(['Cell_Af','Cell_Cf','InterneuronCellf','SOM_Cell','CR_Cell'],current)
