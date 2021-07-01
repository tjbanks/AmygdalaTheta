# Amygdala Theta
#### Code by Tyler Banks, Pete Canfield, and Matthew Stroud. In partnership with Unal Lab (Tuna and Unal)
Modeling Basal Forebrain GABAergic Neuromodulation of the Amygdala Theta Rhythm


## Running the Model

### 1. Build network input files :

| Input file                    | Purpose | Generating code |
|-------------------------------|---------|-----------------|
| **`thalamus_pyr_spikes.h5`**         | Thalamic 2Hz poisson input to Pyramidal A and C cells | [`build_input.py`](./build_input.py)|
| `thalamus_pv_spikes.h5`         | Thalamic 2Hz poisson input to Interneurons (unused - TESTING ONLY) | [`build_input.py`](./build_input.py)|
| **`thalamus_som_spikes.h5`**         | Thalamic 2Hz poisson input to SOM cells| [`build_input.py`](./build_input.py)|
| **`thalamus_cr_spikes.h5`**         | Thalamic 2Hz poisson input CR+ cells| [`build_input.py`](./build_input.py)|
| `vpsi_pyr_spikes.h5`         | VPSI 2Hz poisson input to Pyramidal A and C cells (unused - TESTING ONLY) | [`build_input.py`](./build_input.py)|
| `vpsi_pv_spikes.h5`         | VPSI 2Hz poisson input to Interneurons (unused - TESTING ONLY) | [`build_input.py`](./build_input.py)|
| `vpsi_inh_spikes_nonrhythmic.h5`         | VPSI poisson input to PN A and C for non-rhythmic inhibition | [`build_input.py`](./build_input.py)|
| **`vpsi_inh_spikes.h5`**         | VPSI 8Hz Rhythmic 3Hz poisson input to PN A and C | [`matlab/generatethetainputs.m`](./matlab/generatethetainputs.m) & [`matlab/convert_spikesmatrix.py`](matlab/convert_spikesmatrix.py)|

(primary files **bold**)

To generate the above files first, create Thalamic and VSPI inputs using
```
python generate_input.py

```

To generate 8Hz rhythmic inputs (based on Fink et al 2015)
```
matlab &
generatethetainputs
```
The files will then need to be converted to the Sonata format BMTK understands.
```
python matlab/convert_spikematrix.py
```

Once these steps have been completed, input files **WILL NOT** need to be regenerated again.

### 2. Build network configuration files

Building of the network can be customized by altering `build_network.py` configuration at the beginning of the script.

To generate necessary network specification files in the `[network](./network)` directory, execute

```
python build_network.py
```

### 3. Execute run script

The network can be tested using any one of the simulation configuration files listed below (`python run_bionet.py [configuration file]`)

| Configuration file | Input Details |Notes|
|--------------------|---------------|-----|
| [simulation_config.json](./simulation_config.json) | Thalamic 2Hz to PN, INT, SOM, CR, VPSI 2Hz exc and 8Hz inh to PN/INT | **deprecated** example |
| [simulation_configECP.json](./simulation_configECP.json) | Thalamic 2Hz to PN, SOM, CR, VPSI 2Hz exc and 8Hz inh to PN/INT | **deprecated** example with LFP recording electrodes [linear_electrode.csv](components/recXelectrodes/linear_electrode.csv)|
| **[simulation_configECP_base.json](./simulation_configECP_base.json)** | Thalamic 2Hz to PN, SOM, CR | Quiet-waking state, baseline configuration file |
| [simulation_configECP_base_vclamp.json](./simulation_configECP_base_vclamp.json) | Thalamic 2Hz to PN, SOM, CR | Same as [simulation_configECP_base.json](./simulation_configECP_base.json), 11 voltage clamped PN cells, igaba from SOM+ synapses are recorded and placed in `outputECP/syn_report.h5` - Analysis completed using [ipsc_analysis.m](./matlab/ipsc_analysis.m) notes below |
| [simulation_configECP_gamma.json](./simulation_configECP_gamma.json) | | Gamma testing for original model replicated from Feng et al 2019|
| **[simulation_configECP_vpsi.json](./simulation_configECP_vpsi.json)** | Thalamic 2Hz to PN, SOM, CR, VPSI 8Hz rhythmic inh to PN/INT | Primary VPSP input test|
| [simulation_configECP_vpsi_vclamp.json](./simulation_configECP_vpsi_vclamp.json) | Thalamic 2Hz to PN, SOM, CR, VPSI 8Hz rhythmic inh to PN/INT | Same as [simulation_configECP_vpsi.json](./simulation_configECP_vpsi.json), 11 voltage clamped PN cells, igaba from SOM+ synapses are recorded and placed in `outputECP/syn_report.h5` - Analysis completed using [ipsc_analysis.m](./matlab/ipsc_analysis.m) notes below
| [simulation_configECP_vpsi_vclamp_nonrhythmic.json](./simulation_configECP_vpsi_vclamp_nonrhythmic.json) | Thalamic 2Hz to PN, SOM, CR, VPSI 3Hz Poisson inh to PN/INT | Testing Non-rhythmic inhibition to PN/PV|

Primary tests **bold**

#### Single Core Mode
1000 Cell models typically run for **4-6 hours** on a single core.

```
python run_bionet.py simulation_configECP_base.json
```

#### Parallel Mode
1000 Cell models run for **5-6 minutes** on ~50 cores.
```
mpirun -n 50 nrniv -mpi -python run_network.py simulation_configECP_base.json
```

### Analysis of the model

Analysis of the model is primarily comprised of a spike raster, mean firing rates, raw LFP, and LFP PSD in **MATLAB**. Launch MATLAB using:
```
matlab &
```

#### [analysis.m](./matlab/analysis.m)

Plots a spike raster, mean firing rates, LFP and LFP PSD on 4 separate graphs.
```
analysis('../outputECP/ecp.h5','../outputECP/spikes.h5');
```


#### [ipsc_analysis.m](./matlab/ipsc_analysis.m)

Used in conjuction with [simulation_configECP_base_vclamp.json](./simulation_configECP_base_vclamp.json) and [simulation_configECP_vpsi_vclamp.json](./simulation_configECP_vpsi_vclamp.json) to sum igaba currents and perform a PSD to produce the raw signal and PSD plots.
```
ipsc_analysis('../outputECP/syn_report.h5';
```

#### [connection_info.py](./connection_info.py)

Used to print the connectivity between cell types.

```
> python connection_info.py
        total uni   bi
PN->PN	2.0	  1.96	0.04
PN->PV	26.82	11.24	15.58
PN->SOM	31.19	29.17	2.01
PN->CR	18.43	16.41	2.02
PV->PN	52.0	36.42	15.58
PV->PV	22.92	17.41	5.5
PV->SOM	9.8	  9.8	  0.0
PV->CR	0.0	  0.0	  0.0
SOM->PN	6.57	4.55	2.01
SOM->PV	0.0	  0.0	  0.0
SOM->SOM	0.0	0.0	  0.0
SOM->CR	0.0   0.0   0.0
CR->PN	11.59	9.57	2.02
CR->PV	29.7	29.7	0.0
CR->SOM	75.25	75.25	0.0
CR->CR	0.0	  0.0	  0.0
```

### Other important files

| File/Directory | Description |
|------|-------------|
|[synapses.py](./synapses.py)| When adding new synapse types to be used by `build_network.py` they must be defined here |
|[components/templates/feng.hoc](./components/templates/feng.hoc)| Contains the PN A, PN C and INT hoc cell templates|
|[components/templates/SOM.hoc](./components/templates/SOM.hoc)| Contains the SOM and CR hoc cell templates|
|[components/synaptic_models](./components/synaptic_models)| Contains parameter configuration files for synapses used |
|[components/mechanisms/modfiles](./components/mechanisms/modfiles)| `.mod` definition files for cells and synapses|
|[tuning/current_clamp](./tuning/current_clamp)| Directory of easy to use hoc-based testers/tuners to simulate current injection into cells in the model



## Connectivity Matrix
* \+ : Excitatory
* \- : Inhibitory
* \+/\- : Both

| Source/Target  | Pyr | PV+ | SOM+ | CR+ |
|----------------|-----|-----|------|-----|
| Pyr            | +   | +   | +    | +   |
| PV+            | -   | -   | -    |     |
| SOM+           | -   |     |      |     |
| CR+            | -   | -   | -    |     |
| Thalamic Input | +   |     | +    | +   |
| VP/SI Input    | -   | -   |      |     |

Note: current model has no connection between PV to SOM.

## Single Cell Profiling

Cell templates can be found in:
```
./components/templates/templates.hoc
./components/templates/SOM.hoc
```
Templates used are:
```
Cell_A
Cell_C
PVCell
SOM_Cell
CR_Cell
```



## Network Profiling

### Plot Connection Totals
```
bmtool plot connection total
```


### BMTOOLS

*Version 0.2.1+*

Install BMTools by running

```
pip install bmtool
```

Help Examples:
```
bmtool --help
bmtool util --help
bmtool util cell --help
bmtool util cell fi --help
```

### BMTK

Built using [BMTK](https://github.com/AllenInstitute/bmtk) checkout 52fee and later. Parallel netcon recording was added 6/29/21.

```
git clone https://github.com/AllenInstitute/bmtk
cd bmtk
python setup.py develop
````


## Appendix

### Installing NEURON and BMTK

See [https://github.com/tjbanks/easy-nrn-install](https://github.com/tjbanks/easy-nrn-install)
