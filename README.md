# Amygdala Theta
#### Tyler Banks, Pete Canfield, and Matthew Stroud
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
| `vpsi_inh_spikes_nonrhythmic.h5`         | VPSI poisson input for non-rhythmic inhibition | [`build_input.py`](./build_input.py)|
| **`vpsi_inh_spikes.h5`**         | VPSI 8Hz Rhythmic 3Hz poisson input | [`matlab/generatethetainputs.m`](./matlab/generatethetainputs.m) & [`matlab/convert_spikematrix.py`](matlab/convert_spikematrix.py)|

(primary files **bold**)


```
python generate_input.py

# To regenerate rhythmic inputs (not usually needed)
matlab &
generatethetainputs

python matlab/convert_spikematrix.py
```

### 2. Build network configuration files

```
python build_network.py
```

### 3. Execute run script

The network can be tested using any one of the simulation configuration files listed below (`python run_bionet.py [configuration file]`)

| Configuration file | Details |
|--------------------|---------|
| [simulation_config.json](./simulation_config.json) | 
| [simulation_configECP.json](./simulation_configECP.json) | 
| [simulation_configECP_base.json](./simulation_configECP_base.json) | 
| [simulation_configECP_base_vclamp.json](./simulation_configECP_base_vclamp.json) | 
| [simulation_configECP_gamma.json](./simulation_configECP_gamma.json) | 
| [simulation_configECP_vpsi.json](./simulation_configECP_vpsi.json) | 
| [simulation_configECP_vpsi_vclamp.json](./simulation_configECP_vpsi_vclamp.json) | 
| [simulation_configECP_vpsi_vclamp_nonrhythmic.json](./simulation_configECP_vpsi_vclamp_nonrhythmic.json) | 

```
python run_bionet.py
```

### Analyze the model

Analysis of the model is primarily comprised of a spike raster, mean firing rates, raw LFP, and LFP PSD.
```
matlab &
analysis('../outputECP/ecp.h5','../outputECP/spikes.h5');
```




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
