# Amygdala Theta
#### Tyler Banks and Pete Canfield
Modeling Basal Forebrain GABAergic Neuromodulation of the Amygdala Theta Rhythm

## Connectivity Matrix
* \+ : Excitatory
* \- : Inhibitory

| Source/Target  | Pyr | PV+ | SOM+ | CR+ |
|----------------|-----|-----|------|-----|
| Pyr            | +   | +   | +    | +   |
| PV+            | -   |     | -    |     |
| SOM+           | -   |     |      |     |
| CR+            | -   | -   | -    |     |
| Thalamic Input | +   |     | +    | +   |
| VP/SI Input    | +   | +   |      |     |

Current model has no connection between PV to SOM but does have connectivity between PV and PV.


## Running the Model

```
# Generate inputs
python manage input generate all

# Build the model
python build_network.py

# Run the model
python run_bionet.py
```

## Single Cell Profiling

Cell templates can be found in:
```
./components/templates/templates.hoc
```
Templates used are:
```
Cell_A
Cell_C
PVCell
SOM_Cell
CR_Cell
```

### Tuning

Individual tuning interfaces can be run with the following:
```
python manage tune
```


### FI Curves

FI Curves can be plotted using:

```
python manage tune fi

# Pyramidal 
python manage tune fi A
python manage tune fi C

# PV+
python manage tune fi PV

# CR+
python manage tune fi CR

# SOM+
python manage tune fi SOM

bmtool util cell --template SOM_Cell vhseg --gleak gl_ichan2OLM --eleak el_ichan2OLM --segvars eh,enat,ekf,ek,el_ichan2OLM,elca,gnatbar_ichan2OLM,gkfbar_ichan2OLM,eh,gkhbar_ihOLM,gkAbar_AOLM,catau_ccanlOLM,gsAHbar_sAHPOLM,glcabar_lcaOLM,gcatbar_catOLM,gbar_napOLM --othersec dend[0],dend[1] --fminpa 44 --fmaxpa 305 --fincrement 20

```

### Current Clamps

FI Curves can be plotted using:

```
python manage tune clamp

# Pyramidal 
python manage tune clamp A
python manage tune clamp C

# PV+
python manage tune clamp PV

# CR+
python manage tune clamp CR

# SOM+
python manage tune clamp SOM
```

### Syanpses

FI Curves can be plotted using:

```
python manage tune syn

# Exp2Syn Pyramidal Cell A
python manage tune syn exp2syn_Cell_A


```

## Network Profiling

### Plot Connection Totals
```
bmtool plot connection total
```

## Appendix

### Installing `manage` dependencies
`manage` uses click for an easy to use command line interface. Install the requirements for this using:
```
pip install -r requirements.txt
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

Built using [BMTK](https://github.com/AllenInstitute/bmtk) checkout 52fee.

```
git clone https://github.com/AllenInstitute/bmtk
cd bmtk
git checkout 52fee3b230ceb14a666c46f57f2031c38f1ac5b1 .
python setup.py develop
````
