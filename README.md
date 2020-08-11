# Amygdala Theta
#### Tyler Banks and Pete Canfield
Modeling Basal Forebrain GABAergic Neuromodulation of the Amygdala Theta Rhythm

## Connectivity Matrix
* \+ : Excitatory
* \- : Inhibitory

| Source/Target  | Pyr | PV+ | SOM+ | CR+ |
|----------------|-----|-----|------|-----|
| Pyr            | +   | +   | +    | +   |
| PV+            | -   |     |      |     |
| SOM+           | -   |     |      |     |
| CR+            | -   |     |      |     |
| Thalamic Input | +   |     | +    | +   |
| VP/SI Input    | +   | +   |      |     |

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
python manage tune fi Cell_A
python manage tune fi Cell_C

# PV+
python manage tune fi PV

# CR+
python manage tune fi CR

# SOM+
python manage tune fi SOM

```

### Current Clamps

FI Curves can be plotted using:

```
python manage tune clamp

# Pyramidal 
python manage tune clamp Cell_A
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