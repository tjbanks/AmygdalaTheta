# AmygdalaTheta
Modeling Basal Forebrain GABAergic Neuromodulation of the Amygdala Theta Rhythm

## Running the Model

```
# Generate inputs
python manage.py input generate all

# Build the model
python build_network.py

# Run the model
python run_bionet.py
```


## Appendix

### Installing `manage.py` dependencies
`manage.py` uses click for an easy to use command line interface. Install the requirements for this using:
```
pip install -r requirements.txt
```

### BMTK Version

Built using [BMTK](https://github.com/AllenInstitute/bmtk) checkout 52fee.

```
git clone https://github.com/AllenInstitute/bmtk
cd bmtk
git checkout 52fee3b230ceb14a666c46f57f2031c38f1ac5b1 .
python setup.py develop
````