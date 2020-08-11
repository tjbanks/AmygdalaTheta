from pathlib import Path
import json
import pandas as pd

def convert_json():

    txt_folder = Path('../components/synaptic_models').rglob('*.json')
    files = [x for x in syn_folder]

    for f in files:
        json_synapse = json.read(f)
        df = pd.DataFrame()

