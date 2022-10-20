# generate_node_sets.py
import json
import random 

"""
Stimulation to PN, PV, and SOM via ACH
BF-PN: 89% 
BF-PV+: 7% 
BF-SOM+: 4% 
"""
output_file_name = "node_sets.json"

output_dict = {
  "biophysical_nodes": {
    "model_type": "biophysical"
  },
  "point_nodes": {
    "model_type": "point_process"
  },
  "vclamp_cells": [0,50,100,200,300,400,500,600,700,750]
}



scale = 1

n_pn = 800
n_pv = 93
n_som = 51

pn_percent = 0.89
pv_percent = 0.07
som_percent = 0.04

pn_ids = random.sample(range(0,n_pn),int(pn_percent*n_pn))
pn_ids.sort()
pv_ids = random.sample(range(n_pn,n_pn+n_pv),int(pv_percent*n_pv))
pv_ids.sort()
som_ids = random.sample(range(n_pn+n_pv,n_pn+n_pv+n_som),int(som_percent*n_som))
som_ids.sort()

print(f"Stimulating {len(pn_ids)} PN cells with ACH.")
print(f"Stimulating {len(pv_ids)} PV cells with ACH.")
print(f"Stimulating {len(som_ids)} SOM cells with ACH.")

output_dict["ach_pn"] = pn_ids
output_dict["ach_pv"] = pv_ids
output_dict["ach_som"] = som_ids

with open(output_file_name, "w") as outfile:
    json.dump(output_dict, outfile, indent=2)

print(f"File written to {output_file_name}")
