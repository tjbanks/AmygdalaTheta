import pandas as pd
import pprint

from bmtool.util.util import relation_matrix

from synapses import syn_params_dicts

dparams = {}
rm = {}

def relation(**kwargs):
    edges = kwargs["edges"]
    source_id_type = kwargs["sid"]
    target_id_type = kwargs["tid"]
    source_id = kwargs["source_id"]
    target_id = kwargs["target_id"]
    t_list = kwargs["target_nodes"]
    s_list = kwargs["source_nodes"]
    
    cons = edges[(edges[source_id_type] == source_id) & (edges[target_id_type]==target_id)]
    
    #mean convergence
    vc = t_list.apply(pd.Series.value_counts)
    vc = vc[target_id_type].dropna().sort_index()
    count = vc.loc[target_id]
    total = edges[(edges[source_id_type] == source_id) & (edges[target_id_type]==target_id)].count()
    total = total.source_node_id # may not be the best way to pick
    mean = round(total/count,1)

    #std
    std = cons.apply(pd.Series.value_counts).target_node_id.dropna().std()

    #weights
    params = cons.dynamics_params.unique()
    if len(params) == 0:
        info_matrix = {}     
    else:
        param_json = params[0]

        dp = dparams[param_json]
        #'initW_lognormal_mean': '2', 'initW_lognormal_std': '1', 'bACH': '0.0'    
        weight_mean = float(dp.get('initW_lognormal_mean',0))
        weight_std = float(dp.get('initW_lognormal_std',0))
        weight_ach = float(dp.get('bACH',0))
        info_matrix = {'convergence_mean':round(mean,2),
            'convergence_std':round(std,2),
            'weight_mean':weight_mean,
            'weight_std':weight_std,
            'weight_ach':weight_ach}

    if not rm.get(source_id):
        rm[source_id] = {}

    rm[source_id][target_id] = info_matrix


def main(config):

    nodes = None
    edges = None
    sources = ['BLA']
    targets = ['BLA']
    sids = ['pop_name']
    tids = ['pop_name']
    prepend_pop = True
    
    # Relation matrix is kind of my secret formula for loading up all the edge files and sticking
    # them into a nice and easy to read pandas dataframe. The 'relation_func' is what that dataframe
    # gets sent to and called for each source/dest combination
    # In this case for each shell cell type, and every BLA cell call add_inputs
    relation_matrix(config,nodes,edges,sources,targets,sids,tids,prepend_pop,relation_func=relation)

    #t = PrettyTable()
    #t.field_names = rm.keys()
    columns = ['S/T']+list(rm.keys()) #set the initial values to the targets for top row
    output = []
    for i, (source, targets) in enumerate(rm.items()):
        row = [source]
        for j, (target, properties) in enumerate(targets.items()):
            infostr = f"{properties.get('convergence_mean',0.0)} ({properties.get('convergence_std',0.0)})\n{properties.get('weight_mean',0.0)} ({properties.get('weight_std',0.0)})"
            row.append(infostr)
        output.append(row)

    df = pd.DataFrame(output,columns=columns)
    df.to_csv('connections.csv',index=False)
    print(df)

if __name__ == '__main__':
    dparams = syn_params_dicts('./components_homogenous/synaptic_models')
    main('simulation_configECP_base_homogenous.json')
