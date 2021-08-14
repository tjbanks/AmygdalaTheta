import pandas as pd
import h5py
import numpy as np
import sys

import warnings
warnings.filterwarnings("ignore") # Annoying h5py read warning

def run(df, scale=1, convergence=True):

    p_start = 0 * scale
    i_start = 900 * scale
    i_end = 1000 * scale

    total_p = i_start - p_start
    total_i = i_end - i_start
    

    p2p = df[(df.source < i_start) & (df.target < i_start)]
    p2i = df[(df.source < i_start) & (df.target >= i_start)]
    i2p = df[(df.source >= i_start) & (df.target < i_start)]
    i2i = df[(df.source >= i_start) & (df.target >= i_start)]

    cd = 'source'
    if convergence:
        cd = 'target'

    # if cd is 'source' we're counting the number of occurances for the source or divergence
    # if cd is 'target' we're counting the number of occurances for the target or convergence

    avg_p2p = np.average(list(p2p[cd].value_counts()))
    std_p2p = np.std(list(p2p[cd].value_counts()))

    avg_p2i = np.average(list(p2i[cd].value_counts()))
    std_p2i = np.std(list(p2i[cd].value_counts()))

    avg_i2p = np.average(list(i2p[cd].value_counts()))
    std_i2p = np.std(list(i2p[cd].value_counts()))

    avg_i2i = np.average(list(i2i[cd].value_counts()))
    std_i2i = np.std(list(i2i[cd].value_counts()))

    print()
    print("for P avg. #P received:\t" + str(round(avg_p2p,4)) + " (" + str(round(std_p2p,4)) + ")")
    print("for P avg. #I received:\t" + str(round(avg_i2p,4)) + " (" + str(round(std_i2p,4)) + ")")
    print("for I avg. #P received:\t" + str(round(avg_p2i,4)) + " (" + str(round(std_p2i,4)) + ")")
    print("for I avg. #I received:\t" + str(round(avg_i2i,4)) + " (" + str(round(std_i2i,4)) + ")")
    print()
    print("P2P Connectivity\t" + str(round(avg_p2p/total_p*100,2)) + "%")
    print("I2P Connectivity\t" + str(round(avg_i2p/total_i*100,2)) + "%")
    print("P2I Connectivity\t" + str(round(avg_p2i/total_p*100,2)) + "%")
    print("I2I Connectivity\t" + str(round(avg_i2i/total_i*100,2)) + "%")
    print()
    
    
    df1 = pd.DataFrame(columns=['source','target'])#flipped
    df1['source'] = df.target
    df1['target'] = df.source

    dfs = df.append(df1,ignore_index=True)
    df_rec = dfs[dfs.duplicated(keep='first')]
    
    p2p_rec = df_rec[(df_rec.source < i_start) & (df_rec.target < i_start)]
    #p2i_rec = df_rec[(df_rec.source < i_start) & (df_rec.target >= i_start)]
    p2i_rec = df_rec[((df_rec.source < i_start) & (df_rec.target >= i_start)) | ((df_rec.target < i_start) & (df_rec.source >= i_start))]
    i2i_rec = df_rec[(df_rec.source >= i_start) & (df_rec.target >= i_start)]
    
    avg_p2p_rec = np.average(list(p2p_rec[cd].value_counts()))/2
    avg_p2i_rec = np.average(list(p2i_rec[cd].value_counts()))/2
    avg_i2i_rec = np.average(list(i2i_rec[cd].value_counts()))/2

    print()
    print("reciprocal P2P\t" + str(round(avg_p2p_rec/total_p*100,2)) + "%")
    print("reciprocal P2I\t" + str(round(avg_p2i_rec/total_p*100,2)) + "%")
    print("reciprocal I2I\t" + str(round(avg_i2i_rec/total_i*100,2)) + "%")
    print()
    
    #import pdb;pdb.set_trace()

def bmtk_run(path):
    hdf = h5py.File(path)

    src = pd.DataFrame(hdf['edges']['BLA_to_BLA']['source_node_id'],columns=['source'])
    trg = pd.DataFrame(hdf['edges']['BLA_to_BLA']['target_node_id'],columns=['target'])
    df = pd.concat([src,trg],axis=1)
    run(df)

def feng_run(path,scale=1):
    
    #Iterate line by line in file
    #create a dataframe with source/target columns
    #each line is the source and each item in the line is the target, append to the df
    #call run(df)

    f = open(path, 'r')
    lines = f.readlines()
    df = None

    row = 0
    # Strips the newline character
    for line in lines:
        #print("Line{}: {}".format(row, line.strip()))
        sources = line.strip().split('\t')
        dft = pd.DataFrame(sources,columns=['source'])
        dft.source = dft.source.astype(int)
        dft['target'] = row
        if df is None:
            df = dft
        else:
            df = pd.concat([df,dft],ignore_index=True)
        row += 1
        
    run(df,scale=scale)
    #import pdb;pdb.set_trace()

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        #bmtk_run(sys.argv[-1])
        feng_run('./source_material/active_syn_op_new')
        #feng_run('./source_material/FengEtAl2019/input/active_syn_op',scale=27)
    else:
        bmtk_run('./network/BLA_BLA_edges.h5')
