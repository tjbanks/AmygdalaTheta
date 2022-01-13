import pandas as pd
import h5py
import numpy as np
import sys

import warnings
warnings.filterwarnings("ignore") # Annoying h5py read warning

def run(df, scale=1, convergence=True):

    p_start = 0 * scale
    i_start = 800 * scale #93
    s_start = 893 * scale #51
    c_start = 944 * scale #56
    c_end = 1000 * scale

    total_p = i_start - p_start
    total_i = s_start - i_start
    total_s = c_start - s_start
    total_c = c_end - c_start

    p2p = df[(df.source < i_start) & (df.target < i_start)]
    p2i = df[(df.source < i_start) & (df.target >= i_start) & (df.target < s_start)]
    p2s = df[(df.source < i_start) & (df.target >= s_start) & (df.target < c_start)]
    p2c = df[(df.source < i_start) & (df.target >= c_start)]

    i2p = df[(df.source >= i_start) & (df.source < s_start) & (df.target < i_start)]
    i2i = df[(df.source >= i_start) & (df.source < s_start) & (df.target >= i_start) & (df.target < s_start)]
    i2s = df[(df.source >= i_start) & (df.source < s_start) & (df.target >= s_start) & (df.target < c_start)]
    i2c = df[(df.source >= i_start) & (df.source < s_start) & (df.target >= c_start)]
 
    s2p = df[(df.source >= s_start) & (df.source < c_start) & (df.target < i_start)] 
    s2i = df[(df.source >= s_start) & (df.source < c_start) & (df.target >= i_start) & (df.target < s_start)]
    s2s = df[(df.source >= s_start) & (df.source < c_start) & (df.target >= s_start) & (df.target < c_start)]
    s2c = df[(df.source >= s_start) & (df.source < c_start) & (df.target >= c_start)]

    c2p = df[(df.source >= c_start) & (df.target < i_start)]
    c2i = df[(df.source >= c_start) & (df.target >= i_start) & (df.target < s_start)]
    c2s = df[(df.source >= c_start) & (df.target >= s_start) & (df.target < c_start)]
    c2c = df[(df.source >= c_start) & (df.target >= c_start)]

    cd = 'source'
    if convergence:
        cd = 'target'

    # if cd is 'source' we're counting the number of occurances for the source or divergence
    # if cd is 'target' we're counting the number of occurances for the target or convergence

    avg_p2p = np.average(list(p2p[cd].value_counts()))
    std_p2p = np.std(list(p2p[cd].value_counts()))

    avg_p2i = np.average(list(p2i[cd].value_counts()))
    std_p2i = np.std(list(p2i[cd].value_counts()))

    avg_p2s = np.average(list(p2s[cd].value_counts()))
    std_p2s = np.std(list(p2s[cd].value_counts()))

    avg_p2c = np.average(list(p2c[cd].value_counts()))
    std_p2c = np.std(list(p2c[cd].value_counts()))


    avg_i2p = np.average(list(i2p[cd].value_counts()))
    std_i2p = np.std(list(i2p[cd].value_counts()))

    avg_i2i = np.average(list(i2i[cd].value_counts()))
    std_i2i = np.std(list(i2i[cd].value_counts()))

    avg_i2s = np.average(list(i2s[cd].value_counts()))
    std_i2s = np.std(list(i2s[cd].value_counts()))

    avg_i2c = np.average(list(i2c[cd].value_counts()))
    std_i2c = np.std(list(i2c[cd].value_counts()))

    
    avg_s2p = np.average(list(s2p[cd].value_counts()))
    std_s2p = np.std(list(s2p[cd].value_counts()))

    avg_s2i = np.average(list(s2i[cd].value_counts()))
    std_s2i = np.std(list(s2i[cd].value_counts()))

    avg_s2s = np.average(list(s2s[cd].value_counts()))
    std_s2s = np.std(list(s2s[cd].value_counts()))

    avg_s2c = np.average(list(s2c[cd].value_counts()))
    std_s2c = np.std(list(s2c[cd].value_counts()))
    

    avg_c2p = np.average(list(c2p[cd].value_counts()))
    std_c2p = np.std(list(c2p[cd].value_counts()))

    avg_c2i = np.average(list(c2i[cd].value_counts()))
    std_c2i = np.std(list(c2i[cd].value_counts()))

    avg_c2s = np.average(list(c2s[cd].value_counts()))
    std_c2s = np.std(list(c2s[cd].value_counts()))

    avg_c2c = np.average(list(c2c[cd].value_counts()))
    std_c2c = np.std(list(c2c[cd].value_counts()))

    print()
    print("for P avg. #P received:\t" + str(round(avg_p2p,4)) + " (" + str(round(std_p2p,4)) + ")")
    print("for P avg. #I received:\t" + str(round(avg_i2p,4)) + " (" + str(round(std_i2p,4)) + ")")
    print("for P avg. #S received:\t" + str(round(avg_s2p,4)) + " (" + str(round(std_p2p,4)) + ")")
    print("for P avg. #C received:\t" + str(round(avg_c2p,4)) + " (" + str(round(std_i2p,4)) + ")")
    print()
    print("for I avg. #P received:\t" + str(round(avg_p2i,4)) + " (" + str(round(std_p2i,4)) + ")")
    print("for I avg. #I received:\t" + str(round(avg_i2i,4)) + " (" + str(round(std_i2i,4)) + ")")
    #print("for I avg. #S received:\t" + str(round(avg_s2i,4)) + " (" + str(round(std_s2i,4)) + ")")
    print("for I avg. #C received:\t" + str(round(avg_c2i,4)) + " (" + str(round(std_c2i,4)) + ")")
    print()
    print("for S avg. #P received:\t" + str(round(avg_p2s,4)) + " (" + str(round(std_p2s,4)) + ")")
    print("for S avg. #I received:\t" + str(round(avg_i2s,4)) + " (" + str(round(std_i2s,4)) + ")")
    #print("for S avg. #S received:\t" + str(round(avg_s2s,4)) + " (" + str(round(std_s2s,4)) + ")")
    print("for S avg. #C received:\t" + str(round(avg_c2s,4)) + " (" + str(round(std_c2s,4)) + ")")
    print()
    print("for C avg. #P received:\t" + str(round(avg_p2c,4)) + " (" + str(round(std_p2c,4)) + ")")
    #print("for C avg. #I received:\t" + str(round(avg_i2c,4)) + " (" + str(round(std_i2c,4)) + ")")
    #print("for C avg. #S received:\t" + str(round(avg_s2c,4)) + " (" + str(round(std_s2c,4)) + ")")
    #print("for C avg. #C received:\t" + str(round(avg_c2c,4)) + " (" + str(round(std_c2c,4)) + ")")
    print()
    print("P2P Connectivity\t" + str(round(avg_p2p/total_p*100,2)) + "%")
    print("I2P Connectivity\t" + str(round(avg_i2p/total_i*100,2)) + "%")
    print("S2P Connectivity\t" + str(round(avg_s2p/total_s*100,2)) + "%")
    print("C2P Connectivity\t" + str(round(avg_c2p/total_c*100,2)) + "%")
    print()
    print("P2I Connectivity\t" + str(round(avg_p2i/total_p*100,2)) + "%")
    print("I2I Connectivity\t" + str(round(avg_i2i/total_i*100,2)) + "%")
    #print("S2I Connectivity\t" + str(round(avg_s2i/total_s*100,2)) + "%")
    print("C2I Connectivity\t" + str(round(avg_c2i/total_c*100,2)) + "%")
    print()
    print("P2S Connectivity\t" + str(round(avg_p2s/total_p*100,2)) + "%")
    print("I2S Connectivity\t" + str(round(avg_i2s/total_i*100,2)) + "%")
    #print("S2S Connectivity\t" + str(round(avg_s2s/total_s*100,2)) + "%")
    print("C2S Connectivity\t" + str(round(avg_c2s/total_c*100,2)) + "%")
    print()
    print("P2C Connectivity\t" + str(round(avg_p2c/total_p*100,2)) + "%")
    #print("I2C Connectivity\t" + str(round(avg_i2c/total_i*100,2)) + "%")
    #print("S2C Connectivity\t" + str(round(avg_s2c/total_s*100,2)) + "%")
    #print("C2C Connectivity\t" + str(round(avg_c2c/total_c*100,2)) + "%")
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
    if 'homogenous' in sys.argv:
        bmtk_run('./network_homogenous/BLA_BLA_edges.h5')
    else:
        bmtk_run('./network/BLA_BLA_edges.h5')
