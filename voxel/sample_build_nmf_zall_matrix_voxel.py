#this python script loads in raw data for each metric 
#zscores data for each metric across BOTH subjects and vertices
#ie treats matrix of ct data as one array, zscore across all, repeat for each metric
#then concatenates and writes out matrices for nmf usage

## LOAD MODULES/SOFTWARE
import os
import glob
import pandas as pd
import numpy as np

import sys
import scipy
from scipy import stats
import hdf5storage
options = hdf5storage.Options(oned_as = 'column', matlab_compatible = True, action_for_matlab_incompatible = 'error')

input_list=["wholebrain_t1t2","wholebrain_dbm"] #MODIFY to match the filenames of the wholebrain files created in your sample_extract*py scripts. should be name of the file without the .mat extension

z_dict = {}

for file in input_list:
    fname = file + ".mat"
    res = hdf5storage.loadmat(fname) #load raw data
    x_z = np.asarray(stats.zscore(res['X'],axis=None)) #zscore, across both subjects and vertices
    #x_z_s = x_z - np.min(x_z) #shift so all positive
    z_dict[file] = x_z


#concatenate each zscored shifted matrix together
#forms vertex X subject*n_metrics matrix
file=input_list[0]
#print(file)
wb_z_all = z_dict[file]
for file in input_list[1:]:
    print(file)
    wb_z_all = np.concatenate((wb_z_all, z_dict[file]),axis=1)

wb_z_s_all = wb_z_all - np.min(wb_z_all)

#write out z scored, shifted data fow whole group nmf analysis
fname = sys.argv[1]
print(fname, np.shape(wb_z_s_all), np.min(wb_z_s_all))
hdf5storage.savemat(fname, {'X': wb_z_s_all}, format='7.3', options=options)
