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
import pickle
import scipy
from scipy.io import savemat, loadmat
from scipy import stats

input_list=["wholebrain_ct_raw","wholebrain_sa_raw","wholebrain_curv_raw"] #MODIFY to match the filenames of the wholebrain files created in your sample_extract*py scripts. should be name of the file without the .mat extension

z_dict = {}

for metric in input_list:
    fname = metric + ".mat"
    res = loadmat(fname) #load raw data
    x_z = np.asarray(stats.zscore(res['X'],axis=None)) #zscore, across both subjects and vertices
    z_dict[metric] = x_z


#concatenate each zscored shifted matrix together
#forms vertex X subject*n_metrics matrix
metric=input_list[0]
wb_z_all = z_dict[metric]
#if only one metric, skip the for loop below, go straight to line 38 #
for metric in input_list[1:]:
    print(metric)
    wb_z_all = np.concatenate((wb_z_all, z_dict[metric]),axis=1)

wb_z_s_all = wb_z_all - np.min(wb_z_all)

#write out z scored, shifted data fow whole group nmf analysis
fname = sys.argv[1]
print(fname, np.shape(wb_z_s_all), np.min(wb_z_s_all))
savemat(fname, {'X': wb_z_s_all})
