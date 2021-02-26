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
import argparse
import hdf5storage
options = hdf5storage.Options(oned_as = 'column', matlab_compatible = True, action_for_matlab_incompatible = 'error')


parser=argparse.ArgumentParser(
    description='''This script concatenates voxel x subject matrices into one voxel x subject*num_metrics matrix.
    Each metric matrix is z scored prior to concatenation and final matrix is shifted by min value to obtain non-negativity ''')

parser.add_argument(
    "--inputs",help="metric matrices to concatenate", metavar='list', nargs='+', required=True)

parser.add_argument(
    "--output", help='output .mat filename', default='output.mat')

args=parser.parse_args()


#shift

#initiate matrix with first input
res = hdf5storage.loadmat(args.inputs[0])
z_all = np.asarray(stats.zscore(res['X'],axis=None))

#zscore and concatenate remaining inputs
num_metrics = len(args.inputs)
for m in range(1,num_metrics):
    res = hdf5storage.loadmat(args.inputs[m])
    z_all = np.concatenate((z_all, np.asarray(stats.zscore(res['X'],axis=None))), axis = 1)

#shift by min value
z_all_shift = z_all - np.min(z_all)
    
#write out z scored, shifted data fow whole group nmf analysis
fname = args.output
print('saving',fname, 'with shape',np.shape(z_all_shift), 'and minimum val',np.min(z_all_shift))
hdf5storage.savemat(fname, {'X': z_all_shift}, format='7.3', options=options)

#input_list=["wholebrain_t1t2","wholebrain_dbm"] #MODIFY to match the filenames of the wholebrain files created in your sample_extract*py scripts. should be name of the file without the .mat extension

#z_dict = {}

#for file in input_list:
#    fname = file + ".mat"
#    res = hdf5storage.loadmat(fname) #load raw data
#    x_z = np.asarray(stats.zscore(res['X'],axis=None)) #zscore, across both subjects and vertices
#    #x_z_s = x_z - np.min(x_z) #shift so all positive
#    z_dict[file] = x_z


#concatenate each zscored shifted matrix together
#forms vertex X subject*n_metrics matrix
#file=input_list[0]
#print(file)
#wb_z_all = z_dict[file]
#for file in input_list[1:]:
#    print(file)
#    wb_z_all = np.concatenate((wb_z_all, z_dict[file]),axis=1)

#wb_z_s_all = wb_z_all - np.min(wb_z_all)

#write out z scored, shifted data fow whole group nmf analysis
#fname = sys.argv[1]
#print(fname, np.shape(wb_z_s_all), np.min(wb_z_s_all))
#hdf5storage.savemat(fname, {'X': wb_z_s_all}, format='7.3', options=options)
