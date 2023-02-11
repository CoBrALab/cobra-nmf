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

