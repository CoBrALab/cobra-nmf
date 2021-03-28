#this script takes in nmf results and writes out .txt files
#two files are written out - one for each hemisphere listing the component scores at each vertex, and an additional column with a winner take all clustering
#so, if you have 10 comps, expect 11 columns
#use brainview and these outputs to visualize the spatial patterns

## LOAD MODULES/SOFTWARE
import os
import glob
import pandas as pd
import numpy as np

import sys
import pickle
import scipy
from scipy.io import savemat, loadmat
import scipy.stats
import argparse
parser=argparse.ArgumentParser(
    description='''This script reads in nmf results and outputs a .txt listing component scores and winnner take all labelling''')

group = parser.add_argument_group(title="Execution options")

group.add_argument(
    '--nmf_results', help='.mat file containing nmf results',required=True)
group.add_argument(
    '--mask_file', help='path to CIVET midline mask',required=True)
group.add_argument(
    '--file', help='path to CIVET midline mask',required=True)
group.add_argument(
    '--output_dir', help='directory to store outputs',required=True)

args=parser.parse_args()

#LOAD CIVET MASK TO IDENTIFY MIDLINE/ CORPOS COLLOSUM REGION
midline_mask = np.loadtxt(args.mask_file)

n_vertex=np.shape(midline_mask)[0]

if not os.path.exists(args.output_dir):
    os.mkdir(out_dir)
    
#load spatial components
W = scipy.io.loadmat(args.nmf_results)

compnum=np.shape(W['W'])[1]
left_W=W['W'][0:38561,:]
right_W=W['W'][38561:,:]



left_outarray = np.zeros((n_vertex,compnum+1))

for comp in range(0,compnum):
    valid_idx = 0
    for idx in range(0,n_vertex):
        if midline_mask[idx] == 1:
            left_outarray[idx,comp] = left_W[valid_idx,comp]
            valid_idx +=1

#define winner take all clustering
left_cluster = np.zeros((np.shape(left_W)[0],1))
for vertex in range(0,np.shape(left_W)[0]):
    left_cluster[vertex,0] = np.argmax(left_W[vertex,:])
left_cluster = left_cluster + 1 #plus one for zero based indexing

#now add clustering to out array, in the last column
valid_idx=0    
for idx in range(0,n_vertex):
    if midline_mask[idx] == 1:
        left_outarray[idx,-1] = left_cluster[valid_idx,0]
        valid_idx +=1
left_statmap = out_dir + "/left_k" + str(compnum) + ".txt"
np.savetxt(left_statmap,left_outarray)


right_outarray = np.zeros((n_vertex,compnum+1))
for comp in range(0,compnum):
    valid_idx = 0
    for idx in range(0,n_vertex)):
        if midline_mask[idx] == 1:
            right_outarray[idx,comp] = right_W[valid_idx,comp]
            valid_idx +=1

#define winner take all clustering
right_cluster = np.zeros((np.shape(right_W)[0],1))
for vertex in range(0,np.shape(right_W)[0]):
    right_cluster[vertex,0] = np.argmax(right_W[vertex,:])
right_cluster = right_cluster + 1 #plus one for zero based indexing

#now add clustering to out array, in the last column
valid_idx=0    
for idx in range(0,n_vertex):
    if midline_mask[idx] == 1:
        right_outarray[idx,-1] = right_cluster[valid_idx,0]
        valid_idx +=1
right_statmap = out_dir + "/right_k" + str(compnum) + ".txt"
np.savetxt(right_statmap,right_outarray)
