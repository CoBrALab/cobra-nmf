#this script loads in raw .txt files from each subject and:
#1) concatenates each file to build a vertex X subj matrix for each metric, for both left and right hemisphers
#2) concatenates the left and right hemisphere data to build matrix of vertex X subj for whole brain 
#3) write out .mat files containing left, right, whole brain data

## LOAD MODULES/SOFTWARE
import os
import glob
import pandas as pd
import numpy as np

import sys
import pickle
import scipy
from scipy.io import savemat, loadmat
import argparse

parser=argparse.ArgumentParser(
    description='''This script extracts vertex data from .txt files and outputs a vertex x subject matrix
    in .mat format''')

parser.add_argument('--metric',type=str, nargs='+', action='append')
parser.add_argument('--metric_column',type=str, nargs='+', action='append')
group = parser.add_argument_group(title="Execution options")

group.add_argument(
    '--input_csv', help='demographic spreadsheet, must contain subject id',required=True)
group.add_argument(
    '--mask_file', help='path to CIVET midline mask',required=True)

args=parser.parse_args()

def load_vertex_data(df, column, n_subjects, n_vertex, mask=None):
#helper function to cycle through array of filepaths, load data, build matrix
    filepaths = df[[column]].values
    for f in range(0,n_subjects):
        fname=filepaths[f,0]
        #print(fname)
        if f == 0:
            vertex_data = np.loadtxt(fname).reshape(1,n_vertex)
        else:
            vertex_data = np.concatenate(
            (vertex_data, np.loadtxt(fname).reshape(1,n_vertex)), axis=0)

    vertex_mean = np.mean(vertex_data,axis=0)
    vertex_std = np.std(vertex_data,axis=0)

    if mask is not None:
        return vertex_data[:,mask], vertex_mean, vertex_std
    else:
        return vertex_data, vertex_mean, vertex_std

def save_mat(x,key,fname):
    print("Saving ", np.shape(x), key, "to", fname)
    scipy.io.savemat(fname, {'X': x})

#LOAD CIVET MASK TO IDENTIFY MIDLINE/ CORPOS COLLOSUM REGION
#IDENTIFY 'VALID VERTICES' - IE VERTICES NOT IN THIS REGION
midline_mask = np.loadtxt(args.mask_file)
valid_vertices = np.where(midline_mask==1)[0].tolist() # list of valid indices in civet .txt file
invalid_vertices = np.where(midline_mask==0)[0].tolist() #list of invalid indices in civet .txt file

df_inputs = pd.read_csv(args.input_csv)
n_subjects=df_inputs.shape[0] #num rows in spreadsheet
n_vertex=np.shape(midline_mask)[0]

metric_dict = {}
for m_idx,m in enumerate(args.metric):
    metric = m[0] #args.metric is list of lists, each w length 1. we just want the metric, not whole list
    print('extracting', metric)
    for c_idx,c in enumerate(args.metric_column[m_idx]):
        vertex_data, vertex_mean, vertex_std = load_vertex_data(df_inputs,c,n_subjects,n_vertex, valid_vertices)

        np.savetxt(c + '_mean.txt',vertex_mean.astype('float32'),delimiter='\t',fmt='%f')
        np.savetxt(c + '_stdev.txt',vertex_std.astype('float32'),delimiter='\t',fmt='%f')
        save_mat(np.transpose(vertex_data), c, c + '.mat')

        if c_idx == 0:
            metric_dict[metric] = np.transpose(vertex_data.copy())
        else:
            metric_dict[metric] = np.concatenate(
            (metric_dict[metric], np.transpose(vertex_data.copy())), axis=0)
    save_mat(metric_dict[metric], 'wb_' + metric, 'wb_' + metric + '.mat')
