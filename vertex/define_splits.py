## LOAD MODULES/SOFTWARE
import os
import glob
import pandas as pd
import numpy as np

import sys
import pickle
import scipy
import scipy.stats
from scipy.io import savemat, loadmat
import sklearn
from sklearn.model_selection import StratifiedShuffleSplit
import argparse

parser=argparse.ArgumentParser(
    description='''This script creates stratified input matrices for stability analysis,
    stores outputs as .mat files in stability_splits directory''')

group = parser.add_argument_group(title="Execution options")

group.add_argument(
    '--demo_csv', help='demographic spreadsheet, must contain subject id',required=True)
group.add_argument(
    '--id_col', help='name of subject Id column in demographic sheet',required=True)

group.add_argument('--n_folds', help='number of folds', type=int, default=10)

group.add_argument(
    "--inputs",help="metric matrices to stratify", metavar='list', nargs='+', required=True)

group.add_argument(
    "--stratifyby",help="demographic variables to stratify splits by", metavar='list', nargs='+', required=True)


args=parser.parse_args()

def save_mat(x,key,fname):
    print("Saving ", np.shape(x), key, "to", fname)
    scipy.io.savemat(fname, {'X': x})

#read in demographic spreadsheet with subject ids, age, prisma etc
df_sorted = pd.read_csv(args.demo_csv)

## create demo matrix containing subj id and the variables to stratify by (age)

demo_vars = []
demo_vars.append(args.id_col)
for x in args.stratifyby:
    demo_vars.append(x)
demo = df_sorted[demo_vars].values

#define train data as subj ids (x)
#define categorical vars as vars to stratify by (y, ie labels)
X = demo[:,0:1]
y = demo[:,1:]


#use sklearn Stratify tools to generate stratified splits of data
n_folds=args.n_folds
sss = StratifiedShuffleSplit(n_splits=n_folds, test_size=0.5, random_state=0)
sss.get_n_splits(X, y)

Asplits_indices = {}; Asplits_subjectIDs = {}  #dicts for storing indices and corresponding subj ids
Bsplits_indices = {}; Bsplits_subjectIDs = {}

iter=0
#cycle through train test splits, add to above dictionaries
for train_index, test_index in sss.split(X, y):
    Asplits_indices[str(iter)] = train_index;
    Bsplits_indices[str(iter)] = test_index;
    
    ID_list = []
    s = train_index[0]
    ID_list.append(df_sorted[args.id_col].iloc[s])
    for s in train_index[1:]:
        ID_list.append(df_sorted[args.id_col].iloc[s])
    Asplits_subjectIDs[iter] = ID_list

    ID_list = []
    s = test_index[0]
    ID_list.append(df_sorted[args.id_col].iloc[s])
    for s in test_index[1:]:
        ID_list.append(df_sorted[args.id_col].iloc[s])
    Bsplits_subjectIDs[iter] = ID_list

    iter = iter + 1

#save splits
out_dir = "stability_splits/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

pickle.dump(Asplits_indices, open( out_dir + "Asplits_indices.p", "wb"))
pickle.dump(Bsplits_indices, open( out_dir + "Bsplits_indices.p", "wb" ))

input_list = args.inputs

data_dict={}
for f in input_list:
    data_dict[f] = loadmat(f)['X'] #load raw data
data_dict.keys()

#for each split, build required indices
for split in range(0, args.n_folds):

    #get data from first metric
    metric=input_list[0]
    data_all = data_dict[metric]
    #get data_a and data_b, containing ct data for A indicies and B indices
    data_a = data_all[:,Asplits_indices[str(split)]]; data_b = data_all[:,Bsplits_indices[str(split)]]
    #z score each 
    a_mx_wb = scipy.stats.zscore(data_a,axis=None)
    b_mx_wb = scipy.stats.zscore(data_b,axis=None)

    #repeat for each metric 
    for metric in input_list[1:]:
        data_all = data_dict[metric]
        data_a = data_all[:,Asplits_indices[str(split)]]; data_b = data_all[:,Bsplits_indices[str(split)]]
        data_a_z = scipy.stats.zscore(data_a,axis=None) #zscore
        data_b_z = scipy.stats.zscore(data_b,axis=None)
        a_mx_wb = np.concatenate((a_mx_wb,data_a_z),axis=1) #append z scored data for this metric to the rest
        b_mx_wb = np.concatenate((b_mx_wb,data_b_z),axis=1)

    #shift each to be non negative
    a_mx_shift_wb = a_mx_wb - np.min(a_mx_wb)
    b_mx_shift_wb = b_mx_wb - np.min(b_mx_wb)

    #write out
    save_mat(a_mx_shift_wb, 'split a_' + str(split), out_dir + "a_" + str(split) + ".mat")
    save_mat(b_mx_shift_wb, 'split b_' + str(split), out_dir + "b_" + str(split) + ".mat")

