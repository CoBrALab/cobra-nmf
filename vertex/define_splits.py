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

#read in demographic spreadsheet with subject ids, age, prisma etc
df_sorted = pd.read_csv('../../../../raw_data/sheets/07-04-20-McGillData_WH_Exprodo-Report_IncExc_CR_CRmed_cham_CRtopfdemeduc_civetpass_slopes_sorted.csv') #MODIFY

## create demo matrix containing subj id and the variables to stratify by (age)
#variables to stratify by - probably age and sex, maybe disease group as well
age_var='OX.AGE'; sex_var= 'OX.SEX' #MODIFY
#group_var = ?? #MODIFY if applicable, add others if applicable

demo = df_sorted[['oxmg_id',age_var,sex_var]].values #MODIFY to be subjectID and the variables you want to stratify by (probably age_var, sex_var, plus others defined above)

#sklearn stratify tool needs categorical vars - bin age by median val
demo = np.hstack((demo, np.zeros((len(df_sorted),1)))) #new col for age
median_age = np.median(df_sorted[age_var].values) 
#<= median age -> you are young 
demo[np.where(df_sorted[age_var].values > median_age),-1]="old"
demo[np.where(df_sorted[age_var].values <= median_age),-1]="young"

#define train data as subj ids (x)
#define categorical vars as vars to stratify by (y, ie labels)
X = demo[:,0:1]
y = demo[:,2:] #this skips [:,1] which is numerical age, instead includes categorical age

#use sklearn Stratify tools to generate stratified splits of data
n_splits=10 #MODIFY
sss = StratifiedShuffleSplit(n_splits=n_splits, test_size=0.5, random_state=0)
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
    ID_list.append(df_sorted['oxmg_id'].iloc[s]) #MODIFY oxmg_id to subject id name
    for s in train_index[1:]:
        ID_list.append(df_sorted['oxmg_id'].iloc[s]) #MODIFY oxmg_id to subject id name
    Asplits_subjectIDs[iter] = ID_list
    
    ID_list = []
    s = test_index[0]
    ID_list.append(df_sorted['oxmg_id'].iloc[s]) #MODIFY oxmg_id to subject id name
    for s in test_index[1:]:
        ID_list.append(df_sorted['oxmg_id'].iloc[s]) #MODIFY oxmg_id to subject id name
    Bsplits_subjectIDs[iter] = ID_list
    
    iter = iter + 1
  
#save splits
out_dir = "stability_splits/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

pickle.dump(Asplits_indices, open( out_dir + "Asplits_indices.p", "wb"), protocol=0)
pickle.dump(Asplits_subjectIDs, open( out_dir + "Asplits_subjectIDs.p", "wb" ))
pickle.dump(Bsplits_indices, open( out_dir + "Bsplits_indices.p", "wb" ))
pickle.dump(Bsplits_subjectIDs, open( out_dir + "Bsplits_subjectIDs.p", "wb" ))


input_list=["wholebrain_t1t2","wholebrain_dbm"] #MODIFY to match the filenames of the wholebrain files created in your sample_extract*py scripts. should be name of the file without the .mat extension

data_dict={}
raw_data_dir = '/folder/where/wholebrain/rawdata/stored/' #MODIFY to location of wholebrain .mat files containig raw data
for file in input_list:
    fname = raw_data_dir + '/' + file + ".mat"
    data_dict[file] = scipy.io.loadmat(fname)['X'] #load raw data

#for each split, build required indices
for split in range(0, n_splits):

    #get data from first metric
    metric=input_list[0]
    data_all = data_dict[metric]
    #get data_a and data_b, containing ct data for A indicies and B indices
    data_a = data_all[:,Asplits_indices[str(split)]]; data_b = data_all[:,Bsplits_indices[str(split)]]
    #z score each 
    a_mx_wb = scipy.stats.zscore(data_a,axis=None)
    b_mx_wb = scipy.stats.zscore(data_b,axis=None)

    #repeat for each metric 
    #if only one metric, skip the for loop below, go straight to line 108 #
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
    fname = out_dir + "a_" + str(split) + ".mat"
    print(fname, np.min(a_mx_shift_wb), np.shape(a_mx_shift_wb))
    scipy.io.savemat(fname, {'X': a_mx_shift_wb})

    fname = out_dir + "b_" + str(split) + ".mat"
    print(fname, np.min(b_mx_shift_wb), np.shape(b_mx_shift_wb))
    scipy.io.savemat(fname, {'X': b_mx_shift_wb})

