#this script loads in raw ct .txt files from each subject and:
#1) concatenates each file to build a vertex X subj matrix of raw ct, for both left and right hemisphers
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

#Read in csv with subject demographics 
df_sorted = pd.read_csv('../../raw_data/sheets/20-03-14-McGillData_WH_Exprodo-Report_IncExc_CR_CRmed_cham_CRtopfdemeduc_civetpass_sorted.csv') #MODIFY to replace the .csv filename with the path to your sorted demographics csv file

#LOAD CIVET MASK TO IDENTIFY MIDLINE/ CORPOS COLLOSUM REGION
#IDENTIFY 'VALID VERTICES' - IE VERTICES NOT IN THIS REGION
left_mask = np.loadtxt('../surfsamp/mask_files/CIVET_2.0_mask_left_short.txt') #MODIFY to point to your left mask file
left_valid = np.where(left_mask==1) # list of valid indices in civet .txt file
left_invalid = np.where(left_mask==0) #list of invalid indices in civet .txt file

right_mask = np.loadtxt('../surfsamp/mask_files/CIVET_2.0_mask_right_short.txt') #MODIFY to point to your right mask file
right_valid = (np.where(right_mask==1))
right_invalid = (np.where(right_mask==0))
#38561 valid vertices

n_subjects=df_sorted.shape[0] #num rows in spreadsheet
n_vertex=40962

#Load left thicknesses into matrix left_ct
#do first row to create matrix, concat from there
row = df_sorted['oxmg_id'].tolist()[0] #MODIFY replace oxmg_id with the name of the column that contains subject id in your spreadsheet
fname="../civet/thickness/oxmg_" + str(row) + '_Struct_N_defaced_native_rms_rsl_tlaplace_30mm_left.txt' #MODIFY paths so that this evaluates to your civet thickness files
left_ct = np.loadtxt(fname).reshape(1,n_vertex)
print(np.shape(left_ct)) #1 x 40962

#load thickness file for rest of subjects and concatenate to make subjects X vertices matrix
for row in df_sorted['oxmg_id'].tolist()[1:]: #MODIFY replace oxmg_id with the name of the column that contains subject id in your spreadsheet
    fname="../civet/thickness/oxmg_" + str(row) + '_Struct_N_defaced_native_rms_rsl_tlaplace_30mm_left.txt' #MODIFY paths so that this evaluates to your civet thickness files
    x = np.loadtxt(fname).reshape(1,n_vertex)
    left_ct = np.concatenate((left_ct, x), axis=0)
print("raw left has", np.shape(left_ct)[0], "subjects", np.shape(left_ct)[1], "vertices")

#Repeat for right side
row = df_sorted['oxmg_id'].tolist()[0] #MODIFY replace oxmg_id with the name of the column that contains subject id in your spreadsheet
fname="../civet/thickness/oxmg_" + str(row) + '_Struct_N_defaced_native_rms_rsl_tlaplace_30mm_right.txt' #MODIFY paths so that this evaluates to your civet thickness files
right_ct = np.loadtxt(fname).reshape(1,n_vertex)
print(np.shape(right_ct)) #1 x 40962

#load thickness file for rest of subjects and concatenate to make subjects X vertices matrix
for row in df_sorted['oxmg_id'].tolist()[1:]: #MODIFY replace oxmg_id with the name of the column that contains subject id in your spreadsheet
    fname="../civet/thickness/oxmg_" + str(row) + '_Struct_N_defaced_native_rms_rsl_tlaplace_30mm_right.txt' #MODIFY paths so that this evaluates to your civet thickness files
    x = np.loadtxt(fname).reshape(1,n_vertex)
    right_ct = np.concatenate((right_ct, x), axis=0)
print("raw right has", np.shape(right_ct)[0], "subjects", np.shape(right_ct)[1], "vertices")

#valid ct matrices have dimensions subj x validvertices, masked midline vertices removed
n_rows = np.shape(left_ct)[0] #same for right side
n_cols = np.shape(left_valid)[1] #same for right side

#take all rows (subjects), but only columns which are valid vertices
left_ct_valid = left_ct[:, left_valid].reshape(n_rows,n_cols)
right_ct_valid = right_ct[:, right_valid].reshape(n_rows,n_cols)


#write out in .mat format
#nmf will want vertex X subj, so transpose
out_matrix = np.transpose(np.asmatrix(left_ct_valid))
fname = "left_ct_raw.mat"
print("Saving raw left ct to", fname,  "with shape", np.shape(out_matrix))
scipy.io.savemat(fname, {'X': out_matrix})
del out_matrix

out_matrix = np.transpose(np.asmatrix(right_ct_valid))
fname = "right_ct_raw.mat"
print("Saving raw right ct to", fname, "with shape", np.shape(out_matrix))
scipy.io.savemat(fname, {'X': out_matrix})
del out_matrix

#Concatenate the left and right hemisphere data to get whole brain data
wb_ct_valid = np.concatenate((np.transpose(left_ct_valid), np.transpose(right_ct_valid)),axis=0)

#write out whole brain data in .mat format
#write out z scored, un shifted data for stability analysis (to be shifted within each split later)
out_matrix = np.asmatrix(wb_ct_valid)
fname = "wholebrain_ct_raw.mat"
print("Saving whole brain ct to", fname, "with shape", np.shape(out_matrix))
scipy.io.savemat(fname, {'X': out_matrix})


