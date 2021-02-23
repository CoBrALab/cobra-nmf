#In this script we append the H weightings to the demographics sheet for each subject
#Results in new columns being added to demographics sheet, new cols are
#Comp1_ct....compN_t1t2,comp1_dbm..compn_dbm

import os
import glob
import pandas as pd
import numpy as np

import sys
import scipy
import scipy.stats
import scipy.io
from scipy.io import savemat, loadmat

#Read in csv with subject demographics 
df_sorted = pd.read_csv('../../../raw_data/sheets/06-05-20-McGillData_WH_Exprodo-Report_IncExc_CR_CRmed_cham_CRtopfdemeduc_civetpass_slopes_time.age_sorted.csv') #MODIFY

metrics=["CT", "SA", "GI"] #MODIFY to be a list of the metrics used
                           #the exact text used is up to you, but there should be one string for every metric
                           #the order matters - if the first chunk of data is ct, then put the ct string first

#this function appends NMF weights to a demographics df
#order of new cols is Comp1_ct....compN_t1t2,comp1_dbm..compn_dbm
def append_subjweights_plsstyle(df_demo,nmf_weights, metrics):
    df_demo_nmf=df_demo.copy()
    n_subjects=len(df_demo)
    maxrow=np.shape(nmf_weights)[0]
    
    for comp in range(0,maxrow):
        
        for m in range(0,len(metrics)):
            col='Comp'+str(comp+1)+metrics[m]
            df_demo_nmf[col] = nmf_weights[comp,n_subjects*m:n_subjects*(m+1)]
      
    print( "df had", len(df_demo.columns), "columns")
    print("numcomps is",  np.shape(nmf_weights)[0], "add", (len(metrics)*maxrow), "columns")
    print("df_demo has", len(df_demo_nmf.columns), "columns")
    
    if ((len(df_demo_nmf.columns) - len(df_demo.columns)) == (len(metrics)*maxrow)):
            return df_demo_nmf
    else:
        print("ERROR")
        return
    
#load in the nmf results of interest, check shape
nmf_res_filename = "sample_nmf_res.mat" #MODIFY point to your nmf results
H = loadmat(nmf_res_filename)['H']
print(np.shape(H)) #check shape

#use append_subjweights_plsstyle to concat the demographic df with the nmf weights
df_sorted_nmfweights = append_subjweights_plsstyle(df_sorted, H)
fname_alldata = "demographics_and_nmfweights_k10.csv" #MODIFY filename if you want
print("saving all columns to ", fname_alldata)
df_sorted_nmfweights.to_csv(fname_alldata, index = False) 

#the above writes out your full demographic spreadsheet and the nmf weights
#for pls it would be easier to have, perhaps, just subject ID and nmf weight columns. below does this

#build list of column headers of interest - ie keep subject id and comp weights only for nmf csv
comp_cols = ['Subject'] #set comp_cols to start as the non nmf columns of interest. eg this is my subject id column, so im saying i want to end up with subject ID only and component weights. add age, sex, etc, if you want
for comp in range(1, np.shape(H)[0]+1):
    for m in metrics:
        comp_cols.append("Comp" + str(comp) + "_" + m)

#create a new df to write out which contains only the cols of interest
df_nmfweights_pls = df_sorted_nmfweights[comp_cols].copy()
fname_pls_data="pls_hweights_k10.csv" #MODIFY filename to appropriate name
print("saving select columns", comp_cols, "to", fname_pls_data)
df_nmfweights_pls.to_csv(fname_pls_data, index = False) 

#for pls it is also good to have a csv that is just subject ID and the cognitive data of interest
cog_cols = ['Subject','Age','Sex','Memory','ExecFunction'] #MODIFY to be the column names of subject id and the behavioural variables of interest
df_cog_pls = df_sorted_nmfweights[cog_cols].copy()
fname_pls_data="pls_cogdata.csv" #MODIFY filename to appropriate name
print("saving select columns", comp_cols, "to", fname_pls_data)
df_cog_pls.to_csv(fname_pls_data, index = False) 

