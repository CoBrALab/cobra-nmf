#In this script we append the H weightings to the demographics sheet for each subject
#Results in new columns being added to demographics sheet, new cols are
#Comp1_ct....compN_t1t2,comp1_dbm..compn_dbm

import os
import glob
import pandas as pd
import numpy as np

import sys
import pickle
import scipy
import scipy.stats
import argparse
from scipy.io import savemat, loadmat
#Read in csv with subject demographics 

parser=argparse.ArgumentParser(
    description='''This script reads in nmf results and outputs a .csv containing component weights for each subject''')

group = parser.add_argument_group(title="Execution options")

group.add_argument(
    '--nmf_results', help='.mat file containing nmf results',required=True)
group.add_argument(
    "--metrics",help='''metrics used in nmf analysis. output will have 1 
    column for every component metric pair''', metavar='list', nargs='+', required=True)
group.add_argument(
    '--demo_csv', help='demographic spreadsheet, must contain subject id',required=True)
group.add_argument(
    '--id_col', help='name of subject Id column in demographic sheet',required=True)

args=parser.parse_args()

#Read in csv with subject demographics 
df_sorted = pd.read_csv(args.demo_csv)

#this function appends NMF weights to a demographics df
#order of new cols is Comp1_ct....compN_t1t2,comp1_dbm..compn_dbm
def append_subjweights_plsstyle(df_demo,nmf_weights, metrics):
    df_demo_nmf=df_demo.copy()
    n_subjects=len(df_demo)
    maxrow=np.shape(nmf_weights)[0]
    
    for comp in range(0,maxrow):
        
        for m in range(0,len(metrics)):
            col='Comp'+str(comp+1)+'_'+metrics[m]
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
#nmf_res_filename = "sample_nmf_res.mat" #MODIFY point to your nmf results
H = loadmat(args.nmf_results)['H']
print(np.shape(H)) #check shape

#use append_subjweights_plsstyle to concat the demographic df with the nmf weights
df_sorted_nmfweights = append_subjweights_plsstyle(df_sorted, H, args.metrics)

numcomps=np.shape(H)[0]
fname_alldata = "demographics_and_nmfweights_k" + str(numcomps) + '.csv' 
print("saving all columns to ", fname_alldata)
df_sorted_nmfweights.to_csv(fname_alldata, index = False) 

#the above writes out your full demographic spreadsheet and the nmf weights
#for pls it would be easier to have, perhaps, just subject ID and nmf weight columns. below does this

#build list of column headers of interest - ie keep subject id and comp weights only for nmf csv
comp_cols = [args.id_col] #set comp_cols to start as the non nmf columns of interest. eg this is my subject id column, so im saying i want to end up with subject ID only and component weights. add age, sex, etc, if you want
for comp in range(1, np.shape(H)[0]+1):
    for m in args.metrics:
        comp_cols.append("Comp" + str(comp) + "_" + m)

#create a new df to write out which contains only the cols of interest
df_nmfweights_pls = df_sorted_nmfweights[comp_cols].copy()
fname_pls_data="Hweights_k" + str(numcomps) + '.csv'
print("saving select columns", comp_cols, "to", fname_pls_data)
df_nmfweights_pls.to_csv(fname_pls_data, index = False) 

