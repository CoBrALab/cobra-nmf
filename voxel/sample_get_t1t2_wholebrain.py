#load libraries
import os
import glob
import pandas as pd
import numpy as np
import scipy
from scipy.io import savemat
import argparse
import sys

parser=argparse.ArgumentParser(
    description='''This script extracts voxel data from nifti files and outputs a voxel x subject matrix
    in .mat format''')

parser.add_argument(
    "--metric",help="metric to extract", metavar='T1T2',required=True)
group = parser.add_argument_group(title="Execution options")

group.add_argument(
    '--lookup', metavar='.csv',help='''csv file with two columns (Label_val, Label)
    containing label number and name pairings''',required=True)
group.add_argument(
    '--metric_stem', help='''stem of all t1t2 file names (eg. t1t2_stem =
    ratio.nii.gz for subject1_ratio.nii.gz''',required=True)
group.add_argument(
    '--label_stem', help='stem of all label file names',required=True)
group.add_argument(
    '--mask_label', type=int, help='label value identifying voxels of interest',required=True)
group.add_argument(
    '--demo_csv', help='demographic spreadsheet, must contain subject id',required=True)
group.add_argument(
    '--data_dir', help='directory containing subject data folders', required=True)
 
group.add_argument('--tract_rec', help='path to TractREC library', required=True)

args=parser.parse_args()



#load TractREC 
sys.path.append(args.tract_rec) #MODIFY
import TractREC as tr

working_dir=args.data_dir
#working_dir="/data/scratch/raihaan/hcp-micro/329subjectNMF_singleshell/warped/" #MODIFY

#you'll need a lookup table in .csv format with label number and structure. create one for your labels (left right striatum)
#look at the .csv file below as an example, then modify the line below to point to your .csv lookup
#all_label_seg_file='/data/chamal/projects/raihaan/projects/inprogress/hc-nmf-micro/raw_data/sheets/2015_09_labels_CB_Hipp_Subcort_63.csv' #MODIFY
#all_label_seg_file=args.lookup
all_label_seg=pd.read_csv(args.lookup)
all_label_seg=all_label_seg.set_index(['Label_val']) #make the index column the actual index
#all_lobule_subset_idx=[24] #MODIFY - this should be the label value of the structure of interest eg if this is extracting left striatum, make this the left striatum label val
all_lobule_subset_idx=[args.mask_label] #MODIFY - this should be the label value of the structure of interest eg if this is extracting left striatum, make this the left striatum label val

#df=pd.read_csv('/data/chamal/projects/raihaan/projects/inprogress/hc-nmf-micro/raw_data/sheets/df_sorted_unrestricted_329.csv') #MODIFY to point to your copy of the demographics file
#df_sorted=df
df_sorted = pd.read_csv(args.demo_csv)

#Create two lists, each containing the filenames of the t1t2 filtered files and majority voted labels
metric_files=[]
fused_label_files=[]
for row in df_sorted['Subject']:
    #fname = working_dir  + str(row) + '/lefthc_correctedt1t2.nii.gz' #MODIFY
    fname = working_dir  + str(row) + '/*' + args.metric_stem
    metric_files.append(glob.glob(fname)[0])
    #fname = working_dir  + str(row) + '/majvote_hccorrected_label.nii.gz' #MODIFY
    fname = working_dir  + str(row) + '/*' + args.label_stem 
    fused_label_files.append(glob.glob(fname)[0])
    
#metric_IDs=[os.path.dirname(name.split("t1t2")[0]) for name in t1divt2_files]
metric_IDs=[]
for name in metric_files:
    metric_IDs.append(os.path.dirname(name))

#use TractREC to extract the t1t2 data for each subject, only in voxels overlaying with the label
df_seg,res=tr.extract_quantitative_metric(metric_files,fused_label_files,IDs=metric_IDs,ROI_mask_files=None,label_df=all_label_seg,\
                                      label_subset_idx=all_lobule_subset_idx,thresh_mask_files=None,thresh_val=None,\
                                      max_val=None,thresh_type=None,erode_vox=None,metric='data',DEBUG_DIR=None,\
                                      VERBOSE=True,USE_LABEL_RES=False,volume_idx=0)
#the res variable is now a list of the results for each subject

#create a matrix that is, for now, just containing the first subjects data in the res list (index 0)
metric_stack=res[0].data 

#cycle through the remaining subjects in res (index 1 - end)
#iteratively concatenate to t1t2_stack
#before this loop, t1t2_stack has shape (1, voxels)
#after 1 iteration, t1t2_stack has shape (2, voxels)....ends at (329, voxels)
for img in range(1,len(metric_IDs)):
    metric_stack=np.concatenate((metric_stack,res[img].data),axis=0)    

#Save the matrix in .mat format. nmf wants voxels X subjects, so save the transpose of t1t2_stack
metric_out = np.transpose(metric_stack)
print(np.shape(metric_out))
scipy.io.savemat('raw_' + args.metric + '.mat', mdict={'X': metric_out}) #MODIFY PATH 

