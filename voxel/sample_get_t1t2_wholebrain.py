#load libraries
import os
import glob
import pandas as pd
import numpy as np
import scipy
from scipy.io import savemat

import sys

#load TractREC 
sys.path.append('/data/chamal/projects/raihaan/HCP/HCP900/scripts/TractREC/working/TractREC/TractREC/') #MODIFY
import TractREC as tr

working_dir="/data/scratch/raihaan/hcp-micro/329subjectNMF_singleshell/warped/" #MODIFY

#you'll need a lookup table in .csv format with label number and structure. create one for your labels (left right striatum)
#look at the .csv file below as an example, then modify the line below to point to your .csv lookup
all_label_seg_file='/data/chamal/projects/raihaan/projects/inprogress/hc-nmf-micro/raw_data/sheets/2015_09_labels_CB_Hipp_Subcort_63.csv' #MODIFY
all_label_seg=pd.read_csv(all_label_seg_file)
all_label_seg=all_label_seg.set_index(['Index']) #make the index column the actual index
all_lobule_subset_idx=[24] #MODIFY - this should be the label value of the structure of interest eg if this is extracting left striatum, make this the left striatum label val
all_label_seg.head()
all_lobule_subset_idx
df=pd.read_csv('/data/chamal/projects/raihaan/projects/inprogress/hc-nmf-micro/raw_data/sheets/df_sorted_unrestricted_329.csv') #MODIFY to point to your copy of the demographics file
df_sorted=df

#Create two lists, each containing the filenames of the t1t2 filtered files and majority voted labels
t1divt2_files=[]
fused_label_files=[]
for row in df_sorted['Subject']:
    fname = working_dir  + str(row) + '/lefthc_correctedt1t2.nii.gz' #MODIFY
    t1divt2_files.append(glob.glob(fname)[0])
    fname = working_dir  + str(row) + '/majvote_hccorrected_label.nii.gz' #MODIFY
    fused_label_files.append(glob.glob(fname)[0])
    
t1divt2_IDs=[os.path.dirname(name.split("t1t2")[0]) for name in t1divt2_files]

#use TractREC to extract the t1t2 data for each subject, only in voxels overlaying with the label
df_seg,res=tr.extract_quantitative_metric(t1divt2_files,fused_label_files,IDs=t1divt2_IDs,ROI_mask_files=None,label_df=all_label_seg,\
                                      label_subset_idx=all_lobule_subset_idx,thresh_mask_files=None,thresh_val=None,\
                                      max_val=None,thresh_type=None,erode_vox=None,metric='data',DEBUG_DIR=None,\
                                      VERBOSE=True,USE_LABEL_RES=False,volume_idx=0)
#the res variable is now a list of the results for each subject

#create a matrix that is, for now, just containing the first subjects data in the res list (index 0)
t1t2_stack=res[0].data 

#cycle through the remaining subjects in res (index 1 - end)
#iteratively concatenate to t1t2_stack
#before this loop, t1t2_stack has shape (1, voxels)
#after 1 iteration, t1t2_stack has shape (2, voxels)....ends at (329, voxels)
for img in range(1,len(t1divt2_IDs)):
    t1t2_stack=np.concatenate((t1t2_stack,res[img].data),axis=0)    

#Save the matrix in .mat format. nmf wants voxels X subjects, so save the transpose of t1t2_stack
t1t2_out = np.transpose(t1t2_stack)
print(np.shape(t1t2_out))
scipy.io.savemat('/path/to/save/wholebrain_t1t2.mat', mdict={'X': t1t2_out}) #MODIFY PATH 

