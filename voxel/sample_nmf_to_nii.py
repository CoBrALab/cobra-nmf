#writing/mapping values to nifti file for visualization
#main steps:
#1: load some label/image data into python, as you did for building input matrix. used as a coord space/ref for the output label
#2. load in .mat results from pnmf
#3. map W component scores (or a clustering) to output label - ie take the array loaded in from a nifti file, replace values with component scores/cluster labels, write file out 

#load libraries
import os
import glob
import pandas as pd
import numpy as np
import scipy
import sys
import hdf5storage
import argparse

options = hdf5storage.Options(oned_as = 'column', matlab_compatible = True, action_for_matlab_incompatible = 'error')
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
    '--id_col', help='name of subject Id column in demographic sheet',required=True)
group.add_argument(
    '--data_dir', help='directory containing subject data folders', required=True)

group.add_argument('--tract_rec', help='path to TractREC library', required=True)
group.add_argument('--nmf_components', help='.mat output results from nmf', required=True)
group.add_argument('--out_dir', help='output directory', required=True)

args=parser.parse_args()

#this function takes in component scores and writes out a nifti label file, using the coordinates extracted via tract rec
# voxscores - n_voxels X n_components matrix/array.  function will write out one label for each column
# base - a string specifying folder/path base for output file.
# lb_ref - the majority voted label, used as a reference image for the output label. the fcn makes an image 'like' lb_ref, but replaces values with component scores/cluster lbels
# hemi - a string included in the output filename of the label. useful if you need to specify left/right, less useful otherwise
#refimg_res - a res object from Tractrect. will contain voxel coordinate info in the vox_coord variable. required
def voxelscores_to_label(voxscores, base, lb_ref, hemi,refimg_res):
    voxcoords=refimg_res[0].vox_coord
    compnum  = np.shape(voxscores)[1]
    for c in range(0,compnum):
        outpath = base + hemi + 'hc-' + str(c) + '.nii.gz'
        print(outpath)
        b=np.transpose(voxscores[:,c])
        print(np.shape(b))
        d,a,h = tr.imgLoad(lb_ref,RETURN_HEADER=True)
        d_out = np.zeros_like(d).astype(np.float32)
        for idx, coord in list(enumerate(voxcoords[0])):
            d_out[coord[0]][coord[1]][coord[2]] = b[idx]

        tr.niiSave(outpath,d_out,a,header=h,CLOBBER=True)
        
#load TractREC 
sys.path.append(args.tract_rec)
import TractREC as tr

working_dir=args.data_dir

#you'll need a lookup table in .csv format with label number and structure
all_label_seg=pd.read_csv(args.lookup)
all_label_seg=all_label_seg.set_index(['Label_val']) #make the index column the actual index
all_lobule_subset_idx=[args.mask_label] 
df_sorted = pd.read_csv(args.demo_csv)

#Create two lists, each containing the filenames of the t1t2 filtered files and majority voted labels
metric_files=[]
fused_label_files=[]
for row in df_sorted[args.id_col]:
    fname = working_dir  + str(row) + '/*' + args.metric_stem
    metric_files.append(glob.glob(fname)[0])
    fname = working_dir  + str(row) + '/*' + args.label_stem 
    fused_label_files.append(glob.glob(fname)[0])
    
metric_IDs=[]
for name in metric_files:
    metric_IDs.append(os.path.dirname(name))



#use TractREC to extract the t1t2 data for each subject, only in voxels overlaying with the label
df_seg,res=tr.extract_quantitative_metric(metric_files,fused_label_files,IDs=metric_IDs,ROI_mask_files=None,label_df=all_label_seg,\
                                      label_subset_idx=all_lobule_subset_idx,thresh_mask_files=None,thresh_val=None,\
                                      max_val=None,thresh_type=None,erode_vox=None,metric='data',DEBUG_DIR=None,\
                                      VERBOSE=True,USE_LABEL_RES=False,volume_idx=0)

#the res variable is now a list of the results for each subject

#load in the results. in particular we want W here, which is voxels X components
w = hdf5storage.loadmat(args.nmf_components)['W']
h = hdf5storage.loadmat(args.nmf_components)['H']

num_components = np.shape(w)[1]
res_dir = args.out_dir
if not os.path.exists(res_dir):
    os.makedirs(res_dir)

lb_ref = fused_label_files[0] #assumes all labels the same

#the fcn call on the next line should create/write one nifti file for each component
voxelscores_to_label(w,res_dir,lb_ref,hemi,res) 

#instead of one file per component, the code below uses a winner take all strategy to create cluster labels, and writes out the file

w_label = np.zeros((np.shape(w)[0],1))
np.shape(w_label)
for v in range(0,np.shape(w_label)[0]):
    w_label[v,0] = np.argmax(w[v,:])
print(np.shape(w_label))
print(np.max(w_label))
w_label_out = w_label + 1

voxelscores_to_label(w_label_out,res_dir,lb_ref,hemi + '-cluster' + str(compnum),res) 


