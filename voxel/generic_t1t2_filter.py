import os
import glob
import pandas as pd
import numpy as np
import sys
sys.path.append('/scratch/m/mchakrav/patelmo6/pnmf_bin/TractREC_p3/') #MODIFY to point to your cpy of TractREC
#from TractREC import * #my custom stuff
import TractREC as tr
import t1t2_outlier_correction as outlier
import argparse

parser=argparse.ArgumentParser(
    description='''Filtering of T1T2 values in subcortical regions to correct for outliers. Outputs one file per label, containing corrected values in a given region ''')
group = parser.add_argument_group(title="Execution options")

group.add_argument('--lookup', metavar='hi',help='csv file with two columns (Label_val, Label), containing label number and name pairings',required=True)
group.add_argument('--resn', type=float, help='resolution of t1/t2 images',required=True)
group.add_argument('--subjdir',help='directory containing t1t2 and label images for given subject',required=True)
group.add_argument('--t1t2_stem', help='stem of all t1t2 file names (eg. t1t2_stem = ratio.nii.gz for subject1_ratio.nii.gz',required=True)
group.add_argument('--label_stem', help='stem of all label file names',required=True)
group.add_argument('--left_gp_lb', help='label value for left gp', default=0, type=int)    
group.add_argument('--right_gp_lb', help='label value for right gp', default=0, type=int)

args=parser.parse_args()


lblegend=pd.read_csv(args.lookup)
lblegend=lblegend.set_index(['Label_val'])

resnn=args.resn
t1t2_stem=args.t1t2_stem
label_stem=args.label_stem
folder=os.path.abspath(args.subjdir)
t1t2_path=folder+"/*"+t1t2_stem
label_path=folder+"/*"+label_stem
t1t2_files=tr.natural_sort(glob.glob(t1t2_path))
lb_files=tr.natural_sort(glob.glob(label_path))

mean_dict={}
std_dict={}

t1t2_id=[]
for name in t1t2_files:
    t1t2_id.append(os.path.dirname(name))
print(resnn)
print(t1t2_stem)
print(t1t2_files)
print(lb_files)
print(t1t2_id)

mean_dict['Subject']=[os.path.basename(t1t2_id[0])]
std_dict['Subject']=[os.path.basename(t1t2_id[0])]
cols_list=[]
cols_list.append('Subject')
for lb in lblegend.index:
    roi=lblegend.loc[lb]['Label']
    corrected_fname=t1t2_id[0]+"/"+os.path.basename(t1t2_id[0])+"_"+str(roi)+"_filter.nii.gz"
    if os.path.exists(corrected_fname):
        continue
    roi=lblegend.loc[lb]['Label']
    cols_list.append(roi)
    print(lb,roi)
    print(corrected_fname)
    df_seg,res=tr.extract_quantitative_metric(t1t2_files,lb_files,IDs=t1t2_id[0],ROI_mask_files=None,label_df=lblegend,\
                                      label_subset_idx=lb,thresh_mask_files=None,thresh_val=None,\
                                      max_val=None,thresh_type=None,erode_vox=None,metric='data',DEBUG_DIR=None,\
                                      VERBOSE=False,USE_LABEL_RES=False,volume_idx=0)

    if((lb == args.left_gp_lb) | (lb == args.right_gp_lb)):
        sigma=1.27
    else:
        sigma=2.12
    print(sigma)
    corrected_data = outlier.outlier_correction(res[0].data,res[0].vox_coord, 25, sigma, resnn, t1t2_id[0])
    ref_fname = lb_files[0]
    d,a,h = tr.imgLoad(ref_fname,RETURN_HEADER=True)
    d_out = np.zeros_like(d).astype(np.float32)
    for idx, coord in list(enumerate(res[0].vox_coord[0])):
        d_out[coord[0]][coord[1]][coord[2]] = corrected_data[0,idx]

    tr.niiSave(corrected_fname,d_out,a,header=h,CLOBBER=True)
    mean_dict[str(roi)]=[np.mean(corrected_data)]
    std_dict[str(roi)]=[np.std(corrected_data)]

mean_df = pd.DataFrame.from_dict(mean_dict)
mean_df = mean_df[cols_list]
fname=folder + "/" + os.path.basename(t1t2_id[0]) + "_meant1t2.csv"
mean_df.to_csv(fname, index=False)

std_df = pd.DataFrame.from_dict(std_dict)
std_df = std_df[cols_list]
fname=folder + "/" + os.path.basename(t1t2_id[0]) + "_stdt1t2.csv"
std_df.to_csv(fname, index=False)



