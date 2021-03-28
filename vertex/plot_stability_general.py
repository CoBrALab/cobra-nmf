import os
import glob
import pandas as pd
import numpy as np
import sys
import matplotlib as mpl
from matplotlib import cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedFormatter, FixedLocator
from matplotlib import pyplot as plt
plt.switch_backend('Agg')
import argparse

parser=argparse.ArgumentParser(
    description='''This script takes in computed stability metrics and plots
    number of components on the x axis, stability coefficient and graident of recon error on y''')

group = parser.add_argument_group(title="Execution options")

group.add_argument('--spacing', help='interval between computed ks', type=int, default=1)

group.add_argument(
    "--stability_correlations",help="csv containing computed stability metrics", required=True)

group.add_argument(
    "--output",help="output .png file to store plot", required=True)


args=parser.parse_args()

#plot stability and error gradient on same plot 

#df_stab = pd.read_csv(sys.argv[1])
df_stab = args.stability_correlations

#fix formatting of recon error values - remove [[ and ]] from start/end
reconA_corrected=[] #list of corrected reconA vals
#cycle through recon a column, extract the  characters starting at index 3, and drop the last 2 (ie [2:-2])
#append to corrected list, make float
for val in df_stab['Recon_errorA'].values:
    reconA_corrected.append((float(val[2:-2])))
#repeat for reconB
reconB_corrected=[] #list of corrected reconB vals
for val in df_stab['Recon_errorB'].values:
    reconB_corrected.append((float(val[2:-2])))
#add these columns to your df_stab
df_stab['Recon_errorA_corrected'] = reconA_corrected
df_stab['Recon_errorB_corrected'] = reconB_corrected

max_gran = np.max(df_stab['Granularity'].values)
min_gran = np.min(df_stab['Granularity'].values)
#interval=2 #MODIFY - this should represent the spacing of granularities investigated. ie if 2,4,5,8... interval = 2
interval = args.spacing

split_corr = []
for g in np.arange(min_gran,max_gran+1,interval):
    split_corr.append(df_stab.loc[df_stab['Granularity'] == g][['Corr_mean']].values)
    
plt_arr = np.zeros((1,np.shape(split_corr)[0]))
plt_std_arr = np.zeros((1,np.shape(split_corr)[0]))
for g in range(0,int(((max_gran - min_gran)/interval) + 1)):
    plt_arr[0,g] = np.mean(split_corr[g])
    plt_std_arr[0,g] = np.std(split_corr[g])
    
dict_errorgrad = {'Granularity' : np.arange(min_gran + interval,max_gran+1,interval).flatten()}
dict_errorgrad['Index'] = np.arange(0, np.shape(dict_errorgrad['Granularity'])[0], 1)
for iter in range(1,11):
    dict_errorgrad["A_iter" + str(iter)] = np.diff(df_stab.loc[df_stab['Iteration'] == iter][['Recon_errorA_corrected']].values.flatten(), axis=0).tolist()
    dict_errorgrad["B_iter" + str(iter)] = np.diff(df_stab.loc[df_stab['Iteration'] == iter][['Recon_errorB_corrected']].values.flatten(), axis=0).tolist()
df_errorgrad = pd.DataFrame(data=dict_errorgrad, index = np.arange(1,np.shape(dict_errorgrad['Granularity'])[0]+1).flatten())

error_grad_arr = np.zeros((1,np.shape(np.arange(min_gran + interval,max_gran+1,interval))[0]))
error_grad_std_arr = np.zeros((1,np.shape(np.arange(min_gran + interval,max_gran+1,interval))[0]))
for idx in range(0,int((max_gran - min_gran)/interval)):
    error_grad_arr[0,idx] = np.mean(df_errorgrad.iloc[idx,2:])
    error_grad_std_arr[0,idx] = np.std(df_errorgrad.iloc[idx,2:])
    
stab_x = np.arange(min_gran,max_gran+1,interval)
grad_err_x = np.arange(min_gran+interval,max_gran+1,interval)
fig, ax1 = plt.subplots(figsize=(16, 8), dpi=100)

color = 'tab:red'
ax1.set_xlabel('Number of Components', fontsize = 25)
ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(5))
ax1.set_ylabel('Stability Coefficient', color=color, fontsize = 25)
ax1.errorbar(stab_x,plt_arr.flatten(),yerr=plt_std_arr.flatten(), c=color, marker=".", lw=2, ms = 10)

ax1.tick_params(axis='y', labelcolor=color, labelsize=20)
ax1.tick_params(axis='x', labelsize=20)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Gradient(Reconstruction Error)', color=color, fontsize = 25)  # we already handled the x-label with ax1
ax2.errorbar(grad_err_x,error_grad_arr.flatten(),yerr=error_grad_std_arr.flatten(), c=color, marker=".", lw=2, ms = 10)

ax2.tick_params(axis='y', labelcolor=color, labelsize=20)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.title('NMF Stability', fontsize=30)
plt.savefig(args.output, dpi='figure', bbox_inches='tight')
