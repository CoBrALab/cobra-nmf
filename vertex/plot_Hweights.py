## LOAD MODULES/SOFTWARE
import os
import glob
import pandas as pd
import numpy as np

import sys
import pickle
import hdf5storage
import scipy
from scipy import stats

import matplotlib as mpl
from matplotlib import cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedFormatter, FixedLocator
from matplotlib import pyplot as plt
import argparse
options = hdf5storage.Options(oned_as = 'column', matlab_compatible = True, action_for_matlab_incompatible = 'error')


parser=argparse.ArgumentParser(
    description='''This script outputs a .png file containing a heatmap of input nmf matrix data''')

group = parser.add_argument_group(title="Execution options")

group.add_argument(
    '--nmf_weights', help='.mat file containing nmf results',required=True)
group.add_argument(
    '--output', help='output .png filename',required=True)
group.add_argument(
    '--minimum', type=float, help='min value',required=False,default=-2)
group.add_argument(
    '--maximum', type=float, help='max value',required=False,default=2)
group.add_argument(
    '--width', type=float, help='figure width',required=False,default=16)
group.add_argument(
    '--height', type=float, help='figure height',required=False,default=8)

args=parser.parse_args()

h=hdf5storage.loadmat(args.nmf_weights)['H']

#heat mapping for H matrix
def heatmapping(data,minn,maxx,cbar_tix,fig_width,fig_height,title='',fname=''):
    import matplotlib as mpl
    from matplotlib import cm
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedFormatter, FixedLocator
    plt.rc('figure', titlesize=30)  # fontsize of the figure title
    #Linearly interpoalte a colour gradient 
   
    viridis = cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    cmap = mpl.colors.ListedColormap(newcolors)
    img = plt.imshow(data,interpolation='nearest', \
    cmap = cmap, origin='upper',vmin=minn,vmax=maxx)
    #Set the axis of the plot so it isn't a long rectangle
    ax = plt.gca()
    ax.set_aspect('auto') #use 'auto' for automatic aspect
    ax.tick_params(axis='both',which='both',bottom='off',top='off',labelbottom='on',left='on',labelleft='on', pad = 5)
    ax.set_xlabel('')
    ax.set_ylabel('Voxels', fontsize=40)
    ax.yaxis.set_ticklabels([])
    ax.yaxis.labelpad = 5
    ax.tick_params(axis='y',size=15)
    ax.grid(False)
    fig = plt.gcf()
    fig.set_size_inches(fig_width,fig_height)
    cbar = plt.colorbar(img,cmap=cmap)
    
    cbar.set_ticks(np.arange(minn, maxx, cbar_tix))
    cbar.ax.tick_params(labelsize=30)
    if title:
        plt.title(title, fontsize=30)
    plt.savefig(fname, bbox_inches='tight')
    

h_z=scipy.stats.zscore(h,axis=1)
heatmapping(h_z,args.minimum,args.maximum+0.0001,2,args.width,args.height,title="Hweights",fname=args.output)    
