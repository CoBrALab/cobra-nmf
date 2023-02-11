## LOAD MODULES/SOFTWARE
import os
import glob
import pandas as pd
import numpy as np

import sys
import pickle
import scipy
from scipy.io import savemat, loadmat
from scipy import stats

import matplotlib as mpl
from matplotlib import cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedFormatter, FixedLocator
from matplotlib import pyplot as plt

import argparse

parser=argparse.ArgumentParser(
    description='''This script outputs a .png file containing a heatmap of input nmf matrix data''')

group = parser.add_argument_group(title="Execution options")

group.add_argument(
    '--input', help='.mat file containing voxel x subjects data',required=True)
group.add_argument(
    '--output', help='output .png filename',required=True)
group.add_argument(
    '--minimum', type=float, help='min value',required=False)
group.add_argument(
    '--maximum', type=float, help='max value',required=False)
group.add_argument(
    '--width', type=float, help='figure width',required=False,default=16)
group.add_argument(
    '--height', type=float, help='figure height',required=False,default=8)

args=parser.parse_args()
x=loadmat(args.input)['X']


#heat mapping for input matrix
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
    ax.set_ylabel('Vertices', fontsize=40)
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
    
if args.minimum is None:
    args.minimum = np.min(x)

if args.maximum is None:
    args.maximum = np.percentile(x,99.5)

cbar_spacing=np.floor((args.maximum - args.minimum)/2) #defualt to show cbar ticks at min, max, halfway

heatmapping(x,args.minimum,args.maximum+0.00001,cbar_spacing,args.width,args.height,title="Input",fname=args.output)    
#+0.00001 is for edge case where max = exact integer
