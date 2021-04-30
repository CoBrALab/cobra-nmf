import sys, os
import glob
import pandas as pd
import numpy as np
import scipy
import scipy.io
from scipy.io import savemat, loadmat
from sklearn.metrics.pairwise import cosine_similarity
import errno
import argparse

parser=argparse.ArgumentParser(
    description='''This script computed stability metrics (spatial similarity,
    recon error) for a specified number of components and stores results in a .csv file''')

group = parser.add_argument_group(title="Execution options")

group.add_argument('--n_folds', help='number of folds', type=int, default=10)

group.add_argument(
    "--stability_results_dir",help="parent dir of stability nmf outputs", required=True)

group.add_argument(
    "--k",help="granularity to assess", type=int, required=True)


args=parser.parse_args()

n_splits = args.n_folds
g = args.k

out_dir = "stability_correlations/" 
stab_dir = args.stability_results_dir

if not os.path.exists(out_dir): #make output directory
    try:
        os.makedirs(out_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else:
            raise

cols = ["Granularity","Iteration","Corr_mean","Corr_median","Corr_std","Recon_errorA","Recon_errorB"]
df = pd.DataFrame(columns = cols)
for i in range(0,n_splits):
    
    #load split input, get W mx for each
    fname = stab_dir + "/k" + str(g) + "/a_" + str(i) + ".mat" 
    resA = scipy.io.loadmat(fname)
    Wa = resA['W']
    ea = resA['recon']
        
    fname = stab_dir + "/k" + str(g) + "/b_" + str(i) + ".mat" 
    resB = scipy.io.loadmat(fname)
    Wb = resB['W']
    eb = resB['recon']

    #assess correlation of identified parcel component scores - which parcels vary together?
    c_Wa = cosine_similarity(Wa)
    c_Wb = cosine_similarity(Wb)
        
    corr = np.zeros((1,np.shape(c_Wa)[0]))

    for parcel in range(0,np.shape(c_Wa)[0]):
        corr[0,parcel] = np.corrcoef(c_Wa[parcel,:],c_Wb[parcel,:])[0,1]

    df_curr = pd.DataFrame(data = [[g, i+1, np.mean(corr),np.median(corr),np.std(corr),ea,eb]], columns = cols)

    print(i, np.mean(corr))
    df = df.append(df_curr)
    del Wa,Wb,ea,eb,resA,resB
    
fname = out_dir + "stability_corr_k" + str(g) + ".csv"
print(fname)
df.to_csv(fname)
del df, df_curr

