import pylab as plt
import numpy as np
import pandas as pd
import scipy
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
# matplotlib used plotting. Not required to run the code.
import matplotlib.pyplot as plt
import re
import sys
import glob
import crflux.models as crf
# matplotlib used plotting. Not required to run the code.
import matplotlib.pyplot as plt
import random
try:
    import cPickle as pickle
except ImportError:
    import pickle



r=0
for files in glob.glob('/n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere/SelectionCuts/run*alleventscutinfo_highTH.h5'):
    #print("starting on files",files)
    run=re.search('run(.*)all', files)
    X=str(run.group(1))

    try:
        df1=pd.read_hdf(files,key='dEdx')
        df1['run']=int(X)
    except FileNotFoundError:
        print('file not there')
        continue
    except KeyError:
        print('KeyError')
        continue


    if r==0:
        df=pd.read_hdf(files,key='dEdx')
        df['run']=int(X)
        #print(files)
        r+=1
    else:
        df=df.append(df1, ignore_index=True)
        #print(files,df[(df.perconline<.2)&(df.FiducialPass==True)])
        #print(files)



for files in glob.glob('/n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere2/SelectionCuts/run*alleventscutinfo_highTH.h5'):
    #print("starting on files",files)                                                                                                                         
    run=re.search('run(.*)all', files)
    X=str(run.group(1))

    try:
        df1=pd.read_hdf(files,key='dEdx')
        df1['run']=int(X)+10000
    except FileNotFoundError:
        print('file not there')
        continue
    except KeyError:
        print('KeyError')
        continue
    print(int(X)+10000)
    df=df.append(df1, ignore_index=True)
    
for files in glob.glob('/n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/SmallerSphere3/SelectionCuts/run*alleventscutinfo_highTH.h5'):
    #print("starting on files",files)                                                                                                                                             
    run=re.search('run(.*)all', files)
    X=str(run.group(1))

    try:
        df1=pd.read_hdf(files,key='dEdx')
        df1['run']=int(X)+20000
    except FileNotFoundError:
        print('file not there')
        continue
    except KeyError:
        print('KeyError')
        continue
    print(int(X)+20000)
    df=df.append(df1, ignore_index=True)

df['passed']=False
df[(df.dEdx<.015)&(df.linelength>73.5)&(df.perconline>.79)&(df.muenergy>=2)].passed=True




df.to_hdf('/n/holystore01/LABS/guenette_lab/Users/lrogers/PredictionStuff/data/CombinedFiles.h5',key="muons")
