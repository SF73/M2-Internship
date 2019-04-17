# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:19:21 2019

@author: sylvain.finot
"""

import os
import matplotlib.pyplot as plt
import numpy as np
from mergeNoise import mergeData
import pickle

def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles

files=getListOfFiles(r'C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2594 - Rampe')
mask = [x for x in files if ((x.endswith("TRCL.dat")))]
T = list()
fig, ax = plt.subplots()
for i in range(0,len(mask)):
    try:
        p = mask[i]
        T.append(int(os.path.basename(os.path.dirname(p))[:-1]))
        if (T[-1] in [10,270]):continue
        counts = np.loadtxt(p) #Full histogram
        binNumber = int(counts[0]) #Nombre de bin
        binsize = counts[3] #Taille d'un bin
        counts = counts[-binNumber:] #histogramh
        t = np.arange(binNumber)*binsize #echelle de temps en ns
        countmax = counts[int(1/binsize):binNumber-int(5/binsize)].argsort()[-1]+int(1/binsize)
        tmin = max(0,countmax-int(1/binsize))
        tmax = min(countmax+int(5/binsize),binNumber)
        reduced_time = t[tmin:tmax]
        reduced_time = reduced_time - min(reduced_time)
        merge=True
        if merge=='auto':
            merge = counts.max() < 5E3
        if merge==True:
            reduced_counts = mergeData(counts,binNumber,binsize)
        else:
            reduced_counts = counts[tmin:tmax]
        reduced_counts = reduced_counts / reduced_counts.max() *1e4
        ax.semilogy(reduced_time,reduced_counts,'.',label = '%dK'%T[-1])
    except:
        pass
ax.legend()
ax.set_xlabel("Time (ns)")
ax.set_ylabel(r"TRCL intensity (arb. units)")
ax.set_title('Decay time vs temperature')
#with open('data_temperature.pickle','wb') as f:
#    pickle.dump(fig,f)