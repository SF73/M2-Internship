# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:07:13 2019

@author: sylvain.finot
"""

#merge noisy data
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

def mergeData(counts,binNumber,binsize,name=""):
#    counts = np.loadtxt(path) #Full histogram
#    binNumber = int(counts[0]) #Nombre de bin
#    binsize = counts[3] #Taille d'un bin
#    counts = counts[-binNumber:] #histogram
#    t = np.arange(binNumber)*binsize #echelle de temps en ns
    extrema = np.array(argrelextrema(counts*(counts>(0.8*counts.max())), np.greater)[0])
    peaks = []
    for p in extrema:
        peaks.append(np.argmax(counts[p-50:p+50])+p-50) 
        
    peaks = np.array(peaks)
    peaks = np.unique(peaks)
    tmin = max(0,peaks[0]-int(1/binsize))
    tmax = min(peaks[0]+int(5/binsize),binNumber)
    data0 = counts[tmin:tmax]
    if len(peaks) < 3: return data0
    fig,ax = plt.subplots()
    ax.semilogy(data0,'.')
    for p in peaks[1:]:
        tmin = max(0,p-int(1/binsize))
        tmax = min(p+int(5/binsize),binNumber)
        data = counts[tmin:tmax]
        ax.semilogy(data,'.')
        data0 +=data
        
    data0 = data0/len(peaks)
    ax.semilogy(data0,'.')
    ax.set_title(name)
    return data0