# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:07:13 2019

@author: sylvain.finot
"""

#merge noisy data
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

def mergeData(time,counts,binsize,name="",show=False):
#    counts = np.loadtxt(path) #Full histogram
#    binNumber = int(counts[0]) #Nombre de bin
#    binsize = counts[3] #Taille d'un bin
#    counts = counts[-binNumber:] #histogram
#    t = np.arange(binNumber)*binsize #echelle de temps en ns
    binNumber = len(counts)
    extrema = np.array(argrelextrema(counts*(counts>(0.9*counts.max())), np.greater)[0])
    peaks = []
    for p in extrema:
        peaks.append(np.argmax(counts[p-50:p+50])+p-50) 
        
    peaks = np.array(peaks)
    peaks = np.unique(peaks)
    peaks = peaks[(peaks>int(1/binsize)) & (peaks<(binNumber-int(1/binsize)))]
    tmin = max(0,peaks[0]-int(1/binsize))
    tmax = min(peaks[0]+int(5/binsize),binNumber)
    data0 = counts[tmin:tmax]
    rightBaseline = np.median(data0[-int(1/binsize):])
    leftBaseline  = np.median(data0[:int(1/binsize)])
    baselineError = abs(rightBaseline-leftBaseline)/min(rightBaseline,leftBaseline) > 0.20
    data = data0
    c=1
    if len(peaks) < 3: return time[tmin:tmax],data0
    if show:
        fig,ax = plt.subplots()
        ax.semilogy(data0,'.')
    for p in peaks[1:]:
        tmin = max(0,p-int(1/binsize))
        tmax = min(p+int(5/binsize),binNumber)
        print(max(np.abs(counts[tmin:tmax]-data0)/data0))
        diff= np.abs(counts[tmin:tmax]-data0)/data0
        if baselineError:
            if np.any(diff>0.20):
                continue
        if show:ax.semilogy(counts[tmin:tmax],'.')
        data += counts[tmin:tmax]
        c+=1
#        data0 +=data
        
    data = data/c
    if show:
        ax.semilogy(data0,'.')
        ax.set_title(name)
    return time[tmin:tmax],data