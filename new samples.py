# -*- coding: utf-8 -*-
"""
Created on Thu May  2 14:39:07 2019

@author: sylvain.finot
"""

from expfit import process_fromFile
import os
import matplotlib.pyplot as plt
import numpy as np
import itertools
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


fig, (ax,bx) = plt.subplots(1,2)
markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
#files=getListOfFiles(r"C:\Users\sylvain.finot\Documents\data")
files=getListOfFiles(r"F:\data")
# =============================================================================
# T2628
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2628"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=True,autoclose=True,merge=True)
    if ((1-R)>1e-2):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2628')
ax.plot(T,Tau,m,alpha=0.6,label='T2628')
# =============================================================================
# =============================================================================
# T2629
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2629"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=True,autoclose=True,merge=True)
    if ((1-R)>1e-2):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2629')
ax.plot(T,Tau,m,alpha=0.6,label='T2629')
# =============================================================================
# =============================================================================
# T2630
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2630"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=True,autoclose=True,merge=True)
    if ((1-R)>1e-2):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2630')
ax.plot(T,Tau,m,alpha=0.6,label='T2630')
# =============================================================================
# =============================================================================
# T2630
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2631"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=True,autoclose=True,merge=True)
    if ((1-R)>1e-2):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2631')
ax.plot(T,Tau,m,alpha=0.6,label='T2631')
# =============================================================================

ax.set_xlabel("Temperature (K)")
ax.set_ylabel(r"$\tau$ (ps)")
ax.legend()

bx.set_xlabel("Temperature (K)")
bx.set_ylabel(r"$\tau_{10}$ (ps)")
bx.legend()