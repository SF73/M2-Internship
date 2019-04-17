# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:19:50 2019

@author: sylvain.finot
"""

from expfit import process
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


fig, ax = plt.subplots()
markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
files=getListOfFiles(r"C:\Users\sylvain.finot\Documents\data")
# =============================================================================
# T2594
mask = [x for x in files if ((x.endswith("TRCL.dat"))& ("T2594"in x) & (not ("Rampe"in x)))]
T = list()
Tau = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,A1,A2,tau1,tau2,R = process(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)

T = np.array(T)
Tau = np.array(Tau)*1e3
ax.plot(T,Tau,next(markers),alpha=0.6,label='T2594')
# =============================================================================
# =============================================================================
# T2594
mask = [x for x in files if ((x.endswith("TRCL.dat"))& ("T2594"in x) & ("Rampe"in x))]
T = list()
Tau = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,A1,A2,tau1,tau2,R = process(p,save=False,autoclose=True,merge=True)
#    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)

T = np.array(T)
Tau = np.array(Tau)*1e3
ax.plot(T,Tau,next(markers),alpha=0.6,label='T2594 same wire')
# =============================================================================
# =============================================================================
# T2597
mask = [x for x in files if ((x.endswith("TRCL.dat"))& ("T2597"in x))]
T = list()
Tau = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,A1,A2,tau1,tau2,R = process(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)

T = np.array(T)
Tau = np.array(Tau)*1e3
ax.plot(T,Tau,next(markers),alpha=0.6,label='T2597')
# =============================================================================

# =============================================================================
# T2597
mask = [x for x in files if ((x.endswith("TRCL.dat"))& ("T2601"in x) & (not ("375nm" in x)))]
T = list()
Tau = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,A1,A2,tau1,tau2,R = process(p,save=False,autoclose=False,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)

T = np.array(T)
Tau = np.array(Tau)*1e3
ax.plot(T,Tau,next(markers),alpha=0.6,label='T2601')
# =============================================================================
ax.set_xlabel("Temperature (K)")
ax.set_ylabel(r"$\tau$ (ps)")
ax.legend()