# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 09:50:38 2019

@author: sylvain.finot
"""

from expfit import process
import os
import matplotlib.pyplot as plt
import numpy as np
import itertools
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

#def process_all():
files=getListOfFiles(r'C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2594 - Rampe')
mask = [x for x in files if (("TRCL" in x) & (x.endswith(".dat")))]
taus=list()
A1s=list()
A2s=list()
tau1s=list()
tau2s=list()
T = list()
for p in mask:
    try:
        A,tau,A1,A2,tau1,tau2 = process(p,save=False,autoclose=False,merge=True)
        T.append(int(os.path.basename(os.path.dirname(p))[:-1])) 
        taus.append(-tau)
        A1s.append(A1)
        A2s.append(A2)
        tau1s.append(-tau1)
        tau2s.append(-tau2)
    except:
        pass
taus=np.array(taus)*1e3
tau1s=np.array(tau1s)*1e3
tau2s=np.array(tau2s)*1e3
A1s=np.array(A1s)
A2s=np.array(A2s)
taueffs = (A1s*tau1s+A2s*tau2s)/(A1s+A2s)
markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
fig,(ax,bx) = plt.subplots(1,2)
ax.plot(T,taus,next(markers),label=r'$\tau$')
ax.plot(T,tau1s,next(markers),label=r'$\tau_1$')
ax.plot(T,tau2s,next(markers),label=r'$\tau_2$')
bx.plot(T,A2s/A1s,next(markers),label=r'$A_1/A_2$')
ax.plot(T,taueffs,next(markers),label=r'$\tau_{eff}$')
ax.set_xlabel("Temperature (K)")
ax.set_ylabel(r"$\tau$ (ps)")
fig.suptitle(r"TRCL of T2594 (w/o) from 10K to RT")
bx.set_xlabel("Temperature (K)")
bx.set_ylabel(r"A2/A1")
ax.legend()
