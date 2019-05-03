# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:19:50 2019

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


fig, (ax,bx) = plt.subplots(1,2,sharex=True)
markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
files=getListOfFiles(r"C:\Users\sylvain.finot\Documents\data")
# =============================================================================
# T2594
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2594"in x) & (not ("Rampe"in x)) & (not("Al" in x)) & (not("Ag" in x)))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2594')
ax.plot(T,Tau,m,alpha=0.6,label='T2594')
# =============================================================================
# =============================================================================
# T2594
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2594"in x) & ("Rampe"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2594 same wire')
ax.plot(T,Tau,m,alpha=0.6,label='T2594 same wire')
# =============================================================================
# =============================================================================
# T2597
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2597"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2597')
ax.plot(T,Tau,m,alpha=0.6,label='T2597')
# =============================================================================
# =============================================================================
# T2597
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2601"in x) & (not ("375nm" in x)))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
ax.plot(T,Tau,m,alpha=0.6,label='T2601')
bx.plot(T,T10,m,alpha=0.6,label='T2601')
# =============================================================================


# T2594 Al
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2594"in x) & (not ("Rampe"in x)) & ("Al" in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=False,merge=True)
    if ((1-R)>1e-3):
        plt.close()
        continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
bx.plot(T,T10,m,alpha=0.6,label='T2594 Al')
ax.plot(T,Tau,m,alpha=0.6,label='T2594 Al')

# T2594 Ag
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2594"in x) & (not ("Rampe"in x)) & ("Ag" in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=False,merge=True)
    if ((1-R)>1e-3):
        plt.close()
        continue
    plt.close()
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(taueff)
    T10.append(t10)

T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
ax.plot(T,Tau,m,alpha=0.6,label='T2594 Ag')
bx.plot(T,T10,m,alpha=0.6,label='T2594 Ag')



ax.set_xlabel("Temperature (K)")
ax.set_ylabel(r"$\tau$ (ps)")
ax.legend()

bx.set_xlabel("Temperature (K)")
bx.set_ylabel(r"$\tau_{10}$ (ps)")
bx.legend()