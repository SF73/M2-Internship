# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 09:50:06 2019

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


markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
fig,ax = plt.subplots()

files=getListOfFiles(r"F:\data\2019-03-21 - T2594 - 005K\Fil 3")
mask = [x for x in files if (("TRCL" in x) & (x.endswith(".dat")))]
taus=list()
A1s=list()
A2s=list()
tau1s=list()
tau2s=list()
tau10=[]
L = list()
for p in mask:
    try:
        name = os.path.basename(os.path.dirname(p))
        lenght = float(name[5:name.find("um")])
        A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=False)
        L.append(lenght) 
        taus.append(-tau)
        A1s.append(A1)
        A2s.append(A2)
        tau10.append(t10)
        tau1s.append(-tau1)
        tau2s.append(-tau2)
    except Exception as e:
        print(e)
        print(lenght)
taus=np.array(taus)*1e3
tau10 = np.array(tau10)*1e3
tau1s=np.array(tau1s)*1e3
tau2s=np.array(tau2s)*1e3
A1s=np.array(A1s)
A2s=np.array(A2s)
taueffs = (A1s*tau1s+A2s*tau2s)/(A1s+A2s)
L = np.array(L)+14.42
m=next(markers)
#ax.plot(L,taueffs,m,label=r'$\tau_{eff}$')
ax.plot(L,tau10,m,label=r'$A-\tau_{10}$')
#ax.plot(L,taueffs,m,label=r'$\tau$')

files=getListOfFiles(r"F:\data\2019-03-21 - T2597 - 005K\Fil 1")
mask = [x for x in files if (("TRCL" in x) & (x.endswith(".dat")))]
taus=list()
A1s=list()
A2s=list()
tau1s=list()
tau2s=list()
tau10=[]
L = list()
for p in mask:
    try:
        name = os.path.basename(os.path.dirname(p))
        lenght = float(name[5:name.find("um")])
        A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=False)
        L.append(lenght) 
        taus.append(-tau)
        A1s.append(A1)
        A2s.append(A2)
        tau10.append(t10)
        tau1s.append(-tau1)
        tau2s.append(-tau2)
    except Exception as e:
        print(e)
        print(lenght)
taus=np.array(taus)*1e3
tau10 = np.array(tau10)*1e3
tau1s=np.array(tau1s)*1e3
tau2s=np.array(tau2s)*1e3
A1s=np.array(A1s)
A2s=np.array(A2s)
taueffs = (A1s*tau1s+A2s*tau2s)/(A1s+A2s)
L = np.array(L)+ 18.41
m=next(markers)
#ax.plot(L,taueffs,m,label=r'$\tau_{eff}$')
ax.plot(L,tau10,m,label=r'$B-\tau_{10}$')
#ax.plot(L,taueffs,m,label=r'$\tau$')

files=getListOfFiles(r"F:\data\2019-03-21 - T2601 - 005K\Fil 1")
mask = [x for x in files if (("TRCL" in x) & (x.endswith(".dat")))]
taus=list()
A1s=list()
A2s=list()
tau1s=list()
tau2s=list()
tau10=[]
L = list()
for p in mask:
    try:
        name = os.path.basename(os.path.dirname(p))
        lenght = float(name[5:name.find("um")])
        A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=False)
        L.append(lenght) 
        taus.append(-tau)
        A1s.append(A1)
        A2s.append(A2)
        tau10.append(t10)
        tau1s.append(-tau1)
        tau2s.append(-tau2)
    except Exception as e:
        print(e)
        print(lenght)
taus=np.array(taus)*1e3
tau10 = np.array(tau10)*1e3
tau1s=np.array(tau1s)*1e3
tau2s=np.array(tau2s)*1e3
A1s=np.array(A1s)
A2s=np.array(A2s)
taueffs = (A1s*tau1s+A2s*tau2s)/(A1s+A2s)
L = np.array(L)+16.12
m=next(markers)
#ax.plot(L,taueffs,m,label=r'$\tau_{eff}$')
ax.plot(L,tau10,m,label=r'$C-\tau_{10}$')
#ax.plot(L,taueffs,m,label=r'$\tau$')



#ax.plot(L,tau1s,next(markers),label=r'$\tau_1$')
#ax.plot(L,tau2s,next(markers),label=r'$\tau_2$')
#bx.plot(L,A2s/A1s,next(markers),label=r'$A_1/A_2$')
#ax.plot(L,taueffs,next(markers),label=r'$\tau_{eff}$')
ax.set_xlabel("Distance from the top (µm)")
ax.set_ylabel(r"$\tau_{10} (ps)$")
#fig.suptitle(r"TRCL linescan of T2594 - Wire 1 - 005K")
#ax.set_title(r"TRCL linescan of T2601 - Wire 1 - 005K")
#bx.set_xlabel("Distance (µm)")
#bx.set_ylabel(r"A2/A1")
ax.legend()
plt.show()
#with open("TRCL linescan of T2594 - Wire 1 - 005K_merge_bound.pickle",'wb') as f:
#    pickle.dump(fig,f)
