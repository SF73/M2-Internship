# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:06:41 2019

@author: Sylvain
"""

from expfit import process_fromFile
import os
import matplotlib.pyplot as plt
import numpy as np
import itertools
import pickle

#plt.rc('font', family='serif',size=14)
#plt.rc('text', usetex=False)
#plt.rc('xtick', labelsize=14)
#plt.rc('ytick', labelsize=14)
#plt.rc('axes', labelsize=18)
plt.style.use('Rapport')
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

Tmin = 0
Tmax =10
fig, ax = plt.subplots()
fig.patch.set_alpha(0)
markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
#files=getListOfFiles(r"C:\Users\sylvain.finot\Documents\data")
files=getListOfFiles(r"F:\data")
yerr = 8
# =============================================================================
# T2594
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2594"in x) & (not (("Al" in x) or ("Ag" in x))))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    if not((Temp<Tmax)&(Temp>Tmin)):continue
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(tau)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
#bx.plot(np.repeat(10,len(Tau[T>290])),Tau[T>290],m,alpha=0.6,label='T2594')
ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2594')
#ax.errorbar(Tau[T>290],T10[T>290],yerr=yerr,fmt=m,alpha=0.6,label='T2594')
# =============================================================================
## =============================================================================
## T2597
#mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2597"in x))]
#T = list()
#Tau = list()
#T10 = list()
#for p in mask:
#    idx = p.find('K')
#    Temp = int(p[idx-3:idx])
#    if not((Temp<Tmax)&(Temp>Tmin)):continue
#    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
#    if ((1-R)>1e-3):continue
#    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
#    T.append(Temp)
#    Tau.append(tau)
#    T10.append(t10)
#T = np.array(T)
#Tau = np.array(Tau)*1e3
#T10 = np.array(T10)*1e3
#m = next(markers)
##bx.plot(np.repeat(20,len(Tau[T>290])),Tau[T>290],m,alpha=0.6,label='T2597')
##ax.errorbar(Tau[T>290],T10[T>290],yerr=yerr,fmt=m,alpha=0.6,label='T2597')
#ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2597')
## =============================================================================
##
## =============================================================================
## T2601
#mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2601"in x)& (not ("375nm" in x)))]
#T = list()
#Tau = list()
#T10 = list()
#for p in mask:
#    idx = p.find('K')
#    Temp = int(p[idx-3:idx])
#    if not((Temp<Tmax)&(Temp>Tmin)):continue
#    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
#    if ((1-R)>1e-3):continue
#    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
#    T.append(Temp)
#    Tau.append(tau)
#    T10.append(t10)
#T = np.array(T)
#Tau = np.array(Tau)*1e3
#T10 = np.array(T10)*1e3
#m = next(markers)
##bx.plot(,Tau[T>290],m,alpha=0.6,label='T2601')
##ax.errorbar(Tau[T>290],T10[T>290],yerr=yerr,fmt=m,alpha=0.6,label='T2601')
#ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2601')
## =============================================================================
# =============================================================================
# T2594 Al
#mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2594"in x) & ("Al"in x)& ("450"in x))]
#T = list()
#Tau = list()
#T10 = list()
#for p in mask:
#    idx = p.find('K')
#    Temp = int(p[idx-3:idx])
#    if not((T<Tmax)&(T>Tmin)):continue
#    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
#    if ((1-R)>1e-3):continue
#    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
#    T.append(Temp)
#    Tau.append(tau)
#    T10.append(t10)
#T = np.array(T)
#Tau = np.array(Tau)*1e3
#T10 = np.array(T10)*1e3
#m = next(markers)
##bx.plot(,Tau[T>290],m,alpha=0.6,label='T2601')
##ax.errorbar(Tau[T>290],T10[T>290],yerr=yerr,fmt=m,alpha=0.6,label='T2601')
#ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2594 Al')
# =============================================================================
# =============================================================================
# T2594 Al
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("Al2"in x)& ("450"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    if not((Temp<Tmax)&(Temp>Tmin)):continue
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(tau)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
#bx.plot(,Tau[T>290],m,alpha=0.6,label='T2601')
#ax.errorbar(Tau[T>290],T10[T>290],yerr=yerr,fmt=m,alpha=0.6,label='T2601')
ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2594 Al')
# =============================================================================
## =============================================================================
# T2594 Ag
mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("Ag"in x)& ("T2594"in x))]
T = list()
Tau = list()
T10 = list()
for p in mask:
    idx = p.find('K')
    Temp = int(p[idx-3:idx])
    if not((Temp<Tmax)&(Temp>Tmin)):continue
    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
    if ((1-R)>1e-3):continue
    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
    T.append(Temp)
    Tau.append(tau)
    T10.append(t10)
T = np.array(T)
Tau = np.array(Tau)*1e3
T10 = np.array(T10)*1e3
m = next(markers)
ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2594 Ag')
## =============================================================================


## =============================================================================
## T2628
#mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2628"in x))]
#T = list()
#Tau = list()
#T10 = list()
#for p in mask:
#    idx = p.find('K')
#    Temp = int(p[idx-3:idx])
#    if not((Temp<Tmax)&(Temp>Tmin)):continue
#    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=True,autoclose=True,merge=True)
#    if ((1-R)>1e-2):continue
#    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
#    T.append(Temp)
#    Tau.append(taueff)
#    T10.append(t10)
#T = np.array(T)
#Tau = np.array(Tau)*1e3
#T10 = np.array(T10)*1e3
#m = next(markers)
#ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2628')
## =============================================================================
## =============================================================================
## T2631
#mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2631"in x))]
#T = list()
#Tau = list()
#T10 = list()
#for p in mask:
#    idx = p.find('K')
#    Temp = int(p[idx-3:idx])
#    if not((Temp<Tmax)&(Temp>Tmin)):continue
#    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
#    if ((1-R)>1e-2):continue
#    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
#    T.append(Temp)
#    Tau.append(taueff)
#    T10.append(t10)
#T = np.array(T)
#Tau = np.array(Tau)*1e3
#T10 = np.array(T10)*1e3
#m = next(markers)
#ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2631')
## =============================================================================
## =============================================================================
## T2629
#mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2629"in x))]
#T = list()
#Tau = list()
#T10 = list()
#for p in mask:
#    idx = p.find('K')
#    Temp = int(p[idx-3:idx])
#    if not((Temp<Tmax)&(Temp>Tmin)):continue
#    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
#    if ((1-R)>1e-2):continue
#    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
#    T.append(Temp)
#    Tau.append(taueff)
#    T10.append(t10)
#T = np.array(T)
#Tau = np.array(Tau)*1e3
#T10 = np.array(T10)*1e3
#m = next(markers)
#ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2629')
## =============================================================================
## =============================================================================
## T2630
#mask = [x for x in files if ((("TRCL" in x) & (x.endswith(".dat")))& ("T2630"in x))]
#T = list()
#Tau = list()
#T10 = list()
#for p in mask:
#    idx = p.find('K')
#    Temp = int(p[idx-3:idx])
#    if not((Temp<Tmax)&(Temp>Tmin)):continue
#    A,tau,t10,A1,A2,tau1,tau2,R = process_fromFile(p,save=False,autoclose=True,merge=True)
#    if ((1-R)>1e-2):continue
#    taueff = -(A1*tau1+A2*tau2)/(A1+A2)
#    T.append(Temp)
#    Tau.append(taueff)
#    T10.append(t10)
#T = np.array(T)
#Tau = np.array(Tau)*1e3
#T10 = np.array(T10)*1e3
#m = next(markers)
#ax.plot(Tau[(T<Tmax)&(T>Tmin)],T10[(T<Tmax)&(T>Tmin)],m,alpha=0.6,label='T2630')
## =============================================================================


ax.set_xlabel(r"$\tau$ (ps)")
ax.set_ylabel(r"$\tau_{10}$ (ps)")
ax.legend()
plt.show()
