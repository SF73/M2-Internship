# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 10:30:28 2019

@author: Sylvain
"""

import matplotlib.pyplot as plt
import numpy as np
import os

dirpath = r"D:\M2 Internship\data\2019-07-10 - T2455 Top - 300K"

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
def process_all(path):
    files=getListOfFiles(path)
    mask = [x for x in files if (("TRCL" in x) & ("T2455 Top - 300K" in x)) & x.endswith(".dat")& (not(x.endswith("Hyp.dat")))]
    print(mask)
    Al = []
    GaN = []
#    other = np.loadtxt(r"F:\data\2019-03-08 - T2455 - withUL\TRCL-440nm.dat")[4:]
    for i in range(len(mask)):
        try:
#            process_fromFile(mask[i],save=True,autoclose=True)
            p = mask[i]
            if ('bare'in mask[i]):
                GaN.append(np.loadtxt(p)[3:])
            elif ('etrange' in mask[i]):
                Al.append(np.loadtxt(p)[3:])
            print(100*i/len(mask))
        except Exception as e:
            print(e)
            pass
    Al = np.asarray(Al)
    t=np.arange(len(Al[0])-1)
    try:
        for i in Al[:-1]:
            plt.plot(t*i[0]-100,i[1:]/np.mean(i[3100:3300]),c='C0',alpha=0.6)
        plt.plot(t*Al[-1][0]-100,Al[-1][1:]/np.mean(Al[-1][3100:3300]),c='C0',alpha=1,label='On Al pad')
        for i in GaN[:-1]:
            plt.plot(t*i[0]-100,i[1:]/np.mean(i[3100:3300]),c='C1',alpha=0.6)
        plt.plot(t*GaN[-1][0]-100,GaN[-1][1:]/np.mean(GaN[-1][3100:3300]),c='C1',alpha=1,label='Next to Al pad')
    except Exception as e:
        print(e)
        pass

#    plt.plot(t,np.exp(-(t-19.55)/15.6))
#    plt.plot(t,np.exp(-(t-19.55)/181))
#    plt.plot(t/2,other/(0.94*max(other)))
    plt.xlim((0,550))
    plt.ylim((1e-2,1.1))
    plt.gca().set_yscale('log')
    plt.legend()
    plt.show()
process_all(dirpath)