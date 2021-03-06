# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 10:30:28 2019

@author: Sylvain
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from FullReport import make_linescan
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

    
#path = r"D:\M2 Internship\data\2019-07-10 - T2455 Top - 300K\HYP3-T2455-300K-Vacc5kV-spot4-zoom2000x-gr600-slit0-2-t050ms-cw420nm\Hyp.dat"
path = r"D:\M2 Internship\data\2019-07-10 - T2455 Top - 300K\HYP2-T2455-300K-Vacc5kV-spot4-zoom16000x-gr600-slit0-2-t005ms-cw420nm\Hyp.dat"
hyp, wavelenght, X, Y = make_linescan(path,save=False,autoclose=False,deadPixeltol=200,Linescan=False)
test = np.sum(hyp[:,:,1325:1847],axis=-1)
droite = test[:,45:]
gauche = test[:,:45]
fig,ax = plt.subplots()
profil = np.sum(test,axis=0)
ax.plot(X,profil/profil.max())
print("I_Al/I_Bare : %f" %(gauche.mean()/droite.mean()))
fig2,bx = plt.subplots()
bx.hist(test.flatten(), bins='auto')


#20 nm !!!!
path=r'D:/M2 Internship/data/2019-05-16 - T2455 Al - 300K/FULLSCREEN HYP2-T2455-300K-Vacc5kV-spot5-zoom8000x-gr600-slit0-2-t005ms-cw440nm/Hyp.dat'
hyp, wavelenght, X, Y = make_linescan(path,save=False,autoclose=False,deadPixeltol=200,Linescan=False)
test = np.sum(hyp[:,:,534:1318],axis=-1)
#droite = test[:,45:]
#gauche = test[:,:45]
fig,ax = plt.subplots()
ax.hist(test[33:81,56:79].flatten(), bins=200)

fig,ax = plt.subplots()
ax.pcolormesh(X[56:79],Y[33:81],test[33:81,56:79])
#profil = np.sum(test,axis=0)
#ax.plot(profil/profil.max())

#test = hyp[65:74,52:58,534:1318].sum(axis=-1)
#haut = test[0:4]
#bas = test[5:10]
#haut.mean()/bas.mean()
#print("I_Al/I_Bare : %f" %(gauche.mean()/droite.mean()))
#fig2,bx = plt.subplots()
#bx.hist(test.flatten(), bins='auto')







