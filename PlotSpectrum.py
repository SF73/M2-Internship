# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:58:12 2019

@author: sylvain.finot
"""


import scipy.constants as cst
eV_To_nm = cst.c*cst.h/cst.e*1E9
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost


def plotSpectrum(paths):
    if type(paths)==str:
        paths = [paths]
    fig=plt.figure()
    fig.patch.set_alpha(0)
    ax1=SubplotHost(fig, 111)
    fig.add_subplot(ax1)
    for p in paths:     
        data = np.loadtxt(p,skiprows=9) 
        ax1.plot(data[:,0],data[:,1])
    
    
    ax1.set_ylabel('Intensity')
    ax2=ax1.twin()
    ax2.set_xlabel('Energy (eV)')
     # ax2 is responsible for "top" axis and "right" axis
    tticks=np.round(eV_To_nm/ax1.get_xticks(),2)#np.array([20.,30.,50.,100.,300.])
    ax2.set_xticks( [ eV_To_nm/t for t in tticks ] )
    ax2.set_xticklabels(tticks)
    #ax2.axis["top"].label.set_visible(True)
    ax1.set_xlabel('Wavelength (nm)')
    ax2.set_yticks([])


def sPlot():
    fig=plt.figure()
    fig.patch.set_alpha(0)
    ax1=SubplotHost(fig, 111)
    fig.add_subplot(ax1)
    data1 = np.loadtxt(r"F:\data\2019-04-30 - T2630 - 005K\Fil 1\Spec1-300nm-690nm-005K-slit0-2-t005ms-Vacc5kV-spot5-gr600-zoom8000",skiprows=9)
    data2 = np.loadtxt(r"F:\data\2019-04-30 - T2630 - 005K\Fil 2\Spec1-300nm-690nm-005K-slit0-2-t005ms-Vacc5kV-spot5-gr600-zoom8000",skiprows=9)
    data3 = np.loadtxt(r"F:\data\2019-04-30 - T2630 - 005K\Fil 3\Spec1-300nm-690nm-005K-slit0-2-t005ms-Vacc5kV-spot5-gr600-zoom8000",skiprows=9)
    X = data1[:,0]
    
    Y = np.zeros((3,data1.shape[0]))
    Y[0] = data1[:,1]
    Y[1] = data2[:,1]
    Y[2] = data3[:,1]
    Y = Y.T
    Baseline = np.zeros((1,3))
    Baseline = np.mean(Y[0:150],axis=0)
    Data = np.abs(Y - Baseline)
    #Data = Data/np.max(Data[1288:3888],axis=0)
    
    Xs = np.tile(X,(3,1))
    ax1.plot(Xs.T,Data)
    ax1.set_ylabel('Intensity')
    ax2=ax1.twin() # ax2 is responsible for "top" axis and "right" axis
    ax2.set_xlabel('Energy (eV)')
    
    tticks=np.round(eV_To_nm/ax1.get_xticks(),2)#np.array([20.,30.,50.,100.,300.])
    ax2.set_xticks( [ eV_To_nm/t for t in tticks ] )
    ax2.set_xticklabels(tticks)
    #ax2.axis["top"].label.set_visible(True)
    ax1.set_xlabel('Wavelength (nm)')
    ax2.set_yticks([])

def main():
    path = input("Enter the path of your file: ")
    path=path.replace('"','')
    path=path.replace("'",'')
#    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 5K/Fil3/TRCL-cw455nm/TRCL.dat'
    plotSpectrum(path)
    
#if __name__ == '__main__':
#    main()