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



def main():
    path = input("Enter the path of your file: ")
    path=path.replace('"','')
    path=path.replace("'",'')
#    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 5K/Fil3/TRCL-cw455nm/TRCL.dat'
    plotSpectrum(path)
    
if __name__ == '__main__':
    main()