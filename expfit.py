# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 12:50:23 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import scipy.special as sse
from mergeNoise import mergeData
from stats import R2,R2adj,rChi2
from models import *


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



def process_fromArray(t,counts,save=False,autoclose=False,merge=True,fig=None,show=False):
    binsize = t[1]-t[0]
    binNumber = len(counts)
    name = ''
    if show:
        if fig is None:
            fig, ax = plt.subplots()
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
            ax.set_xlabel("t (ns)")
            ax.set_ylabel("Intensity (arb. unit)")
            ax.set_title("Carrier lifetime : %s"%name)
        else:
            ax = fig.gca()
    countmax = counts[int(1/binsize):binNumber-int(5/binsize)].argsort()[-1]+int(1/binsize)
    tmin = max(0,countmax-int(1/binsize))
    tmax = min(countmax+int(5/binsize),binNumber)
    if merge=='auto':
        merge = counts.max() < 5E3
    if merge==True:
        reduced_time,reduced_counts = mergeData(t,counts,binsize,name)
    else:
        reduced_time = t[tmin:tmax]
        reduced_counts = counts[tmin:tmax]
    reduced_time = reduced_time - min(reduced_time)
    rightBaseline = np.median(reduced_counts[-int(1/binsize):])
    leftBaseline  = np.median(reduced_counts[:int(1/binsize)])
    baselineError = abs(rightBaseline-leftBaseline)/rightBaseline > 0.50
#    ax.plot(reduced_time,reduced_counts,'.',c='k',label='data')
    baseline = rightBaseline
    
    #calcul de la limite de droite pour fit
    leftl = reduced_counts.argmax()+2
    rightl = 0
    c=1e-2#np.exp(-1)
    while(rightl<leftl or np.isnan(rightl)):        
        threshold = (max(reduced_counts)-baseline)*c #(max(reduced_counts)-baseline)*np.exp(-3)
        mask = np.convolve(np.sign(reduced_counts-threshold),[-1,1],'same') #detect le chgt de sign de reduced-threshold
        mask[0] = 0
        rightl = np.argmax(mask)
        c += 0.005
        
    t0 = reduced_time[leftl]#temps correspondant au max
    t10 = reduced_time[rightl] - t0
    if show:
        ax.plot(reduced_time,reduced_counts,'.',c='k',label='data | t%.2d = %.2e'%(np.round((c-0.005)*100),t10))
        ax.hlines(reduced_counts[rightl],t0,t0+t10)
    SNR = max(reduced_counts)/baseline
    if show:
        print("SNR : %.4f"%SNR)
    #Fit exponential decay
        print("------------------simple decay------------------")
    #seuleement decay (pas de montée)
    fit_time = reduced_time[leftl:rightl]
    fit_count = reduced_counts[leftl:rightl]
    
    popt,pcov= fit(fit_time,fit_count,baseline)
    A_lin=np.exp(popt[0])
    K_lin = popt[1]
    p_sigma = np.sqrt(np.diag(pcov))
    R = R2(fit_count,decay_func(fit_time,t0,A_lin,K_lin)+baseline)
    Radj = R2adj(len(fit_count),2,fit_count,decay_func(fit_time,t0,A_lin,K_lin)+baseline)
    rChi = rChi2(len(fit_count),2,fit_count,decay_func(fit_time,t0,A_lin,K_lin)+baseline)
#    ax.plot(fit_time,decay_func(fit_time,t0,A_lin,K_lin)+baseline,label=r'decay fit $R^2 =$%.4f%s$\tau_{eff} =$ %.2e ns'%(R,'\n',-1/K_lin))
    if show:
        print("R2 : %.4f"%R)
        print("R2adj : %.4f"%Radj)
        print("rChi2 : %.4f"%rChi)
    
    fit_time = reduced_time[leftl:-300]
    fit_count = reduced_counts[leftl:-300]
    
    #calcul de la limite de gauche
    c=np.exp(-3)
    while(np.isnan(leftl)):        
        threshold = (max(reduced_counts)-leftBaseline)*c #(max(reduced_counts)-baseline)*np.exp(-3)
        mask = np.convolve(np.sign(reduced_counts-threshold),[1,-1],'same') #detect le chgt de sign de reduced-threshold
        mask[0] = 0
        leftl = np.argmax(mask)
        c += 0.01
    
    if show:print("------------------model1-----------------------")
    #fit simple exp convoluée avec heaviside
#    if not baselineError:
    fit_time = reduced_time[leftl-50:]
    fit_count = reduced_counts[leftl-50:]
    init = [A_lin,K_lin,0.02,t0]
    popt,pcov= model_fit(fit_time,fit_count-baseline,init)
    if show:print(popt)
    A=popt[0]
    K = popt[1]
    sig = popt[2]
    t0 = popt[3]

    R = R2(fit_count,model_func(fit_time,*popt)+baseline)
    Radj = R2adj(len(fit_count),3,fit_count,model_func(fit_time,*popt)+baseline)
    rChi = rChi2(len(fit_count),3,fit_count,model_func(fit_time,*popt)+baseline)
#    ax.plot(fit_time,decay_func(fit_time,t0,A_lin,K_lin)+baseline,label=r'decay fit $R^2 =$%.4f%s$\tau_{eff} =$ %.2e ns'%(R,'\n',-1/K_lin))
    if show:
        print("R2 : %.4f"%R)
        print("R2adj : %.4f"%Radj)
        print("rChi2 : %.4f"%rChi)
        ax.plot(fit_time,model_func(fit_time,*popt)+baseline,c='orange',label=r'simple_fit $1-R^2 =$ %.2e %s $\tau_{eff} =$ %.2e ns %s  $\sigma=$%.2e ns'%((1-R),'\n',-1/K,'\n',sig))

    
    
        print("------------------model2-----------------------")
    init = [A,K,1,1,sig,t0]
    popt,pcov= model2_fit(fit_time,fit_count-baseline,init)
    if show:print(popt)
    A1 = popt[0]
    K1 = popt[1]
    A2 = popt[2]
    K2 = popt[3]
    sig = popt[4]
    t0 = popt[5]
    R = R2(fit_count,model2_func(fit_time,*popt)+baseline)
    Radj = R2adj(len(fit_count),3,fit_count,model2_func(fit_time,*popt)+baseline)
    rChi = rChi2(len(fit_count),3,fit_count,model2_func(fit_time,*popt)+baseline)
    taueff = (A1*(-1/K1)+A2*(-1/K2))/(A1+A2)
    if show:
        print("R2 : %.4f"%R)
        print("R2adj : %.4f"%Radj)
        print("rChi2 : %.4f"%rChi)
        print("A1/A2 : %.4e"%(A1/A2))
    
        print("taueff : %.4e ns"%taueff)
        ax.plot(fit_time,model2_func(fit_time,*popt)+baseline,c='red',label=r'double_fit $1-R^2 =$ %.2e %s$A_{1} =$ %.2e %s$A_{2} =$ %.2e %s$\tau_{1} =$ %.2e ns %s$\tau_{2} =$ %.2e ns %s$\sigma =$ %.2e ns %s$\tau_{eff} =$ %.2e ns'%((1-R),'\n',A1,'\n',A2,'\n',-1/K1,'\n',-1/K2,'\n',sig,'\n',taueff))

        ax.legend()
        fig.tight_layout()
        ax.set_yscale('log')
        if autoclose==True:
            plt.close(fig)
    return A,1/K,A1,A2,1/K1,1/K2,R
def process_fromFile(path,save=False,autoclose=False,merge=True,fig=None):
#    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 005K/Fil3/TRCL-cw455nm/TRCL.dat'
#    path = r"C:\Users\sylvain.finot\Documents\data\Probleme TRCL\TRCL_baseline.dat"
#    path = r"C:\Users\sylvain.finot\Documents\data\Probleme TRCL\TRCL_sync.dat"
#    path=r"C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2594 - Rampe\300K\TRCL.dat"
#    path=r'C:/Users/sylvain.finot/Documents/data/2019-03-08 - T2597 - 300K/Fil2/TRCL-L12-5um-cw440nm/TRCL.dat'
#    path = r"C:\Users\sylvain.finot\Documents\data\2019-03-21 - T2594 - 005K\Fil 3\TRCL 4_-04.42um"
    name = path[path.find('2019'):]
    if len(name)>100:name=''
    if fig is None:
        fig, ax = plt.subplots()
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        ax.set_xlabel("t (ns)")
        ax.set_ylabel("Intensity (arb. unit)")
        ax.set_title("Carrier lifetime : %s"%name)
    else:
        ax = fig.gca()
    
    
    counts = np.loadtxt(path) #Full histogram
    binNumber = int(counts[0]) #Nombre de bin
    binsize = counts[3] #Taille d'un bin
    counts = counts[-binNumber:] #histogramh
    t = np.arange(binNumber)*binsize #echelle de temps en ns
    countmax = counts[int(1/binsize):binNumber-int(5/binsize)].argsort()[-1]+int(1/binsize)
    tmin = max(0,countmax-int(1/binsize))
    tmax = min(countmax+int(5/binsize),binNumber)
    if merge=='auto':
        merge = counts.max() < 5E3
    if merge==True:
        reduced_time,reduced_counts = mergeData(t,counts,binsize,name)
    else:
        reduced_time = t[tmin:tmax]
        reduced_counts = counts[tmin:tmax]
    reduced_time = reduced_time - min(reduced_time)
    rightBaseline = np.median(reduced_counts[-int(1/binsize):])
    leftBaseline  = np.median(reduced_counts[:int(1/binsize)])
    baselineError = abs(rightBaseline-leftBaseline)/rightBaseline > 0.50
#    ax.plot(reduced_time,reduced_counts,'.',c='k',label='data')
    baseline = rightBaseline
    
    #calcul de la limite de droite pour fit
    leftl = reduced_counts.argmax()+2
    rightl = 0
    c=1e-2#np.exp(-1)
    while(rightl<leftl or np.isnan(rightl)):        
        threshold = (max(reduced_counts)-baseline)*c #(max(reduced_counts)-baseline)*np.exp(-3)
        mask = np.convolve(np.sign(reduced_counts-threshold),[-1,1],'same') #detect le chgt de sign de reduced-threshold
        mask[0] = 0
        rightl = np.argmax(mask)
        c += 0.005
        
    t0 = reduced_time[leftl]#temps correspondant au max
    t10 = reduced_time[rightl] - t0
    ax.plot(reduced_time,reduced_counts,'.',c='k',label='data | t%.2d = %.2e'%(np.round((c-0.005)*100),t10))
    SNR = max(reduced_counts)/baseline
    print("SNR : %.4f"%SNR)
    #Fit exponential decay
    print("------------------simple decay------------------")
    #seuleement decay (pas de montée)
    fit_time = reduced_time[leftl:rightl]
    fit_count = reduced_counts[leftl:rightl]
    
    popt,pcov= fit(fit_time,fit_count,baseline)
    A_lin=np.exp(popt[0])
    K_lin = popt[1]
    p_sigma = np.sqrt(np.diag(pcov))
    R = R2(fit_count,decay_func(fit_time,t0,A_lin,K_lin)+baseline)
    Radj = R2adj(len(fit_count),2,fit_count,decay_func(fit_time,t0,A_lin,K_lin)+baseline)
    rChi = rChi2(len(fit_count),2,fit_count,decay_func(fit_time,t0,A_lin,K_lin)+baseline)
#    ax.plot(fit_time,decay_func(fit_time,t0,A_lin,K_lin)+baseline,label=r'decay fit $R^2 =$%.4f%s$\tau_{eff} =$ %.2e ns'%(R,'\n',-1/K_lin))
    print("R2 : %.4f"%R)
    print("R2adj : %.4f"%Radj)
    print("rChi2 : %.4f"%rChi)
    
    fit_time = reduced_time[leftl:-300]
    fit_count = reduced_counts[leftl:-300]
    
    #calcul de la limite de gauche
    c=np.exp(-3)
    while(np.isnan(leftl)):        
        threshold = (max(reduced_counts)-leftBaseline)*c #(max(reduced_counts)-baseline)*np.exp(-3)
        mask = np.convolve(np.sign(reduced_counts-threshold),[1,-1],'same') #detect le chgt de sign de reduced-threshold
        mask[0] = 0
        leftl = np.argmax(mask)
        c += 0.01
    
    print("------------------model1-----------------------")
    #fit simple exp convoluée avec heaviside
#    if not baselineError:
    fit_time = reduced_time[leftl-50:]
    fit_count = reduced_counts[leftl-50:]
    init = [A_lin,K_lin,0.02,t0]
    popt,pcov= model_fit(fit_time,fit_count-baseline,init)
    print(popt)
    A,K,sig,t0=popt
    R = R2(fit_count,model_func(fit_time,*popt)+baseline)
    Radj = R2adj(len(fit_count),3,fit_count,model_func(fit_time,*popt)+baseline)
    rChi = rChi2(len(fit_count),3,fit_count,model_func(fit_time,*popt)+baseline)
#    ax.plot(fit_time,decay_func(fit_time,t0,A_lin,K_lin)+baseline,label=r'decay fit $R^2 =$%.4f%s$\tau_{eff} =$ %.2e ns'%(R,'\n',-1/K_lin))
    print("R2 : %.4f"%R)
    print("R2adj : %.4f"%Radj)
    print("rChi2 : %.4f"%rChi)
    ax.plot(fit_time,model_func(fit_time,*popt)+baseline,c='orange',label=r'simple_fit $1-R^2 =$ %.2e %s $\tau_{eff} =$ %.2e ns %s  $\sigma=$%.2e ns'%((1-R),'\n',-1/K,'\n',sig))

    
    
    print("------------------model2-----------------------")
    init = [A,K,1,1,sig,t0]
    popt,pcov= model2_fit(fit_time,fit_count-baseline,init)
    print(popt)
    A1,K1,A2,K2,sig,t0 = popt
    R = R2(fit_count,model2_func(fit_time,*popt)+baseline)
    Radj = R2adj(len(fit_count),3,fit_count,model2_func(fit_time,*popt)+baseline)
    rChi = rChi2(len(fit_count),3,fit_count,model2_func(fit_time,*popt)+baseline)
    print("R2 : %.4f"%R)
    print("R2adj : %.4f"%Radj)
    print("rChi2 : %.4f"%rChi)
    print("A1/A2 : %.4e"%(A1/A2))
    taueff = (A1*(-1/K1)+A2*(-1/K2))/(A1+A2)
    print("taueff : %.4e ns"%taueff)
    if R>0.8:
        ax.plot(fit_time,model2_func(fit_time,*popt)+baseline,c='red',label=r'double_fit $1-R^2 =$ %.2e %s$A_{1} =$ %.2e %s$A_{2} =$ %.2e %s$\tau_{1} =$ %.2e ns %s$\tau_{2} =$ %.2e ns %s$\sigma =$ %.2e ns %s$\tau_{eff} =$ %.2e ns'%((1-R),'\n',A1,'\n',A2,'\n',-1/K1,'\n',-1/K2,'\n',sig,'\n',taueff))
    
    ax.legend()
    fig.tight_layout()
    if save==True:  
        fig.savefig('%s_double_decay.pdf'%os.path.splitext(path)[0])
        fig.savefig('%s_double_decay.png'%os.path.splitext(path)[0])
        ax.set_yscale('log')
        fig.savefig('%s_double_decay_log.pdf'%os.path.splitext(path)[0])
        fig.savefig('%s_double_decay_log.png'%os.path.splitext(path)[0])
    if autoclose==True:
        plt.close(fig)
    ax.set_yscale('log')
    return A,1/K,A1,A2,1/K1,1/K2,R
def process_all():
    files=getListOfFiles(r'C:/Users/sylvain.finot/Documents/data/')
    mask = [x for x in files if (("TRCL" in x) & (x.endswith(".dat")))]
    for p in mask:
        try:
            process_fromFile(p,save=True,autoclose=True)
        except:
            pass

def main():
    path = input("Enter the path of your file: ")
    path=path.replace('"','')
    path=path.replace("'",'')
#    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 5K/Fil3/TRCL-cw455nm/TRCL.dat'
    process_fromFile(path,save=False,autoclose=False,merge=True)
    
if __name__ == '__main__':
    main()