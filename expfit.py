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
from scipy.stats import pearsonr
def decay_func(t,t0,A,K):
    return A*np.exp(K*(t-t0))

def double_decay(t,A1,A2,K1,K2,c):
    return A1*np.exp(K1*t)+A2*np.exp(K2*t)+c

def model_func(t,A,K,sig,t0):
    return 1/2*(1+sse.erf((t-t0)/(np.sqrt(2)*sig)))*A*np.exp(K*(t-t0))

def model2_func(t,A1,K1,A2,K2,sig,t0):
    return 1/2*(1+sse.erf((t-t0)/(np.sqrt(2)*sig)))*(A1*np.exp(K1*(t-t0))+A2*np.exp(K2*(t-t0)))

def stretched_exp(t,A,K,B):
    return A*np.exp(-np.power((K*(t)),B))#*1/2*(1+sse.erf((t)/(np.sqrt(2)*sig)))

def s_fit(t,y,p0=None):
    popt,pcov = curve_fit(stretched_exp,t,y,p0)#,p0,bounds=(lowerbound,upperbound))
    return popt, pcov
def decay_fit(t,y):
    p0 = [min(t),max(y),1]
    popt, pcov = curve_fit(decay_func, t, y,p0)
    return popt,pcov


def fit(t,y,C=0):
    y = np.abs(y - C)
    y = np.log(y)
    t=t-min(t)
    def func(t, A,K,t0):
        return A + K*(t)
    popt, pcov = curve_fit(func, t, y)
    return popt,pcov

def model_fit(t,y,p0):
    A = p0[0]
    K = p0[1]
    sig = p0[2]
    t0 = p0[3]
    lowerbound = (A*0.1,K*1.5,0,0)
    upperbound = (A*10,K*0.6,1,np.inf)

#    print(p0)
#    print(lowerbound)
#    print(upperbound)
    popt,pcov = curve_fit(model_func,t,y,p0,bounds=(lowerbound,upperbound))
    return popt, pcov

def model2_fit(t,y,p0):
    A1 = p0[0]
    K1 = p0[1]
    A2 = p0[2]
    K2 = p0[3]
    sig = p0[4]
    t0 = p0[5]
    lowerbound = (A1*0.5,K1*1.5,-np.inf,-np.inf,0,0)
    upperbound = (A1*2,K1*0.3,np.inf,+np.inf,1,np.inf)

#    print(p0)
#    print(lowerbound)
#    print(upperbound)
    popt,pcov = curve_fit(model2_func,t,y,p0,bounds=(lowerbound,upperbound))
    return popt, pcov

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

def R2(x,y,y_model):
    num = np.sum((y-y_model)**2)
    denom = np.sum((y-y.mean())**2)
#    return 1-num/denom
    return pearsonr(y,y_model)[0]
    
def correl(i,j,cov):
    return cov[i,j]/(cov[i,i]*cov[j,j])**0.5    
def process(path,save=True,autoclose=False,merge=True):
#    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 005K/Fil3/TRCL-cw455nm/TRCL.dat'
#    path = r"C:\Users\sylvain.finot\Documents\data\Probleme TRCL\TRCL_baseline.dat"
#    path = r"C:\Users\sylvain.finot\Documents\data\Probleme TRCL\TRCL_sync.dat"
    path=r"C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2594 - Rampe\300K\TRCL.dat"
    name = path[path.find('2019'):]
    fig, ax = plt.subplots()
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    ax.set_xlabel("t (ns)")
    ax.set_ylabel("Intensity (arb. unit)")
    ax.set_title("Carrier lifetime : %s"%name)
    counts = np.loadtxt(path) #Full histogram
    binNumber = int(counts[0]) #Nombre de bin
    binsize = counts[3] #Taille d'un bin
    counts = counts[-binNumber:] #histogram
    t = np.arange(binNumber)*binsize #echelle de temps en ns
    countmax = counts.argsort()[-1]
    tmin = max(0,countmax-int(1/binsize))
    tmax = min(countmax+int(5/binsize),binNumber)
    reduced_time = t[tmin:tmax]
    reduced_time = reduced_time - min(reduced_time)
    if merge=='auto':
        merge = counts.max() < 5E3
    if merge==True:
        reduced_counts = mergeData(counts,binNumber,binsize,name)
    else:
        reduced_counts = counts[tmin:tmax]
    rightBaseline = np.median(reduced_counts[-int(1/binsize):])
    leftBaseline  = np.median(reduced_counts[:int(1/binsize)])
    baselineError = abs(rightBaseline-leftBaseline)/rightBaseline > 0.50
    ax.plot(reduced_time,reduced_counts,'.',c='k',label='data')
    baseline = rightBaseline
    
    #calcul de la limite de droite pour fit
    leftl = reduced_counts.argmax()+2
    rightl = 0
    c=np.exp(-3)
    while(rightl<leftl or np.isnan(rightl)):        
        threshold = (max(reduced_counts)-baseline)*c #(max(reduced_counts)-baseline)*np.exp(-3)
        mask = np.convolve(np.sign(reduced_counts-threshold),[-1,1],'same') #detect le chgt de sign de reduced-threshold
        mask[0] = 0
        rightl = np.argmax(mask)
        print(rightl)
        c += 0.01
        
    t0 = reduced_time[leftl]#temps correspondant au max
    
    
    #Fit exponential decay
    print("simple decay")
    fit_time = reduced_time[leftl:rightl]
    fit_count = reduced_counts[leftl:rightl]
    popt,pcov= fit(fit_time,fit_count,baseline)
    A_lin=np.exp(popt[0])
    K_lin = popt[1]
    print(popt)
    p_sigma = np.sqrt(np.diag(pcov))
#    print(pcov)
    R = R2(fit_time-t0,fit_count,decay_func(fit_time,t0,A_lin,K_lin)+baseline)
#    ax.plot(fit_time,decay_func(fit_time,t0,A_lin,K_lin)+baseline,label=r'decay fit $R^2 =$%.4f%s$\tau_{eff} =$ %.2e ns'%(R,'\n',-1/K_lin))
    print(R)
     
    #Stretched
    #    popt,pcov = s_fit(fit_time-t0,fit_count-baseline,[A,abs(K),1])
    #    print(popt)
    #    ax.plot(reduced_time,stretched_exp(reduced_time-t0,*popt)+baseline,c='orange',label=r'stretched')
    #    ax.plot(fit_time,stretched_exp(fit_time,A,K,sig,t0,B)+baseline,c='pink',label=r'stretched')
    
    
        #calcul de la limite de gauche
    c=np.exp(-3)
    while(np.isnan(leftl)):        
        threshold = (max(reduced_counts)-leftBaseline)*c #(max(reduced_counts)-baseline)*np.exp(-3)
        mask = np.convolve(np.sign(reduced_counts-threshold),[1,-1],'same') #detect le chgt de sign de reduced-threshold
        mask[0] = 0
        leftl = np.argmax(mask)
        print(leftl)
        c += 0.01
    
    print("------------------model1-----------------------")
    #fit simple exp convoluée avec heaviside
#    if not baselineError:
    fit_time = reduced_time[leftl-20:]
    fit_count = reduced_counts[leftl-20:]
#    else:
#        fit_time = reduced_time[leftl:]
#        fit_count = reduced_counts[leftl:]
    init = [A_lin,K_lin,0.02,t0]
    popt,pcov= model_fit(fit_time,fit_count-baseline,init)
    print(popt)
    print(pcov)
    p_sigma = np.sqrt(np.diag(pcov))
#    print(p_sigma)
    A=popt[0]
    K = popt[1]
    sig = popt[2]
    t0 = popt[3]
    R = R2(fit_time,fit_count,model_func(fit_time,*popt)+baseline)
    print(R)
    ax.plot(fit_time,model_func(fit_time,*popt)+baseline,c='orange',label=r'simple_fit $R^2 =$ %.4f %s $\tau_{eff} =$ %.2e ns %s  $\sigma=$%.2e ns'%(R,'\n',-1/K,'\n',sig))
    
    
    
    
    print("------------------model2-------------------")
    init = [A,K,1,1,sig,t0]
    popt,pcov= model2_fit(fit_time,fit_count-baseline,init)
    print(popt)
    print(pcov)
    A1 = popt[0]
    K1 = popt[1]
    A2 = popt[2]
    K2 = popt[3]
    sig = popt[4]
    t0 = popt[5]
    R = R2(fit_time,fit_count,model2_func(fit_time,*popt)+baseline)
    print(R)
    if R>0.8:
        taueff = (A1*(-1/K1)+A2*(-1/K2))/(A1+A2)
        print(taueff)
        ax.plot(fit_time,model2_func(fit_time,*popt)+baseline,c='red',label=r'double_fit $R^2 =$ %.4f %s$A_{1} =$ %.2e %s$A_{2} =$ %.2e %s$\tau_{1} =$ %.2e ns %s$\tau_{2} =$ %.2e ns %s$\sigma =$ %.2e ns %s$\tau_{eff} =$ %.2e ns'%(R,'\n',A1,'\n',A2,'\n',-1/K1,'\n',-1/K2,'\n',sig,'\n',taueff))
    
    
    
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
    return A,1/K,A1,A2,1/K1,1/K2
def process_all():
    files=getListOfFiles(r'C:/Users/sylvain.finot/Documents/data/')
    mask = [x for x in files if (("TRCL" in x) & (x.endswith(".dat")))]
    for p in mask:
        try:
            process(p,autoclose=False)
        except:
            pass

def main():
    path = input("Enter the path of your file: ")
    path=path.replace('"','')
    path=path.replace("'",'')
#    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 5K/Fil3/TRCL-cw455nm/TRCL.dat'
    process(path,save=False,autoclose=False,merge=True)
    
if __name__ == '__main__':
    main()