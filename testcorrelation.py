# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:14:14 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import scipy.special as sse
from expfit import model_func, R2,model_fit,fit,correl
from mergeNoise import mergeData

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
merge="auto"
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

#calcul de la limite de gauche
c=np.exp(-3)
while(np.isnan(leftl)):        
    threshold = (max(reduced_counts)-leftBaseline)*c #(max(reduced_counts)-baseline)*np.exp(-3)
    mask = np.convolve(np.sign(reduced_counts-threshold),[1,-1],'same') #detect le chgt de sign de reduced-threshold
    mask[0] = 0
    leftl = np.argmax(mask)
    print(leftl)
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
ax.plot(fit_time,model_func(fit_time,*popt)+baseline,c='orange',label=r'simple_fit $R^2 =$ %.4f %s $\tau_{eff} =$ %.2e ns %s  $\sigma=$%.2e ns'%(R,'\n',-1/K,'\n',sig))


As=list()
Ks=list()
Rs=list()
fig,bx=plt.subplots()
for a in np.linspace(A*0.5,A*2,200):
    for k in np.linspace(K*0.5,K*2,200):
        popt[0]=a
        popt[1]=k
        As.append(a)
        Ks.append(k)
        R = R2(fit_time,fit_count,model_func(fit_time,*popt)+baseline)
        if R > 0.99:
            ax.plot(fit_time,model_func(fit_time,*popt)+baseline,label=r'$R^2 =$ %.4f %s A =$ %.4f %s $\tau_{eff} =$ %.2e ns'%(R,'\n',a,'\n',-1/k))
        Rs.append(R)
#ax.legend()
im=bx.scatter(As,np.array(Ks),c=Rs,cmap='jet',vmin=0.99,vmax=1)
fig.colorbar(im)
