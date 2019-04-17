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
from expfit import model_func, R2,model_fit,fit,correl,model2_fit,model2_func
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
merge=True
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
ax.plot(model_func(fit_time,*popt)+baseline,c='orange',label=r'simple_fit $R^2 =$ %.4f %s $\tau_{eff} =$ %.2e ns %s  $\sigma=$%.2e ns'%(R,'\n',-1/K,'\n',sig))

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
R = R2(fit_count,model2_func(fit_time,*popt)+baseline)
print(R)
if R>0.8:
    taueff = (A1*(-1/K1)+A2*(-1/K2))/(A1+A2)
    print(taueff)
    ax.plot(fit_time,model2_func(fit_time,*popt)+baseline,c='red',label=r'double_fit $R^2 =$ %.4f %s$A_{1} =$ %.2e %s$A_{2} =$ %.2e %s$\tau_{1} =$ %.2e ns %s$\tau_{2} =$ %.2e ns %s$\sigma =$ %.2e ns %s$\tau_{eff} =$ %.2e ns'%(R,'\n',A1,'\n',A2,'\n',-1/K1,'\n',-1/K2,'\n',sig,'\n',taueff))
minR =R
A1s=list()
K1s=list()
A2s=list()
K2s=list()
Rs=list()
fig,bx=plt.subplots()
for a1 in np.linspace(A1*0.3,A1*3,200):
    for a2 in np.linspace(A2*0.3,A2*3,200):
        popt[0]=a1
        popt[2]=a2
        A1s.append(a1)
        A2s.append(a2)
        R = R2(fit_time,fit_count,model2_func(fit_time,*popt)+baseline)
        if R >= minR:
            ax.plot(fit_time,model2_func(fit_time,*popt)+baseline)#,label=r'$R^2 =$ %.4f %s A =$ %.4f %s $\tau_{eff} =$ %.2e ns'%(R,'\n',a,'\n',-1/k))
        Rs.append(R)
#ax.legend()
im=bx.scatter(np.array(A1s),np.array(A2s),c=Rs,cmap='jet',vmin=0.99,vmax=1)
bx.scatter(A1,A2,c='cyan',marker='*')
fig.colorbar(im)
