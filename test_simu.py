# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 09:28:41 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from models import gaussian_heaviside,gaussian_decay,DLT,LT
from PadeLaplace import prepareData,find_exp_decay, decay_gen, convolvedecay
from scipy.optimize import curve_fit
from scipy.stats import lognorm, norm
from stats import R2, R2adj
from mergeNoise import mergeData
from expfit import process_fromArray
import pickle

gaussian_density = lambda x:norm.pdf(x,5,1.7)*(np.sign(x)+1)/2
log_norm = lambda x:lognorm()
gaussian = LT(gaussian_density)

def test_dlt(t,t0,gam,A):
    A=A/np.sum(A)
    return DLT(t,A,gam,t0)
A1s = []
A2s = []
tau1 = []
tau2 = []
Aths = []
Gamths = []
taueffTh = []
taueff = []
for k in range(20):
    try:
        print(k/20)
        binsize = 0.004
        t=np.arange(2500)*binsize
        
        size = np.random.randint(0,5)
        Gamth = np.random.uniform(0,20,4001)
        rtime = 0.026#0.026#np.random.uniform(0.025,0.06)
        
#        #Gaussian
#        mu = 6#np.random.uniform(4,8)
#        sig = 5/3.5#np.random.uniform(0.5,mu/2.5)
#        Ath = norm.pdf(Gamth,mu,sig)*(Gamth>0)
        
##        lognorm
        mu = 6#np.random.uniform(1,5)
        shape = np.random.uniform(0.5,2)
        sig = np.random.uniform(0.5,2)
        Ath = lognorm(s=shape,loc=mu,scale=sig).pdf(Gamth)
        
#                Gaussian
#        mu1,mu2 = 6.5,2.37#np.random.uniform(4,8)
#        sig1,sig2 = np.random.uniform(0,mu1/5),np.random.uniform(0,mu2/5)
#        Ath = (1.6*norm.pdf(Gamth,mu1,sig1)+0.6*norm.pdf(Gamth,mu2,sig2))*(Gamth>0)
        
        
        
        t0=5
        true_taueff = np.sum(Ath/Gamth)/np.sum(Ath)
        simu_data = np.zeros(2500)
        for i in range(size):
            t0=t0+np.random.normal(0,0.003)
        #    model_data = gaussian_heaviside(t,0.027,t0+i*periode)*decay(t,Ath,Gamth,t0+i*periode,0)
            model_data = gaussian_heaviside(t,rtime,t0)*test_dlt(t,(t0),Gamth,Ath)*np.random.uniform(5e3,2e4)+30
        #    model_data = gaussian_heaviside(t,0.027,t0+i*periode)*gaussian_decay(t-(t0+i*periode),Ath[0],sig,mu)
#            simu_data += (model_data).astype(np.int)
            simu_data += (np.round(np.random.poisson(model_data))).astype(np.int)
        simu_data = simu_data/size
        #plt.semilogy(t,model_data+30,label="model")
        #plt.semilogy(t,simu_data,'.',label="Simu wo gaussian")
        ntime,ndata = mergeData(t,simu_data,binsize)
        prepared_data = prepareData(t,simu_data)[:-1]
        ##plt.plot(ntime,ndata)
        fig, (ax,bx) = plt.subplots(1,2)
#        ax.plot(*prepared_data,'.',c='k',label='data')
        bx.bar(Gamth,Ath,width=0.1,label='true distribution')
        for i in range(1,100):
            A,Gam,t0,bg,error = find_exp_decay(*prepared_data,i)
            if not np.any(np.imag(A)>0)and error>0.9:
                print("A",A)
                print("Gamma",Gam)
                bx.bar(Gam,A/np.sum(A),width=0.1,label='%1d poles found'%i)
                print("error %.4e"%error)
#                taueff = np.sum(A/Gam)/np.sum(A)
#                print('taueff %.4e'%taueff)
                ax.semilogy(prepared_data[0],DLT(prepared_data[0],A,Gam,t0)+bg,label='%1d poles found'%i)
                ax.legend()
                bx.legend()
        taueff.append(np.sum(A/Gam)/np.sum(A))
#        popt,pcov = curve_fit(decay_gen,prepared_data[0],prepared_data[1],p0=(*A,*Gam,t0,bg))
#        plt.semilogy(t,decay_gen(t,*popt),label="After curve fit")
#        popt,pcov = curve_fit(convolvedecay,prepared_data[0],prepared_data[1],p0=(*A,*Gam,t0,bg,0))
#        plt.semilogy(prepared_data[0],convolvedecay(prepared_data[0],*popt),label="After")

        _,_,A1,A2,t1,t2,R=process_fromArray(t,simu_data,show=False)
        if R>0.99:
            s = A1+A2
            K1,K2 = -1/t1,-1/t2
    #        ax.bar([K1,K2],[A1/s,A2/s],width=0.1,label='Found decay rates')
    #        ax.bar(1/true_taueff,1,width=0.1,label='True mean')
    #        ax.bar(s/(A1*(-t1)+A2*(-t2)),1,width=0.1,label=r'$1/<\tau>$')
    #        ax.bar((A1*(K1)+A2*(K2))/s,1,width=0.1,label=r'$<1/\tau>$')
    #        ax.legend()
            Aths.append(Ath)
            Gamths.append(Gamth)
            A1s.append(A1)
            A2s.append(A2)
            tau1.append(-t1)
            tau2.append(-t2)
            taueffTh.append(true_taueff)
    except:
        pass
#A1s = np.array(A1s)
#A2s = np.array(A2s)
#tau1 = np.array(tau1)
#tau2 = np.array(tau2)
taueffTh = np.array(taueffTh)
taueff=np.array(taueff)
fig,(ax,bx)=plt.subplots(1,2)
#ax.plot(A1s/A2s,'.')
#bx.plot((A1s*tau1+A2s*tau2)/(A1s+A2s),'.',label="Fitted taueff")
#bx.plot(tau1,'.',label="Fitted tau1")
#bx.plot(tau2,'.',label="Fitted tau2")
ax.plot(taueff/taueffTh,'.',label="Simulated")


#fig,ax=plt.subplots()
#ax.bar(Gamths[10],Aths[10]/np.sum(Aths[10]),width=0.1)
#ax.bar([1/tau1[10],1/tau2[10]],[A1s[10]/(A1s[10]+A2s[10]),A2s[10]/(A1s[10]+A2s[10])],width=0.1)
#with(open('gaussian_simu.pickle','wb')) as f:
#    pickle.dump(fig,f)
#with(open('gaussian_simu.pickle', 'rb'))as f:
#    figx = pickle.load(f)
#figx = pickle.load(open('gaussian_simu.pickle', 'rb'))
#figx.show()