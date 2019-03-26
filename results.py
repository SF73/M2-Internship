# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:52:08 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools



T = [300,300,300,5,5,5]
TauwithUL = np.array([0.151,0.145,0.149,3.05E-01,0.367,0.298])
TauwithSl = np.array([0.169,0.155,0.147,0.329,0.312,4.07E-01])
Tauwithout = np.array([0.189,0.175,0.195,2.38E-01,0.247,0.252])

fig,ax=plt.subplots()
markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
ax.plot(T,Tauwithout*1e3,next(markers),label='w/o UL')
ax.plot(T,TauwithSl*1e3,next(markers),label='SL')
ax.plot(T,TauwithUL*1e3,next(markers),label='UL')
ax.set_xlabel("T (K)")
ax.set_ylabel(r"$\tau$ (ps)")
ax.set_title(r"Carrier lifetime @5K & 300K")
ax.legend()

#fig,ax=plt.subplots()
#ax.plot(T,Tauwithout,'*',label='w/o UL')
#ax.plot(T,TauwithSl,'o',label='SL')
#ax.plot(T,TauwithUL,'D',label='UL')
#ax.legend()

#--------------Position------------------
test = np.asarray([[2.5,1.08E-01,0.102,0.0766,0.2,4340,1120],[5,1.22E-01,0.115,8.10E-02,0.199,5.80E+03,2320],[7.5,1.44E-01,0.139,0.11,0.238,7.61E+03,2.31E+03],[10,1.45E-01,0.141,0.116,0.248,11300,2730],[12.5,1.25E-01,0.124,0.112,0.281,1.31E+04,9.64E+02],[15,1.44E-01,0.141,0.096,0.394,4.09E+03,731],[12.25,1.31E-01,0.13,0.113,0.25,1.44E+04,2.01E+03]])
tau_eff1 = test[:,1]*1e3
tau_eff2 = test[:,2]*1e3
tau1 = test[:,3]*1e3
tau2 = test[:,4]*1e3
A1 = test[:,5]
A2 = test[:,6]
x = test[:,0]
fig,(ax,bx)=plt.subplots(1,2,sharex=True)
markers = itertools.cycle(["o",'D',"s",'h','H','8','*'])
ax.plot(x,tau1,next(markers),label=r'$\tau_1$')
ax.plot(x,tau2,next(markers),label=r'$\tau_2$')
ax.plot(x,tau_eff1,next(markers),label=r'$\tau_{eff1}$')
ax.plot(x,tau_eff2,next(markers),label=r'$\tau_{eff2}$')
bx.plot(x,A2/A1,next(markers),label='A2/A1')
ax.set_xlabel("distance from the top (µm)")
ax.set_ylabel(r"$\tau$ (ps)")
ax.set_title(r"linescan T2597 (UL) @300K")
bx.set_xlabel("distance from the top (µm)")
bx.set_ylabel(r"A2/A1")
ax.legend()