# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 18:34:46 2019

@author: Sylvain
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
import matplotlib.patches as patches
# definition des variables
l=5E-3 #longueur du beam blanker
d=20E-2 #distance entre sortie du BB et l'échantillon
w=100E-6 #entrefer du BB
Ec = 5E3* cst.electron_volt #energie des electrons
Emax=700E-3/w#/10**(20/10) #champ appliqué 

Wsample = 2E-6 #taille de l'echantillon
tm = 400E-12 #temps de montée

total_time = (l+d)*np.sqrt(cst.electron_mass/(2*Ec))
tl = (l)*np.sqrt(cst.electron_mass/(2*Ec)) #temps pour traverser le BB
t = np.linspace(0,total_time,10000)

z=(l+d)-np.sqrt((2*Ec)/cst.electron_mass)*t
fig,ax = plt.subplots()
endx = []
xf = []
for E in np.linspace(-Emax,Emax,100):
    x1 = -cst.electron_volt*E/(2*cst.electron_mass)*t**2
    x1[t>tl] = 0
    x2 = -cst.electron_volt*E/cst.electron_mass * tl*(-tl/2+t)
    x2[t<tl]=0
    x= x1+x2        
    ax.plot(x,z,c='k',alpha=0.1)
    endx.append(x[-1])
    #xf = - cst.electron_volt*E/(2*Ec)*l*(l/2+d)
    xf.append(- cst.electron_volt*E/(2*Ec)*l*(l/2+d))

endx = np.asarray(endx)
xf = np.asarray(xf)
#if(max(xf)>w):
ax.vlines(-w/(2),d,d+l)
ax.vlines(+w/(2),d,d+l)

beamblankerL = patches.Rectangle((-w/(2)-2E-6,d),2E-6,l,linewidth=1,edgecolor='green',facecolor='green') 
beamblankerR = patches.Rectangle((w/(2),d),2E-6,l,linewidth=1,edgecolor='green',facecolor='green')
ax.add_patch(beamblankerL)
ax.add_patch(beamblankerR)
# Create a Rectangle patch
sample = patches.Rectangle((-Wsample/2,0),Wsample,Wsample*5000,linewidth=1,edgecolor='r',facecolor='r')
ax.add_patch(sample)
ax.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
print('amplitude de deviation : %.2e'%(xf[0]-xf[-1]))
print('temps sur echantillon : %.2e'%(Wsample*tm/(xf[0]-xf[-1])))
#fig,ax = plt.subplots()
#ax.plot(endx,'.')