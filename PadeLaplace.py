# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:06:34 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
from models import gaussian_heaviside,gaussian_decay,DLT
from scipy.optimize import curve_fit
from stats import R2, R2adj
from mergeNoise import mergeData
from expfit import process_fromArray
def LaplaceTrapezeIntegral(t,data,i,p0): #p0
    #good p0 = 1/t0.5 inverse time it takes for the data to decay to 1/2 of its initial value
    dt = t[1]-t[0]
    extTerms = 0.5*(np.power(-t[0],i)*np.exp(-p0*t[0])*data[0]+np.power(-t[-1],i)*np.exp(-p0*t[-1])*data[-1])
    sumTerms = np.sum(np.power(-t[1:-1],i)*np.exp(-p0*t[1:-1])*data[1:-1])
    return dt*(extTerms+sumTerms) #d(i)L/dp(i)

def Calculate_di(t,data,p0,n):
    dis = np.zeros(2*n)
    for i in np.arange(0,2*n):
        dis[i] = 1/np.math.factorial(i) * LaplaceTrapezeIntegral(t,data,i,p0)
    return dis

def MakeMatrix(dis):
    n=len(dis)//2
    test=np.flip(dis[0:n])
    for m in range(1,n):
        test = np.row_stack((test,np.flip(dis[m:n+m])))
    return test

def Calculate_bi(dis):
    matrix = MakeMatrix(dis)
    n=len(dis)//2
    vect = -1*dis[n:2*n]
    if matrix.shape[0] >1:
        bis = np.linalg.inv(matrix).dot(vect)
    else:
        bis = vect/matrix[0]
    return bis
    
def Calculate_ai(bis,dis):
    n = len(dis)//2
    ais = np.zeros(n)
    for i in range(0,n):
        d = np.flip(dis[0:i+1])
        b = np.concatenate(([1],bis[0:i]))
        ais[i] = np.sum(d*b)
    return ais

def decay_gen(x,*args):
    n = len(args)
    baseline =  args[-1]
    t0 = args[-2]
    A = args[:(n-2)//2]
    gam = args[(n-2)//2:-2]
    return decay(x,A,gam,t0,baseline)
def convolvedecay(x,*args):
    n = len(args)
    sig =  args[-1]
    t0 = args[-2]
    baseline = args[-3]
    A = args[:(n-3)//2]
    gam = args[(n-3)//2:-3]
    return gaussian_heaviside(x,sig,t0)*decay(x,A,gam,t0,baseline)
def find_exp_decay(t,data,n,p0=None,t0=None,removebg=True):
    background = 0
    if removebg:
        background = np.median(data[-1000:])
        data = data - background
    if p0 is None:
        p0 = 1/(t[np.argmin(abs(data-max(data)/2))]-t[0])
        print('p0 : %.4e'%p0)
    if t0 is None:
        t0=t[0]
    dis = Calculate_di(t-t0,data,p0,n)
    bis = Calculate_bi(dis)
    ais = Calculate_ai(bis,dis)
    A,gamma,_= signal.residue(np.flip(ais),np.concatenate((np.flip(bis),[1])))
    gamma = abs(gamma) - p0
#    error = np.sum((decay(t,A,gamma,t0,background)-(data+background))**2/(data+background))
    error = R2(data+background,decay(t,A,gamma,t0,background))
    return A,gamma,t0,background,error

def Construct_pade(x,p0,ais,bis):
    n = len(ais)
    pis = np.zeros(n)
    for i in range(0,n+1):
        pis[i] = np.power(x-p0,i)
    num = np.sum(ais*pis)
    bs = np.concatenate(([1],bis))
    denom = np.sum(bs*pis)
    return num/denom
def decay(t,A,Gam,t0,bg):
    y=0
    for i in range(len(A)):
        y+= A[i]*np.exp(-Gam[i]*(t-t0))
    return y+bg
def prepareData(t,data):
    idmax = np.argmax(data)
    t0 = t[idmax]
#    cmax = data[idmax]
    return t[idmax:],data[idmax:],t0
