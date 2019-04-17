# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 15:15:00 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os.path
import time
def  scaleSEMimage(file):
    with open(file,'r',encoding="Latin-1") as myfile:
        data =myfile.read()
    PixelWidth=eval(data[data.find('PixelWidth=')+11:data.find('PixelWidth=')+23])
    PixelHeight=eval(data[data.find('PixelHeight=')+12:data.find('PixelHeight=')+24])
    width=eval(data[data.find('ResolutionX=')+12:data.find('ResolutionX=')+16])
    height=eval(data[data.find('ResolutionY=')+12:data.find('ResolutionY=')+16])
    Acceleration=eval(data[data.find('HV=')+3:data.find('HV=')+8]) 
    #FullImage =  image.crop((0,0,width,height))
    Totallength_x = PixelWidth *width
    Totallength_y = height*PixelHeight
    xscale = np.linspace(-Totallength_x/2,Totallength_x/2,width)/1e-6
    yscale = np.linspace(-Totallength_y/2,Totallength_y/2,height)/1e-6
    im = Image.open(file).convert('I')
    im = im.crop((0,0,width,height))
    return xscale,yscale,Acceleration,im

def merge(path1,path2,save=True,autoclose=False,log=False,threshold=0,aspect="auto",EnergyRange = []):
    path1=r'C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2581 - 005K\Fil 3\HYP1-T2581-005K-Vacc5kV-spot7-zoom4000x-gr600-slit0-2-t5ms-cw460nm/Hyp.dat'
    path2=r'C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2581 - 005K\Fil 3\HYP1-T2581-005K-Vacc5kV-spot7-zoom4000x-gr600-slit0-2-t5ms-cw540nm/Hyp.dat'
    dirname = os.path.dirname(os.path.dirname(path1))
    dirname1 = os.path.dirname(path1)
    hyppath1 = path1
    specpath1 = os.path.join(dirname1,'Hyp_X_axis.asc')
    filepath1 = os.path.join(dirname1,'Hyp_SEM image after carto.tif')
    data = np.loadtxt(hyppath1)
    xlen1 = int(data[0,0])
    ylen1 = int(data[1,0])
    wavelenght1 = np.loadtxt(specpath1)
    wavelenght1 = wavelenght1[:2048]
    xcoord1 = data[0,1:]
    ycoord1 = data[1,1:]
    CLdata1 = data[2:,1:]
    
    dirname2 = os.path.dirname(path2)
    hyppath2 = path2
    specpath2 = os.path.join(dirname2,'Hyp_X_axis.asc')
    filepath2 = os.path.join(dirname2,'Hyp_SEM image after carto.tif')
    data = np.loadtxt(hyppath2)
    xlen2 = int(data[0,0])
    ylen2 = int(data[1,0])
    wavelenght2 = np.loadtxt(specpath2)
    wavelenght2 = wavelenght2[:2048]
    xcoord2 = data[0,1:]
    ycoord2 = data[1,1:]
    CLdata2 = data[2:,1:]
    
    
    WaveStep = wavelenght2[-1]-wavelenght2[-2] 
    ridx = abs(wavelenght1-wavelenght2.min()).argmin()
    patch = np.arange(wavelenght1[ridx],wavelenght2[0],WaveStep)
    wavelenght = np.concatenate((wavelenght1[:ridx+1],patch,wavelenght2))
    patch = np.ones((len(patch),CLdata2.shape[1])) * np.min((CLdata1,CLdata2)) -10
    CLdata = np.concatenate((CLdata1[:ridx+1],patch,CLdata2))
    
    hypSpectrum = np.transpose(np.reshape(np.transpose(CLdata),(ylen1,xlen1 ,len(wavelenght))), (0, 1, 2))
    
    if len(EnergyRange)==2:
        lidx = np.argmin(np.abs(wavelenght-EnergyRange[0]))
        ridx = np.argmin(np.abs(wavelenght-EnergyRange[1]))
        wavelenght = wavelenght[lidx:ridx]
        hypSpectrum = hypSpectrum[:,:,lidx:ridx]
    average_axis = 0 #1 on moyenne le long du fil, 0 transversalement
    linescan = np.sum(hypSpectrum,axis=average_axis)
    linescan -= linescan.min()
#    linescan = np.log10(linescan)
    xscale_CL,yscale_CL,acc,image = scaleSEMimage(filepath1)
    fig,(ax,bx,cx)=plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios': [1,1, 3]})
    fig.subplots_adjust(top=0.98,bottom=0.11,left=0.1,right=0.82,hspace=0.05,wspace=0.05)
    newX = np.linspace(xscale_CL[int(xcoord1.min())],xscale_CL[int(xcoord1.max())],len(xscale_CL))
    newY = np.linspace(yscale_CL[int(ycoord1.min())],yscale_CL[int(ycoord1.max())],len(yscale_CL))
    
    nImage = np.array(image.crop((xcoord1.min(),ycoord1.min(),xcoord1.max(),ycoord1.max())))#-np.min(image)
#    minC = nImage.min()
#    maxC=nImage.max()
    ax.imshow(nImage,cmap='gray',vmin=0,vmax=65535,extent=[np.min(newX),np.max(newX),np.max(newY),np.min(newY)])
    ax.set_ylabel("distance (µm)")
    hypSpectrum = hypSpectrum-threshold
    hypSpectrum = hypSpectrum*(hypSpectrum>=0)
    hypimage=np.sum(hypSpectrum,axis=2)
    hypimage -= hypimage.min()
    if log:
        hypimage=np.log10(hypimage+1)
        linescan = np.log10(linescan+1)
    lumimage = bx.imshow(hypimage,cmap='jet',extent=[np.min(newX),np.max(newX),np.max(newY),np.min(newY)])
    if average_axis==1:
        extent = [1239.842/wavelenght.max(),1239.842/wavelenght.min(),np.max(newY),np.min(newY)]
        im=bx.imshow(linescan,cmap='jet',extent=extent)
        bx.set_xlabel("energy (eV)")
        bx.set_ylabel("distance (µm)")
    else:
        extent = [np.min(newX),np.max(newX),1239.842/wavelenght.max(),1239.842/wavelenght.min()]
        im=cx.imshow(linescan.T,cmap='jet',extent=extent)
        cx.set_ylabel("Energy (eV)")
        cx.set_xlabel("distance (µm)")

    cx.set_aspect('auto')
    bx.set_aspect(aspect)
    ax.set_aspect(aspect)
    ax.get_shared_y_axes().join(ax, bx)
    pos = cx.get_position().bounds
    cbar_ax = fig.add_axes([0.85, pos[1], 0.05, pos[-1]*0.9])  
    fig.colorbar(im, cax=cbar_ax)
    cbar_ax.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
    
    pos = bx.get_position().bounds
    cbar_ax = fig.add_axes([0.85, pos[1], 0.05, pos[-1]*0.9])
    fig.colorbar(lumimage,cax=cbar_ax)
    cbar_ax.ticklabel_format(axis='both',style='sci',scilimits=(0,0))
    if save==True:
        savedir = os.path.join(dirname,'Saved')
        if not os.path.exists(savedir):
            os.mkdir(savedir)
        fig.savefig(os.path.join(dirname,'Saved','SEM+Hyp+Linescan_merged.png'),dpi=300)
    if autoclose==True:
        plt.close(fig)