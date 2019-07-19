# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:57:46 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path
import time
import scipy.constants as cst
eV_To_nm = cst.c*cst.h/cst.e*1E9
from detect_dead_pixels import correct_dead_pixel
from SEMimage import scaleSEMimage
from FileHelper import getListOfFiles
from matplotlib.widgets import MultiCursor, Cursor
plt.style.use('Rapport')
plt.rc('axes', labelsize=16)

def process_all(path):
    files=getListOfFiles(path)
    mask = [x for x in files if (("Hyp" in x) & (x.endswith(".dat")))]
    print(mask)
    for i in range(len(mask)):
        try:
            make_linescan(mask[i],save=True,deadPixeltol=200,autoclose=True,Linescan=True)
            print(100*i/len(mask))
        except Exception as e:
            print(e)
            pass
        
        
class MyMultiCursor(MultiCursor):
    def __init__(self, canvas, axes, useblit=True, horizOn=[], vertOn=[], **lineprops):
        super(MyMultiCursor, self).__init__(canvas, axes, useblit=useblit, horizOn=False, vertOn=False, **lineprops)

        self.horizAxes = horizOn
        self.vertAxes = vertOn

        if len(horizOn) > 0:
            self.horizOn = True
        if len(vertOn) > 0:
            self.vertOn = True

        xmin, xmax = axes[-1].get_xlim()
        ymin, ymax = axes[0].get_ylim()
        xmid = 0.5 * (xmin + xmax)
        ymid = 0.5 * (ymin + ymax)

        self.vlines = [ax.axvline(xmid, visible=True, **lineprops) for ax in self.vertAxes]
        self.hlines = [ax.axhline(ymid, visible=True, **lineprops) for ax in self.horizAxes]
        
def make_linescan(path,save=False,autoclose=False,log=False,threshold=0,deadPixeltol = 2,aspect="auto",EnergyRange = [],normalise=False,Linescan=False):
#    path=r'D:/M2 Internship/data/2019-03-08 - T2601 - 300K/Fil2/HYP1-T2601-300K-Vacc5kV-spot7-zoom6000x-gr600-slit0-2-t5ms-Fil1-cw440nm/Hyp.dat'
#    path = r"C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2594 - Rampe\300K\HYP1-T2594-310K-Vacc5kV-spot7-zoom6000x-gr600-slit0-2-t5ms-cw440nm\Hyp.dat"
#    path=r'C:/Users/sylvain.finot/Documents/data/2019-04-29 - T2594 Al - 300K/Fil 2/HYP-T2594Al-300K-Vacc5kV-spot5-zoom10000x-gr600-slit0-2-t005ms-cw450nm/Hyp.dat'
#    path=r'D:/M2 Internship/data/2019-05-16 - T2455 Al - 300K/FULLSCREEN HYP2-T2455-300K-Vacc5kV-spot5-zoom8000x-gr600-slit0-2-t005ms-cw440nm/Hyp.dat'
    if autoclose:
        plt.ioff()
    else:
        plt.ion()
    dirname,filename = os.path.split(path)
    filename, ext = os.path.splitext(filename)
    try:
        infos = os.path.basename(dirname).split("-")
        T = '' if not [x for x in infos if ('K' in x)] else [x for x in infos if ('K' in x)][0]
        Sample = infos[1]
        wire = '' if not [x for x in infos if ('fil' in x.lower())] else [x for x in infos if ('fil' in x.lower())][0]
        if not wire:
            wire = os.path.basename(os.path.dirname(dirname))
        wire = wire.lower().replace('fil','Wire')
    except:
        T=''
        Sample = ''
        wire = ''
    hyppath = path
    specpath = os.path.join(dirname,filename+'_X_axis.asc')
    filepath = os.path.join(dirname,filename+'_SEM image after carto.tif')
    data = np.loadtxt(hyppath)
    xlen = int(data[0,0])   #nbr de pts selon x
    ylen = int(data[1,0])   #nbr de pts selon y
    wavelenght = np.loadtxt(specpath)
    wavelenght = wavelenght[:2048] #bins du spectro
    xcoord = data[0,1:]
    ycoord = data[1,1:]
    CLdata = data[2:,1:] #tableau de xlen * ylen points (espace) et 2048 longueur d'onde CLdata[:,n] n = numero du spectr
    hypSpectrum = np.transpose(np.reshape(np.transpose(CLdata),(ylen,xlen ,len(wavelenght))), (0, 1, 2))
    
    
    #correct dead / wrong pixels
    hypSpectrum, hotpixels = correct_dead_pixel(hypSpectrum,tol=deadPixeltol)
    
    #reduce the spectrum if wanted
    if len(EnergyRange)==2:
        lidx = np.argmin(np.abs(wavelenght-EnergyRange[0]))
        ridx = np.argmin(np.abs(wavelenght-EnergyRange[1]))
        wavelenght = wavelenght[lidx:ridx]
        hypSpectrum = hypSpectrum[:,:,lidx:ridx]
    
    linescan = np.sum(hypSpectrum,axis=0)
    if normalise:
        linescan = linescan/np.max(linescan,axis=1,keepdims=True)
#    linescan -= linescan.min()
    xscale_CL,yscale_CL,acceleration,image = scaleSEMimage(filepath)
    if Linescan:
        fig,(ax,bx,cx)=plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios': [1,1, 3]})
    else:
        fig,(ax,bx,cx)=plt.subplots(3,1,sharex=True,sharey=True)
    fig.patch.set_alpha(0) #Transparency style
    fig.subplots_adjust(top=0.9,bottom=0.12,left=0.15,right=0.82,hspace=0.1,wspace=0.05)
    newX = np.linspace(xscale_CL[int(xcoord.min())],xscale_CL[int(xcoord.max())],len(xscale_CL))
    newY = np.linspace(yscale_CL[int(ycoord.min())],yscale_CL[int(ycoord.max())],len(yscale_CL))
    convX = np.linspace(newX.min(),newX.max(),hypSpectrum.shape[1])
    convY = np.linspace(newY.min(),newY.max(),hypSpectrum.shape[0])
    
    nImage = np.array(image.crop((xcoord.min(),ycoord.min(),xcoord.max(),ycoord.max())))
    SEM = ax.imshow(nImage,cmap='gray',vmin=0,vmax=65535,extent=[np.min(newX),np.max(newX),np.max(newY),np.min(newY)])
    hypSpectrum = hypSpectrum-threshold
    hypSpectrum = hypSpectrum*(hypSpectrum>=0)
    hypimage=np.sum(hypSpectrum,axis=2)
    hypimage -= hypimage.min()
    if log:
        hypimage=np.log10(hypimage+1)
        linescan = np.log10(linescan+1)
    lumimage = bx.imshow(hypimage,cmap='jet',extent=[np.min(newX),np.max(newX),np.max(newY),np.min(newY)])
    if Linescan:
#        im=cx.imshow(linescan.T,cmap='jet',extent=[np.min(newX),np.max(newX),eV_To_nm/wavelenght.max(),eV_To_nm/wavelenght.min()])
        im=cx.pcolormesh(np.linspace(np.min(newX),np.max(newX),hypSpectrum.shape[1]),eV_To_nm/wavelenght,linescan.T,cmap='jet')
        cx.set_ylabel("Energy (eV)")
        cx.set_xlabel("distance (µm)")
        cx.set_aspect('auto')
    else:
        im = cx.imshow(wavelenght[np.argmax(hypSpectrum,axis=2)],cmap='viridis',extent=[np.min(newX),np.max(newX),np.max(newY),np.min(newY)])    
        cx.set_aspect(aspect)
    bx.set_aspect(aspect)
    ax.set_aspect(aspect)
    ax.get_shared_y_axes().join(ax, bx)
    ax.set_ylabel("distance (µm)")
#    fig.text(ax.get_position().bounds[0]-0.11, ax.get_position().bounds[1],'distance (µm)',fontsize=16, va='center', rotation='vertical')
    
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
        savename = 'SEM+Hyp+Linescan'
        if ((len(T)>0) & (len(Sample)>0) & (len(wire)>0)):
            savename='%s_%s_%s_%s'%(savename,Sample,T,wire)
        fig.savefig(os.path.join(dirname,'Saved',savename+".png"),dpi=300)
            
    if autoclose==True:
        plt.close(fig)
        return None
    else:
        spec_fig, spec_ax = plt.subplots()
        spec_data, = spec_ax.plot(eV_To_nm/wavelenght,np.zeros(wavelenght.shape[0]))
        if Linescan:    
            multi = MyMultiCursor(fig.canvas,(ax, bx,cx), color='r',lw=1,horizOn=[ax,bx], vertOn=[ax,bx,cx],useblit=True)#, lw=1, horizOn=[ax,bx], vertOn=[ax,bx,cx])
#            multi = Cursor(bx, color='r',lw=1,useblit=True)#, lw=1, horizOn=[ax,bx], vertOn=[ax,bx,cx])
        def onclick(event):
            x = event.xdata
            y = event.ydata
            indx = np.argmin(abs(x-convX))
            indy=np.argmin(abs(y-convY))
#            spec_data.set_xdata(
            spec_data.set_ydata(hypSpectrum[indy,indx])
            spec_ax.relim()
            spec_ax.autoscale_view(scalex=False)
            spec_fig.canvas.draw_idle()
            plt.draw()
        fig.canvas.mpl_connect('button_press_event', onclick)
        return spec_fig, spec_ax, multi

def main():
    path = input("Enter the path of your file: ")
    path=path.replace('"','')
    path=path.replace("'",'')
    global spec_fig,spec_ax,cursor
#    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 5K/Fil3/TRCL-cw455nm/TRCL.dat'
    cursor = make_linescan(path,save=False,deadPixeltol=200,normalise=False,Linescan=True)
    
if __name__ == '__main__':
    main()