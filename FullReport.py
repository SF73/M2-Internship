# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:57:46 2019

@author: sylvain.finot
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os.path
import time
from detect_dead_pixels import correct_dead_pixel
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
def process_network():
    dirName = "\\\\srv-echange\echange\Sylvain"
    while(True):
        print('Processing')
        dirs = os.listdir(dirName)
        selection = [x for x in dirs if (time.strftime("%Y-%m-%d") in x)]
        for entry in selection:
            fullPath = os.path.join(dirName, entry)
            toProcess = os.path.exists(os.path.join(fullPath,'process.dat')) & (not os.path.exists(os.path.join(fullPath,'done.dat')))
            if toProcess==True:
                print(fullPath)
                process_all(fullPath)
                os.path.exists(os.path.join(fullPath,'process.dat'))
                np.savetxt(os.path.join(fullPath,'done.dat'),[])
        time.sleep(30)
def process_all(path):
    files=getListOfFiles(path)
    mask = [x for x in files if (("Hyp" in x) & (x.endswith(".dat"))& ("T2594" in x))]
    print(mask)
    for i in range(len(mask)):
        try:
            make_linescan(mask[i],save=True,autoclose=True)
            print(100*i/len(mask))
        except:
            pass
        
def make_linescan(path,save=False,autoclose=False,log=False,threshold=0,deadPixeltol = 2,aspect="auto",EnergyRange = []):
#    path=r'C:/Users/sylvain.finot/Documents/data/2019-03-08 - T2601 - 300K/Fil2/HYP1-T2601-300K-Vacc5kV-spot7-zoom6000x-gr600-slit0-2-t5ms-Fil1-cw380nm/Hyp.dat'
#    path = r"C:\Users\sylvain.finot\Documents\data\2019-03-22 - T2594 - Rampe\300K\HYP1-T2594-310K-Vacc5kV-spot7-zoom6000x-gr600-slit0-2-t5ms-cw440nm\Hyp.dat"
    if autoclose:
        plt.ioff()
    else:
        plt.ion()
    dirname = os.path.dirname(path)
    hyppath = path
    specpath = os.path.join(dirname,'Hyp_X_axis.asc')
    filepath = os.path.join(dirname,'Hyp_SEM image after carto.tif')
    data = np.loadtxt(hyppath)
    xlen = int(data[0,0])
    ylen = int(data[1,0])
    wavelenght = np.loadtxt(specpath)
    wavelenght = wavelenght[:2048]
    xcoord = data[0,1:]
    ycoord = data[1,1:]
    CLdata = data[2:,1:] #tableau de xlen * ylen points (espace) et 2048 longueur d'onde CLdata[:,n] n = numero du spectr
    hypSpectrum = np.transpose(np.reshape(np.transpose(CLdata),(ylen,xlen ,len(wavelenght))), (0, 1, 2))
    
    
    #correct dead / wrong pixels
    hypSpectrum = correct_dead_pixel(hypSpectrum,tol=deadPixeltol)
    
    #reduce the spectrum if wanted
    if len(EnergyRange)==2:
        lidx = np.argmin(np.abs(wavelenght-EnergyRange[0]))
        ridx = np.argmin(np.abs(wavelenght-EnergyRange[1]))
        wavelenght = wavelenght[lidx:ridx]
        hypSpectrum = hypSpectrum[:,:,lidx:ridx]
    
    average_axis = 0 #1 on moyenne le long du fil, 0 transversalement
    
    linescan = np.sum(hypSpectrum,axis=average_axis)
#    linescan -= linescan.min()
    xscale_CL,yscale_CL,acceleration,image = scaleSEMimage(filepath)
    
    fig,(ax,bx,cx)=plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios': [1,1, 3]})
    fig.patch.set_alpha(0) #Transparency style
    fig.subplots_adjust(top=0.98,bottom=0.11,left=0.1,right=0.82,hspace=0.05,wspace=0.05)
    newX = np.linspace(xscale_CL[int(xcoord.min())],xscale_CL[int(xcoord.max())],len(xscale_CL))
    newY = np.linspace(yscale_CL[int(ycoord.min())],yscale_CL[int(ycoord.max())],len(yscale_CL))
    nImage = np.array(image.crop((xcoord.min(),ycoord.min(),xcoord.max(),ycoord.max())))
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
        #pas utilise pour le moment
        extent = [1239.842/wavelenght.max(),1239.842/wavelenght.min(),np.max(newY),np.min(newY)]
        im=bx.imshow(linescan,cmap='jet',extent=extent)
        cx.set_xlabel("energy (eV)")
        cx.set_ylabel("distance (µm)")
    else:
        #par defaut
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
        fig.savefig(os.path.join(dirname,'Saved','SEM+Hyp+Linescan.png'),dpi=300)
    if autoclose==True:
        plt.close(fig)

#def main():
#    path = input("Enter the path of your file: ")
#    path=path.replace('"','')
#    path=path.replace("'",'')
##    path = r'C:/Users/sylvain.finot/Documents/data/2019-03-11 - T2597 - 5K/Fil3/TRCL-cw455nm/TRCL.dat'
#    make_linescan(path,save=True)
#    
#if __name__ == '__main__':
#    main()