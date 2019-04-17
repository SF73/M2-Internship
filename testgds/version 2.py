# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:05:56 2019

@author: sylvain.finot
"""

import os.path
from gdsCAD import templates,shapes,core,utils
import string
import numpy as np
# Create some things to draw:

substrat = [int(8e3),int(8e3)]

ALPHA = list(string.ascii_uppercase)
alpha = list(string.ascii_lowercase)

points =np.array([(3,0), (3,1), (1,1), (1,3),(0,3),(0,0)])
points = np.concatenate((points,points*-1))
smallCross = core.Boundary(points).scale(2)
bigCross = smallCross.copy()
bigCross.scale(10)



def getChar(char,height,x,y):
    police = height/0.77777
    return shapes.Label(char,police,position=(x,y-police*0.22222))

#taille lettre
#hauteur = police*0.777777
#largeur = police * 0.555555
#yOffset = 0.2222*police
#fontAspectRatio = np.array([0.5555555,0.7777777,0.222222]) #width, height, offset
smallCharsize = (2/3)*60
bigCharsize = 800
cell = core.Cell('EXAMPLE')
xOffset = 1000
yOffset = -1000
for x in range(0,substrat[0]+1,100):
    for y in range(0,-(substrat[1]+1),-100):
        if (x%1000==0) & (y%1000==0):
#            Mark = core.CellReference(amarks,(x,y),magnification=1)
            Mark = utils.translate(bigCross, (x+xOffset,y+yOffset))
#            letter = shapes.Label('%s'%(ALPHA[x//1000]),smallPolice[1],position=(x,y+smallPolice[2]))
            letter = getChar(ALPHA[x//1000],smallCharsize,x-60+(smallCharsize/0.77777*0.55555)/2+xOffset,y+60-smallCharsize+yOffset)
            number = getChar('%d'%((-y)//1000),smallCharsize,x+(smallCharsize/0.77777*0.55555)/2+xOffset,y-(3/2)*smallCharsize+yOffset)
#            number = shapes.Label('%d'%((-y)//1000),smallPolice[1],position=(x,y))
            cell.add(letter)
            cell.add(number)
        else:
#            Mark = core.CellReference(amarks,(x,y),magnification=0.01)
            Mark = utils.translate(smallCross, (x+xOffset,y+yOffset))
#            text = shapes.Label('%s%d'%(ALPHA[i//5],j//5),100,position=(x-300,y+50))
        cell.add(Mark)

for x in range(0,substrat[0],1000):
#    letter = getChar(ALPHA[x//1000],bigCharsize,x+(1000-bigCharsize/0.77777*0.55555)/2,100)
    letter = getChar(ALPHA[x//1000],bigCharsize,x-bigCharsize/0.77777*0.55555/2+xOffset,100+yOffset)
    #shapes.Label('%s'%(ALPHA[x//1000]),bigCharsize,position=(x+(bigCharsize/0.77777*0.55555)/2,0))
    cell.add(letter)
for y in range(0,-(substrat[1]),-1000):
#    number = shapes.Label('%d'%((-y)//1000),1000,position=(-700,y))
    number = getChar('%d'%((-y)//1000),bigCharsize,-(3/2)*bigCharsize/0.77777*0.55555-100+xOffset,y-bigCharsize/2+yOffset)
    cell.add(number)
# Create two copies of the Cell
#core.CellReference(cell, origin=(0, 40), magnification=1.5)

cell.add(utils.translate(bigCross.scale(10/6), (200,-200)))




#top = core.Cell('TOP')
#cell_array = core.CellArray(cell, 5, 5, (1200, 1200))
#top.add(cell_array)

# Add the copied cell to a Layout and save
layout = core.Layout('LIBRARY')
layout.add(cell)
layout.save('output.gds')

layout.show()