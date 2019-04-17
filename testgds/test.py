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

substrat = [int(0.9e4),int(0.9e4)]

ALPHA = list(string.ascii_uppercase)
alpha = list(string.ascii_lowercase)

amarks = templates.AlignmentMarks('C') #600x400um
# Create a Cell to hold the objects

#points=[(2,0), (4,0), (4,2), (6,2), (6,4), (4,4),(4,6),(2,6),(2,4),(0,4),(0,2),(2,2)]
points =np.array([(3,0), (3,1), (1,1), (1,3),(0,3),(0,0)])
points = np.concatenate((points,points*-1))
smallCross = core.Boundary(points).scale(10/6)

bigCross = smallCross.copy()
bigCross.scale(10)

#taille lettre
#hauteur = police*0.777777
#largeur = police * 0.555555
cell = core.Cell('EXAMPLE')
#cell.add(amarks, origin=(-200, 0))
#cell.add(amarks, origin=(1200, 0))
for x in range(0,substrat[0]+1,100):
    for y in range(0,substrat[1]+1,100):
        if (x%1000==0) & (y%1000==0):
#            Mark = core.CellReference(amarks,(x,y),magnification=1)
            Mark = utils.translate(bigCross, (x,y))
            letter = shapes.Label('%s'%(ALPHA[x//1000]),40,position=(x-50,y+10))
            number = shapes.Label('%d'%((substrat[1]+1-y)//1000),40,position=(x+25,y-60))
            cell.add(letter)
            cell.add(number)
        else:
#            Mark = core.CellReference(amarks,(x,y),magnification=0.01)
            Mark = utils.translate(smallCross, (x,y))
#            text = shapes.Label('%s%d'%(ALPHA[i//5],j//5),100,position=(x-300,y+50))
        cell.add(Mark)

for x in range(0,substrat[0],1000):
    letter = shapes.Label('%s'%(ALPHA[x//1000]),1000,position=(x+222,substrat[1]-100))
#    number = shapes.Label('%d'%(y//1000),40,position=(x+25,y-60))
    cell.add(letter)
for y in range(0,substrat[1],1000):
    number = shapes.Label('%d'%((substrat[1]-1-y)//1000),1000,position=(-700,y-100))
#    number = shapes.Label('%d'%(y//1000),40,position=(x+25,y-60))
    cell.add(number)
# Create two copies of the Cell
#core.CellReference(cell, origin=(0, 40), magnification=1.5)







#top = core.Cell('TOP')
#cell_array = core.CellArray(cell, 5, 5, (1200, 1200))
#top.add(cell_array)

# Add the copied cell to a Layout and save
layout = core.Layout('LIBRARY')
layout.add(cell)
layout.save('output.gds')

layout.show()