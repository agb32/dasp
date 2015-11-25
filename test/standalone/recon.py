# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:26:50 2015

@author: yhz, nab
@email: n.a.bharmal@durham.ac.uk
"""

from __future__ import print_function

# aosim modules
import base.aobase
import base.readConfig
import science.iatmos
import science.iscrn
import science.tomoRecon
import science.wfscent
import science.xinterp_dm
import util.Ctrl

import numpy
import sys
try:
   import pylab
except ImportError:
   no_pylab=None



#
# Connections between modules are::
# 
#  ______       ______        ______        ______         ______ 
# | scrn | --> | atmos|  --> | dm   |  --> | wfs  |   --> |recon |
# |______|     |______|      |______|      |______|       |______|
#                              /\                            |
#                               |____________________________|


## ----------------------------------------------------------
## start of instance creation and setup
## 
## 


## start of instance creation and setup

config=base.readConfig.AOXml("recon.xml") # read the configuration file

scrn=science.iscrn.iscrn(None,config,idstr="L0") # create one phase screen
scrn.finalInitialisation()

atmos=science.iatmos.iatmos({"L0":scrn},config,idstr="a") # create one atmosphere, connected to one phase screen
atmos.finalInitialisation()

dm=science.xinterp_dm.dm(atmos,config,idstr="dma") # create one dm
dm.finalInitialisation()
 
wfscent=science.wfscent.wfscent(dm,config,idstr="a") # create one WFS, connected to the DM
wfscent.finalInitialisation()

recon=science.tomoRecon.recon({"a":wfscent,},config,idstr="ngs") # create one reconstructor, connected to the WFS
recon.finalInitialisation()

dm.newParent({"1":atmos,"2":recon,},"dma") # connect the DM to the atmosphere *and* the reconstructor
# then screen->atmos->DM<-recon, & DM->WFS->recon


## ----------------------------------------------------------
## Create a poke and reconstruction matrix
## 
## 

class ctrl(object): pass

ctrl.compList = [scrn,atmos,dm,wfscent,recon]

for obj in ctrl.compList:
 if obj.control.has_key("poke"):
  obj.control["poke"]=1
 if obj.control.has_key("science_integrate"):
  obj.control["science_integrate"]=0
 if obj.control.has_key("zero_dm"):
  obj.control["zero_dm"]=1
 if obj.control.has_key("cal_source"):
  obj.control["cal_source"]=1



for i in range(50): # iterate 50 times to make the PMX
    for obj in ctrl.compList:
        obj.doNextIter()


## ----------------------------------------------------------
## Prepare for closed loop
## 
## 


for obj in ctrl.compList:
    if obj.control.has_key("zero_dm"):
        obj.control["zero_dm"]=1
    if obj.control.has_key("zero_science"):
        obj.control["zero_science"]=10
    if obj.control.has_key("science_integrate"):
        obj.control["science_integrate"]=1
    if obj.control.has_key("cal_source"):
        obj.control["cal_source"]=0
    if obj.control.has_key("close_dm"):
        obj.control["close_dm"]=1 # set to zero to not close the loop (!)

## ----------------------------------------------------------
## Do closed loop iterations
## 
## 

dmShapes = []
print("[ ",end="")
for i in range(1000):
    print("\r[ "+"-"*(int(i/1000.0*70)) + " "*(70-int(i/1000.0*70)) + " ]",end="") # print a progress bar
    sys.stdout.flush()
    scrn.doNextIter()
    #phaseScreen = scrn.unwrapPhase("L0")#optional: phase in the screen
    atmos.doNextIter()
    #pupilPhase = atmos.outputData#optional: phase in the pupil
    dm.doNextIter()
    if i%100==0: # every 100th iteration, save the atmosphere and DM surface
        dmShapes.append(
                [ i, atmos.outputData.copy(), dm.mirrorSurface.phsOut.copy() ]
            )
    
    wfscent.doNextIter()
    #shsImg = wfscent.drawCents(0)#optional: SH spot images
    recon.doNextIter()

print()

## ----------------------------------------------------------
## Plot the first 10 input pupil phase and associated DM shape
## 
## 

import pylab
pylab.ioff()
pylab.figure(1)
for i in range(10):
    pylab.subplot(4,6,2*i+1)
    pylab.imshow( dmShapes[i][1], interpolation='nearest')
    pylab.title("Iteration no. {0:d}".format(dmShapes[i][0]))
    #
    pylab.subplot(4,6,2*i+1+1)
    pylab.imshow( dmShapes[i][2], interpolation='nearest')
    
pylab.show()
