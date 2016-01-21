# -*- coding: utf-8 -*-
"""
Using DASP modules, build a simple closed-loop model and alter after
forming the interaction matrix to introduce a spider to examine the
effects.

Created on Mon Jan 18 13:06:00 2016

@author: nab
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

## variables
## \/
spiderWidth = 4
nLoops = 1000
## /\
## (end)

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
#wfsOverview_Cspider = config.getVal("wfsOverview")
#wfsObject = config.getVal("wfsObject")
config.setVal("wfsOverview", config.getVal("wfsObject") ) # alter the wfsOverview value

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


nPokeIntegrations = dm.dmflag.sum()+1
print("\t>>>\tPOKING:: Expected for {0:d} integrations".format(
      nPokeIntegrations)
   ) 
shsPokingImg = 0 
for i in range(nPokeIntegrations): # iterate over all actuators to make the PMX
    for obj in ctrl.compList:
        obj.doNextIter()
        shsPokingImg += wfscent.drawCents(0)*nPokeIntegrations**-1.0 

## alter the centroid object's pupil to have a spider
##
pupil = numpy.array( wfscent.wfsobj.pupil.copy() )
pupilShape = pupil.shape[0]
pupil[pupilShape/2-spiderWidth/2:pupilShape/2+spiderWidth/2,:]=0
pupil[:,pupilShape/2-spiderWidth/2:pupilShape/2+spiderWidth/2]=0

## update the wfscent object to use this new spider
##
   # ???
   # ?? WHICH ONE OF THE THREE BELOW IS THE CORRECT ONE TO CHANGE ??
   # ?? -> which is necessary and which are optional ?? 
   # ???
wfscent.wfscentObj.pup   = pupil
wfscent.wfscentObj.pupfn = pupil
wfscent.wfscentObj.pupil = pupil
##??wfscent.wfscentObj.initialiseCmod()
wfscent.wfscentObj.finishInit() ## but this does not change the subap fluxes
wfscent.wfscentObj.subarea = wfscent.wfscentObj._calculateSubAreas() ## this should



## ----------------------------------------------------------
## Prepare for closed loop
## 
## 

print("\t>>>\tCLOSED-LOOP:: Entering for {0:d} integrations".format(
      nLoops)
   ) 

for obj in ctrl.compList:
    if obj.control.has_key("zero_dm"):
        obj.control["zero_dm"]=0
    if obj.control.has_key("zero_science"):
        obj.control["zero_science"]=10
    if obj.control.has_key("science_integrate"):
        obj.control["science_integrate"]=1
    if obj.control.has_key("cal_source"):
        obj.control["cal_source"]=0
    if obj.control.has_key("close_dm"):
        obj.control["close_dm"]=1 # set to zero to not close the loop (!)
    if obj.control.has_key("dm_update"):
        obj.control["dm_update"]=1


## ----------------------------------------------------------
## Do closed loop iterations
## 
## 

dmShapes, shsImg = [],0
print("[ ",end="")
for i in range( nLoops ):
    print("\r[ "+
         "-"*(int(i*nLoops**-1.0*70)) +
         " "*(70-int(i*nLoops**-1.0*70)) +
         " ]",end="") # print a progress bar
    sys.stdout.flush()
    scrn.doNextIter()
    atmos.doNextIter()
    dm.doNextIter()
    wfscent.doNextIter()
    recon.doNextIter()
    #
    if i%10==0: # every 10th iteration, keep: atmosphere, DM sfc & centroids
        dmShapes.append(
                [ i, atmos.outputData.copy(), dm.mirrorSurface.phsOut.copy(),
                  wfscent.outputData.copy() ]
            )
    #phaseScreen = scrn.unwrapPhase("L0")#optional: phase in the screen
    #pupilPhase = atmos.outputData#optional: phase in the pupil
    shsImg += wfscent.drawCents(0)*nLoops**-1.0 #optional: SH spot images

print()

## ----------------------------------------------------------
## Plot the first 10 input pupil phase and associated DM shape
## 
## 

print("\t>>>\tPLOTTING:: ".format()) 

import pylab
pylab.ion()
pylab.figure(1)
for i in range(min(4*6/2,len(dmShapes))):
    pylab.subplot(4,6,2*i+1)
    pylab.imshow( dmShapes[i][1], interpolation='nearest')
    pylab.title("Iteration no. {0:d}".format(dmShapes[i][0]))
    #
    pylab.subplot(4,6,2*i+1+1)
    pylab.imshow( dmShapes[i][2], interpolation='nearest')

pylab.figure(2)
pylab.plot(
      [ dmShapes[j][0] for j in range(len(dmShapes)) ],
      [ dmShapes[j][3].var() for j in range(len(dmShapes)) ],
   )
pylab.xlabel("Iteration #")
pylab.ylabel("WFS centroid variance")
pylab.show()
