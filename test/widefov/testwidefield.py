#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Tests the new iscrn module.  
Note - this may not be kept up to date with the module."""
import numpy
import science.iscrn
import science.xinterp_dm
import science.wideField
import base.readConfig
import util.Ctrl
class Parent:
    dataValid=1
    outputData=numpy.zeros((120),numpy.float32)
recon=Parent()
ctrl=util.Ctrl.Ctrl(globals=globals())
iscrn=science.iscrn.iscrn(None,ctrl.config,idstr="L0-2")
dm=science.xinterp_dm.dm({"recon":recon},ctrl.config,idstr="dma")
wf=science.wideField.WideField({"L0-2":iscrn,"dma":dm},ctrl.config,idstr="a")

execOrder=[iscrn,dm,wf]
ctrl.mainloop(execOrder)

