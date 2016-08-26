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
import science.wfscent
import science.tomoRecon
import science.iatmos
import science.science
import base.readConfig
import base.saveOutput
import util.Ctrl
import sys
ctrl=util.Ctrl.Ctrl(globals=globals())
ctrl.initialCommand("wf.control['cal_source']=1",freq=-1,startiter=0)
ctrl.initialCommand("wf.control['cal_source']=0",freq=-1,startiter=1)
ctrl.initialCommand("c.newCorrRef();print 'Done new corr ref'",freq=-1,startiter=1)
if "--user=poke" in sys.argv:
    ctrl.doInitialPokeThenRun(startiter=2)
else:
    ctrl.initialCommand("ctrl.doSciRun()",freq=-1,startiter=2)
iscrn=science.iscrn.iscrn(None,ctrl.config,idstr="L0-2")
iatmos=science.iatmos.iatmos({"L0-2":iscrn},ctrl.config,idstr="b")
dm=science.xinterp_dm.dm(None,ctrl.config,idstr="dma")#this one (with no phase) for the widefield object (which adds the phase)
dm2=science.xinterp_dm.dm(None,ctrl.config,idstr="dmNFb")#this one for the science.
wf=science.wideField.WideField({"L0-2":iscrn,"dma":dm},ctrl.config,idstr="a")
c=science.wfscent.wfscent(wf,ctrl.config,idstr="acent")
r=science.tomoRecon.recon({"acent":c},ctrl.config,idstr="recon")
dm.newParent({"recon":r},"dma")
dm2.newParent({"recon":r,"atmos":iatmos},"dmNFb")
s=science.science.science(dm2,ctrl.config,idstr="b")
#save the solar images
save=base.saveOutput.saveOutput(wf,ctrl.config,idstr="solar")
#save the slopes
save2=base.saveOutput.saveOutput(c,ctrl.config,idstr="corrslopes")
nFieldX=ctrl.config.getVal("nFieldX")

execOrder=[iscrn,iatmos,dm,dm2,wf,c,r,s,save,save2]
#generate the truth, and save these...
for i in range(nFieldX):
    for j in range(nFieldX):
        execOrder.append(science.iatmos.iatmos({"L0-2":iscrn},ctrl.config,idstr="%d"%(i*nFieldX+j)))
        execOrder.append(science.wfscent.wfscent(execOrder[-1],ctrl.config,idstr="%d"%(i*nFieldX+j)))
        execOrder.append(base.saveOutput.saveOutput(execOrder[-1],ctrl.config,idstr="slopes%d"%(i*nFieldX+j)))
ctrl.mainloop(execOrder)

