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
import util.Ctrl
ctrl=util.Ctrl.Ctrl(globals=globals())
ctrl.initialCommand("wf.control['cal_source']=1",freq=-1,startiter=0)
ctrl.initialCommand("wf.control['cal_source']=0",freq=-1,startiter=1)
ctrl.initialCommand("c.newCorrRef();print 'Done new corr ref'",freq=-1,startiter=1)
ctrl.doInitialPokeThenRun(startiter=2)
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

execOrder=[iscrn,iatmos,dm,dm2,wf,c,r,s]
ctrl.mainloop(execOrder)

