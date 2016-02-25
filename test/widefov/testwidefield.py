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

