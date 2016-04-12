"""Tests the new iscrn module.  
Note - this may not be kept up to date with the module."""
import science.iscrn
import science.wideField
#import science.iatmos
import base.readConfig
import util.Ctrl
ctrl=util.Ctrl.Ctrl(globals=globals())
iscrn=science.iscrn.iscrn(None,ctrl.config,idstr="L0-2")
wf=science.wideField.WideField({"L0-2":iscrn},ctrl.config,idstr="a")
#atmos=science.iatmos.iatmos({"L0-2":iscrn},ctrl.config,idstr="a")
execOrder=[iscrn,wf]#,atmos]
ctrl.mainloop(execOrder)

