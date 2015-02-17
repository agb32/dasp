"""Tests the new iscrn module.  
Note - this may not be kept up to date with the module."""

import science.iscrn
import science.iatmos
import base.readConfig
import util.Ctrl
ctrl=util.Ctrl.Ctrl(globals=globals())
iscrn01=science.iscrn.iscrn(None,ctrl.config,idstr="L0-1")
iscrn2=science.iscrn.iscrn(None,ctrl.config,idstr="L2")
iatm=science.iatmos.iatmos({"L0-1":iscrn01,"L2":iscrn2},ctrl.config,idstr="a")

execOrder=[iscrn01,iscrn2,iatm]
ctrl.mainloop(execOrder)

