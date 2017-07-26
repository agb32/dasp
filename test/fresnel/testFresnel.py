import numpy
import science.infScrn
import science.physProp
import science.xinterp_dm
import science.wfscent
import science.tomoRecon
import science.science
import util.Ctrl
#Demonstration of physical optics propagation for the atmosphere.

ctrl=util.Ctrl.Ctrl(globals=globals())
ctrl.doInitialPokeThenRun()
nlayer=ctrl.config.getVal("nlayer")
scrnList=[]
lDict={}
#create the layers.
for i in range(nlayer):
    scrn=science.infScrn.infScrn(None,ctrl.config,idstr="L%d"%i)
    scrnList.append(scrn)
    lDict["L%d"%i]=scrn
#propagate along a line of sight (Fresnel)
fres=science.physProp.PhysProp(lDict,ctrl.config,idstr="1")
#the DM (ground conjugate)
dm=science.xinterp_dm.dm(None,ctrl.config,idstr="dm1")
#The WFS (Fourier, since from focal plane to image plane)
wfs=science.wfscent.wfscent(dm,ctrl.config,args={},idstr="1")
#the reconstructor
recon=science.tomoRecon.recon({"1":wfs,},ctrl.config,idstr="recon")
#The science object (Fourier, since from focal plane to image plane)
sci=science.science.science(dm,ctrl.config,idstr="sci1")
#set the feedback loop:
dm.newParent({"1":fres,"2":recon,},"dm1")


execOrder=scrnList+[fres,dm,wfs,recon,sci]

ctrl.mainloop(execOrder)
