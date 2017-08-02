import numpy
import science.infScrn
import science.physProp
import science.xinterp_dm
import science.wfscent
import science.tomoRecon
import science.science
import util.Ctrl
#Demonstration of physical optics propagation for the atmosphere.
#No AO, but many science directions.

ctrl=util.Ctrl.Ctrl(globals=globals())
nlayer=ctrl.config.getVal("nlayer")
nsci=ctrl.config.getVal("nsci")
scrnList=[]
lDict={}
#create the layers.
for i in range(nlayer):
    scrn=science.infScrn.infScrn(None,ctrl.config,idstr="L%d"%i)
    scrnList.append(scrn)
    lDict["L%d"%i]=scrn
execOrder=scrnList
for i in range(nsci):
    #propagate along a line of sight (Fresnel)
    execOrder.append(science.physProp.PhysProp(lDict,ctrl.config,idstr="sci%d"%(i+1)))
    #The science objects (Fourier, since from focal plane to image plane)
    execOrder.append(science.science.science(execOrder[-1],ctrl.config,idstr="sci%d"%(i+1)))


ctrl.mainloop(execOrder)
