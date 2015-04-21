#mpirun -np 1 -hostlist n1-c437  /usr/local/bin/mpipython $PWD/thisfile.py
#Python code created using the simulation setup GUI...
#Order of execution may not be quite optimal - you can always change by hand
#for large simulations - typically, the order of sends and gets may not be
#quite right.  Anyway, enjoy...
import util.Ctrl
import science.infScrn
import science.infAtmos
import science.xinterp_dm
import science.wfscent
import science.tomoRecon
import science.science
ctrl=util.Ctrl.Ctrl(globals=globals())
ctrl.doInitialPokeThenRun()
scrn=science.infScrn.infScrn(None,ctrl.config,idstr="L0")
atmos=science.infAtmos.infAtmos({"L0":scrn,},ctrl.config,idstr="a")
dm=science.xinterp_dm.dm(None,ctrl.config,idstr="dma")
wfs=science.wfscent.wfscent(dm,ctrl.config,idstr="a")
recon=science.tomoRecon.recon({"a":wfs,},ctrl.config,idstr="ngs")
sci=science.science.science(dm,ctrl.config,idstr="m")
dm.newParent({"1":atmos,"2":recon,},"dma")
ctrl.mainloop([scrn,atmos,dm,wfs,recon,sci])
print "Simulation finished..."


