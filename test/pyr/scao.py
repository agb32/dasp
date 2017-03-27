#mpirun -np 1 -hostlist n1-c437  /usr/local/bin/mpipython $PWD/thisfile.py
#Python code created using the simulation setup GUI...
#Order of execution may not be quite optimal - you can always change by hand
#for large simulations - typically, the order of sends and gets may not be
#quite right.  Anyway, enjoy...
import numpy
import util.Ctrl
import base.mpiGet
import base.mpiSend
import base.shmGet
import base.shmSend
import science.iscrn
import science.iatmos
import science.xinterp_dm
import science.pyramid
import science.tomoRecon
import science.science
ctrl=util.Ctrl.Ctrl(globals=globals())
print "Rank %d imported modules"%ctrl.rank
#Set up the science modules...
newMPIGetList=[]
newMPISendList=[]
newSHMGetList=[]
newSHMSendList=[]
iscrnList=[]
iatmosList=[]
dmList=[]
pyramidList=[]
reconList=[]
scienceList=[]
if not "nopoke" in ctrl.userArgList:ctrl.doInitialPokeThenRun()
#Add any personal code after this line and before the next, and it won't get overwritten
if ctrl.rank==0:
    iscrnList.append(science.iscrn.iscrn(None,ctrl.config,args={},idstr="allLayers"))
    iatmosList.append(science.iatmos.iatmos({"allLayers":iscrnList[0],},ctrl.config,args={},idstr="1"))
    dmList.append(science.xinterp_dm.dm(None,ctrl.config,args={},idstr="dm1"))
    pyramidList.append(science.pyramid.Pyramid(dmList[0],ctrl.config,args={},idstr="1"))
    reconList.append(science.tomoRecon.recon({"1":pyramidList[0],},ctrl.config,args={},idstr="recon"))
    scienceList.append(science.science.science(dmList[0],ctrl.config,args={},idstr="sci1"))
    dmList[0].newParent({"1":iatmosList[0],"2":reconList[0],},"dm1")
    execOrder=[iscrnList[0],iatmosList[0],dmList[0],pyramidList[0],reconList[0],scienceList[0],]
    ctrl.mainloop(execOrder)
print "Simulation finished..."
#Add any personal code after this, and it will not get overwritten
