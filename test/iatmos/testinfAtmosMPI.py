#mpirun -np 3 -hostlist n2-c437 n3-c437 n1-c437  /usr/local/bin/mpipython $PWD/thisfile.py
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
#import Scientific.MPI
import science.infAtmos
import science.infScrn
ctrl=util.Ctrl.Ctrl(globals=globals())
print "Rank %d imported modules"%ctrl.rank
#Set up the science modules...
newMPIGetList=[]
newMPISendList=[]
newSHMGetList=[]
newSHMSendList=[]
infAtmosList=[]
infScrnList=[]
#Add any personal code after this line and before the next, and it won't get overwritten
if ctrl.rank==1:
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L0").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,1,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L1").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,3,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L2").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,5,ctrl.mpiComm))
    infAtmosList.append(science.infAtmos.infAtmos({"L0":newMPIGetList[0],"L1":newMPIGetList[1],"L2":newMPIGetList[2],},ctrl.config,args={},idstr="b"))
    execOrder=[newMPIGetList[0],newMPIGetList[1],newMPIGetList[2],infAtmosList[0],]
    ctrl.mainloop(execOrder)
if ctrl.rank==2:
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L0"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[0],0,2,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[0],1,1,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L1"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[1],0,4,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[1],1,3,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L2"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[2],0,6,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[2],1,5,ctrl.mpiComm))
    execOrder=[infScrnList[0],newMPISendList[0],newMPISendList[1],infScrnList[1],newMPISendList[2],newMPISendList[3],infScrnList[2],newMPISendList[4],newMPISendList[5],]
    ctrl.mainloop(execOrder)
if ctrl.rank==0:
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L0").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,2,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L1").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,4,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L2").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,6,ctrl.mpiComm))
    infAtmosList.append(science.infAtmos.infAtmos({"L0":newMPIGetList[0],"L1":newMPIGetList[1],"L2":newMPIGetList[2],},ctrl.config,args={},idstr="a"))
    execOrder=[newMPIGetList[0],newMPIGetList[1],newMPIGetList[2],infAtmosList[0],]
    ctrl.mainloop(execOrder)
print "Simulation finished..."
#Add any personal code after this, and it will not get overwritten
#Scientific.MPI.world.abort(0)
ctrl.config.abort()

