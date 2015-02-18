#mpirun -np 3 -hostlist n1-c437 n2-c437 n3-c437  /usr/local/bin/mpipython $PWD/thisfile.py
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
import science.iscrn
import science.iatmos
ctrl=util.Ctrl.Ctrl(globals=globals())
print "Rank %d imported modules"%ctrl.rank
#Set up the science modules...
newMPIGetList=[]
newMPISendList=[]
newSHMGetList=[]
newSHMSendList=[]
iscrnList=[]
iatmosList=[]
#Add any personal code after this line and before the next, and it won't get overwritten
if ctrl.rank==2:
    dims,dtype=science.iscrn.iscrn(None,ctrl.config,args={},forGUISetup=1,idstr="L0-1").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,0,1,ctrl.mpiComm))
    dims,dtype=science.iscrn.iscrn(None,ctrl.config,args={},forGUISetup=1,idstr="L2").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,0,3,ctrl.mpiComm))
    iatmosList.append(science.iatmos.iatmos({"L0-1":newMPIGetList[0],"L2":newMPIGetList[1],},ctrl.config,args={},idstr="b"))
    execOrder=[newMPIGetList[0],newMPIGetList[1],iatmosList[0],]
    ctrl.mainloop(execOrder)
if ctrl.rank==0:
    iscrnList.append(science.iscrn.iscrn(None,ctrl.config,args={},idstr="L0-1"))
    newMPISendList.append(base.mpiSend.newMPISend(iscrnList[0],1,2,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(iscrnList[0],2,1,ctrl.mpiComm))
    iscrnList.append(science.iscrn.iscrn(None,ctrl.config,args={},idstr="L2"))
    newMPISendList.append(base.mpiSend.newMPISend(iscrnList[1],1,4,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(iscrnList[1],2,3,ctrl.mpiComm))
    execOrder=[iscrnList[0],newMPISendList[0],newMPISendList[1],iscrnList[1],newMPISendList[2],newMPISendList[3],]
    ctrl.mainloop(execOrder)
if ctrl.rank==1:
    dims,dtype=science.iscrn.iscrn(None,ctrl.config,args={},forGUISetup=1,idstr="L0-1").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,0,2,ctrl.mpiComm))
    dims,dtype=science.iscrn.iscrn(None,ctrl.config,args={},forGUISetup=1,idstr="L2").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,0,4,ctrl.mpiComm))
    iatmosList.append(science.iatmos.iatmos({"L0-1":newMPIGetList[0],"L2":newMPIGetList[1],},ctrl.config,args={},idstr="a"))
    execOrder=[newMPIGetList[0],newMPIGetList[1],iatmosList[0],]
    ctrl.mainloop(execOrder)
print "Simulation finished..."
#Add any personal code after this, and it will not get overwritten
#Scientific.MPI.world.abort(0)
ctrl.config.abort()

