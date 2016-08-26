#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
    dims,dtype=science.iscrn.iscrn(None,ctrl.config,args={},forGUISetup=1,idstr="L0-34").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,0,1,ctrl.mpiComm))
    iatmosList.append(science.iatmos.iatmos({"L0-34":newMPIGetList[0],},ctrl.config,args={},idstr="b"))
    execOrder=[newMPIGetList[0],iatmosList[0],]
    ctrl.mainloop(execOrder)
if ctrl.rank==0:
    iscrnList.append(science.iscrn.iscrn(None,ctrl.config,args={},idstr="L0-34"))
    newMPISendList.append(base.mpiSend.newMPISend(iscrnList[0],1,2,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(iscrnList[0],2,1,ctrl.mpiComm))
    execOrder=[iscrnList[0],newMPISendList[0],newMPISendList[1],]
    ctrl.mainloop(execOrder)
if ctrl.rank==1:
    dims,dtype=science.iscrn.iscrn(None,ctrl.config,args={},forGUISetup=1,idstr="L0-34").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,0,2,ctrl.mpiComm))
    iatmosList.append(science.iatmos.iatmos({"L0-34":newMPIGetList[0],},ctrl.config,args={},idstr="a"))
    execOrder=[newMPIGetList[0],iatmosList[0],]
    ctrl.mainloop(execOrder)
print "Simulation finished..."
#Add any personal code after this, and it will not get overwritten
