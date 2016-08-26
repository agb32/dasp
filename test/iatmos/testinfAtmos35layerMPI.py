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
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L3").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,7,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L4").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,9,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L5").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,11,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L6").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,13,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L7").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,15,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L8").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,17,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L9").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,19,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L10").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,21,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L11").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,23,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L12").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,25,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L13").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,27,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L14").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,29,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L15").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,31,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L16").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,33,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L17").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,35,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L18").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,37,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L19").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,39,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L20").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,41,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L21").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,43,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L22").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,45,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L23").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,47,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L24").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,49,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L25").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,51,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L26").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,53,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L27").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,55,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L28").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,57,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L29").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,59,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L30").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,61,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L31").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,63,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L32").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,65,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L33").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,67,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L34").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,69,ctrl.mpiComm))
    infAtmosList.append(science.infAtmos.infAtmos({"L0":newMPIGetList[0],"L1":newMPIGetList[1],"L2":newMPIGetList[2],"L3":newMPIGetList[3],"L4":newMPIGetList[4],"L5":newMPIGetList[5],"L6":newMPIGetList[6],"L7":newMPIGetList[7],"L8":newMPIGetList[8],"L9":newMPIGetList[9],"L10":newMPIGetList[10],"L11":newMPIGetList[11],"L12":newMPIGetList[12],"L13":newMPIGetList[13],"L14":newMPIGetList[14],"L15":newMPIGetList[15],"L16":newMPIGetList[16],"L17":newMPIGetList[17],"L18":newMPIGetList[18],"L19":newMPIGetList[19],"L20":newMPIGetList[20],"L21":newMPIGetList[21],"L22":newMPIGetList[22],"L23":newMPIGetList[23],"L24":newMPIGetList[24],"L25":newMPIGetList[25],"L26":newMPIGetList[26],"L27":newMPIGetList[27],"L28":newMPIGetList[28],"L29":newMPIGetList[29],"L30":newMPIGetList[30],"L31":newMPIGetList[31],"L32":newMPIGetList[32],"L33":newMPIGetList[33],"L34":newMPIGetList[34],},ctrl.config,args={},idstr="b"))
    execOrder=[newMPIGetList[0],newMPIGetList[1],newMPIGetList[2],newMPIGetList[3],newMPIGetList[4],newMPIGetList[5],newMPIGetList[6],newMPIGetList[7],newMPIGetList[8],newMPIGetList[9],newMPIGetList[10],newMPIGetList[11],newMPIGetList[12],newMPIGetList[13],newMPIGetList[14],newMPIGetList[15],newMPIGetList[16],newMPIGetList[17],newMPIGetList[18],newMPIGetList[19],newMPIGetList[20],newMPIGetList[21],newMPIGetList[22],newMPIGetList[23],newMPIGetList[24],newMPIGetList[25],newMPIGetList[26],newMPIGetList[27],newMPIGetList[28],newMPIGetList[29],newMPIGetList[30],newMPIGetList[31],newMPIGetList[32],newMPIGetList[33],newMPIGetList[34],infAtmosList[0],]
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
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L3"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[3],0,8,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[3],1,7,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L4"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[4],0,10,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[4],1,9,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L5"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[5],0,12,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[5],1,11,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L6"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[6],0,14,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[6],1,13,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L7"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[7],0,16,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[7],1,15,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L8"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[8],0,18,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[8],1,17,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L9"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[9],0,20,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[9],1,19,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L10"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[10],0,22,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[10],1,21,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L11"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[11],0,24,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[11],1,23,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L12"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[12],0,26,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[12],1,25,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L13"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[13],0,28,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[13],1,27,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L14"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[14],0,30,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[14],1,29,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L15"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[15],0,32,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[15],1,31,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L16"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[16],0,34,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[16],1,33,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L17"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[17],0,36,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[17],1,35,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L18"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[18],0,38,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[18],1,37,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L19"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[19],0,40,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[19],1,39,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L20"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[20],0,42,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[20],1,41,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L21"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[21],0,44,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[21],1,43,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L22"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[22],0,46,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[22],1,45,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L23"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[23],0,48,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[23],1,47,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L24"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[24],0,50,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[24],1,49,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L25"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[25],0,52,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[25],1,51,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L26"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[26],0,54,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[26],1,53,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L27"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[27],0,56,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[27],1,55,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L28"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[28],0,58,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[28],1,57,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L29"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[29],0,60,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[29],1,59,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L30"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[30],0,62,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[30],1,61,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L31"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[31],0,64,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[31],1,63,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L32"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[32],0,66,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[32],1,65,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L33"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[33],0,68,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[33],1,67,ctrl.mpiComm))
    infScrnList.append(science.infScrn.infScrn(None,ctrl.config,args={},idstr="L34"))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[34],0,70,ctrl.mpiComm))
    newMPISendList.append(base.mpiSend.newMPISend(infScrnList[34],1,69,ctrl.mpiComm))
    execOrder=[infScrnList[0],newMPISendList[0],newMPISendList[1],infScrnList[1],newMPISendList[2],newMPISendList[3],infScrnList[2],newMPISendList[4],newMPISendList[5],infScrnList[3],newMPISendList[6],newMPISendList[7],infScrnList[4],newMPISendList[8],newMPISendList[9],infScrnList[5],newMPISendList[10],newMPISendList[11],infScrnList[6],newMPISendList[12],newMPISendList[13],infScrnList[7],newMPISendList[14],newMPISendList[15],infScrnList[8],newMPISendList[16],newMPISendList[17],infScrnList[9],newMPISendList[18],newMPISendList[19],infScrnList[10],newMPISendList[20],newMPISendList[21],infScrnList[11],newMPISendList[22],newMPISendList[23],infScrnList[12],newMPISendList[24],newMPISendList[25],infScrnList[13],newMPISendList[26],newMPISendList[27],infScrnList[14],newMPISendList[28],newMPISendList[29],infScrnList[15],newMPISendList[30],newMPISendList[31],infScrnList[16],newMPISendList[32],newMPISendList[33],infScrnList[17],newMPISendList[34],newMPISendList[35],infScrnList[18],newMPISendList[36],newMPISendList[37],infScrnList[19],newMPISendList[38],newMPISendList[39],infScrnList[20],newMPISendList[40],newMPISendList[41],infScrnList[21],newMPISendList[42],newMPISendList[43],infScrnList[22],newMPISendList[44],newMPISendList[45],infScrnList[23],newMPISendList[46],newMPISendList[47],infScrnList[24],newMPISendList[48],newMPISendList[49],infScrnList[25],newMPISendList[50],newMPISendList[51],infScrnList[26],newMPISendList[52],newMPISendList[53],infScrnList[27],newMPISendList[54],newMPISendList[55],infScrnList[28],newMPISendList[56],newMPISendList[57],infScrnList[29],newMPISendList[58],newMPISendList[59],infScrnList[30],newMPISendList[60],newMPISendList[61],infScrnList[31],newMPISendList[62],newMPISendList[63],infScrnList[32],newMPISendList[64],newMPISendList[65],infScrnList[33],newMPISendList[66],newMPISendList[67],infScrnList[34],newMPISendList[68],newMPISendList[69],]
    ctrl.mainloop(execOrder)
if ctrl.rank==0:
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L0").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,2,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L1").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,4,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L2").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,6,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L3").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,8,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L4").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,10,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L5").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,12,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L6").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,14,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L7").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,16,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L8").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,18,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L9").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,20,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L10").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,22,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L11").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,24,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L12").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,26,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L13").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,28,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L14").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,30,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L15").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,32,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L16").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,34,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L17").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,36,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L18").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,38,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L19").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,40,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L20").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,42,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L21").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,44,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L22").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,46,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L23").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,48,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L24").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,50,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L25").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,52,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L26").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,54,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L27").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,56,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L28").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,58,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L29").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,60,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L30").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,62,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L31").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,64,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L32").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,66,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L33").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,68,ctrl.mpiComm))
    dims,dtype=science.infScrn.infScrn(None,ctrl.config,args={},forGUISetup=1,idstr="L34").outputData
    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,2,70,ctrl.mpiComm))
    infAtmosList.append(science.infAtmos.infAtmos({"L0":newMPIGetList[0],"L1":newMPIGetList[1],"L2":newMPIGetList[2],"L3":newMPIGetList[3],"L4":newMPIGetList[4],"L5":newMPIGetList[5],"L6":newMPIGetList[6],"L7":newMPIGetList[7],"L8":newMPIGetList[8],"L9":newMPIGetList[9],"L10":newMPIGetList[10],"L11":newMPIGetList[11],"L12":newMPIGetList[12],"L13":newMPIGetList[13],"L14":newMPIGetList[14],"L15":newMPIGetList[15],"L16":newMPIGetList[16],"L17":newMPIGetList[17],"L18":newMPIGetList[18],"L19":newMPIGetList[19],"L20":newMPIGetList[20],"L21":newMPIGetList[21],"L22":newMPIGetList[22],"L23":newMPIGetList[23],"L24":newMPIGetList[24],"L25":newMPIGetList[25],"L26":newMPIGetList[26],"L27":newMPIGetList[27],"L28":newMPIGetList[28],"L29":newMPIGetList[29],"L30":newMPIGetList[30],"L31":newMPIGetList[31],"L32":newMPIGetList[32],"L33":newMPIGetList[33],"L34":newMPIGetList[34],},ctrl.config,args={},idstr="a"))
    execOrder=[newMPIGetList[0],newMPIGetList[1],newMPIGetList[2],newMPIGetList[3],newMPIGetList[4],newMPIGetList[5],newMPIGetList[6],newMPIGetList[7],newMPIGetList[8],newMPIGetList[9],newMPIGetList[10],newMPIGetList[11],newMPIGetList[12],newMPIGetList[13],newMPIGetList[14],newMPIGetList[15],newMPIGetList[16],newMPIGetList[17],newMPIGetList[18],newMPIGetList[19],newMPIGetList[20],newMPIGetList[21],newMPIGetList[22],newMPIGetList[23],newMPIGetList[24],newMPIGetList[25],newMPIGetList[26],newMPIGetList[27],newMPIGetList[28],newMPIGetList[29],newMPIGetList[30],newMPIGetList[31],newMPIGetList[32],newMPIGetList[33],newMPIGetList[34],infAtmosList[0],]
    ctrl.mainloop(execOrder)
print "Simulation finished..."
#Add any personal code after this, and it will not get overwritten
