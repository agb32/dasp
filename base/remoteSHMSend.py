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

import types
import cmod.shmem,base.fwdMsg
print "remoteSHMSend is depreciated.  Please use shmSend instead"
class remoteSHMSend:
    """
    DEPRECIATED
    Class to implement the sending end of shared memory communications.
    If you plan to alter this, then it needs some more thought when it
    comes to sending shm data - block before writing the data.
    When a splitter is needed::
    
       P
       |
       o   - SHMSend (can act as a splitter)
      / \  - SHM connections here
     o   o - SHMGet
     |   |
     C   C

    With a joiner is needed::
    
     P   P
     |   |
     o   o - SHMSend
     |   | - SHM connections here
     o   o - SHMGet
      \ / 
       o - joiner object
       |
       C
      
    Class variables (important to simulation programmer):
     - parent - object, the predecessor (parent) science object
     - shmInfo - SHMInfo instance, providing information about the shared memory connection
     - useFwdMsg - int flag, whether to use fwdMsg passing
    @cvar parent: the predecessor (parent) science object
    @type parent: object
    @cvar shmInfo: information about shared memory connection
    @type shmInfo: SHMInfo instance
    @cvar useFwdMsg: whether fwdMsg passing is used
    @type useFwdMsg: Int
    """
    def __init__(self,parent,shmInfo,args={}):
        """Initialise the SHM connection object.
        @param parent: predecessor (parent) object for retrieving science data
        @type parent: Instance
        @param shmInfo: shared memory information
        @type shmInfo: SHMInfo
        """
        #Lock object for fwdMsg can be a standard lock (not rwlock).
        #Lock for data needs to be an array of bit counts.
        self.objID=self.__module__.split(".")[-1]
        if args.has_key("idstr"):
            self.objID=self.objID+"_"+args["idstr"]
        self.parent=parent
        self.shmInfo=shmInfo
        self.useFwdMsg=shmInfo.useFwdMsg

    def next(self,childName,msg=None):
        """Obtain data from the parent object, and place in the shared memory
        array.\n
        First wait until all children have finished reading the data and get
        the write lock.  Note, the write lock can only be obtained once all
        children have requested and freed the read lock.
        Then, call parent.next with the fwdMsgList.  Then write the returned
        data to the shared memory array and free the writelock.
        @param msg: None
        @type msg: None
        @return: The data computed by parent object
        @rtype: None or Numeric array.
        """
        if msg!=None:
            print "ERROR: remoteSHMSend received msg in call to next()",msg
            raise "ERROR: remoteSHMSend received msg in call to next()"
        print "remoteSHMSend: getting writelock"
        self.shmInfo.getWriteLock()#blocks waiting for childs to read data. ie waits for all bitCnt=0.
        print "remoteSHMSend: writelock obtained"
        l=None
        if self.useFwdMsg:
            l=[]
            for key in self.shmInfo.fwdMsgDict.keys():
                fms=self.shmInfo.fwdMsgDict[key]
                f=fms.fwdMsg
                cmod.shmem.getlock(fms.semid)
                f.fromArray()
                l.append(f)
        data=self.parent.next(self.objID,l)
        if self.useFwdMsg:
            for key in self.shmInfo.fwdMsgDict.keys():
                cmod.shmem.freelock(self.shmInfo.fwdMsgDict[key].semid)
        if type(data)==types.NoneType:
            self.shmInfo.setNoneFlag()
        else:
            self.shmInfo.freeNoneFlag()
            self.copyToSHM(data)
        print "remoteSHMSend: freeing writelock"
        self.shmInfo.freeWriteLock()#set bitCnts to 1.
        print "remoteSHMSend: writeock freed"
        return data#not really needed...infact could be dangerous to use (lock)

    def copyToSHM(self,data):
        """Copy data into the shared memory, ready for passing to clients.
        @param data: The data to be shared
        @type data: None or Numeric array
        """
        if data is not self.shmInfo.data:#may already be the same objects...
            if type(self.shmInfo.data)==types.NoneType:
                self.shmInfo.dims=data.shape
                self.shmInfo.dtype=data.typecode()
                self.shmInfo.openShmFile()#open the shm file
            self.shmInfo.data[:]=data#copy data to shm region.
        else:#data is already in the shared memory array.
            pass
    
class fwdMsgStruct:
    """A structure for holding a shared memory fwdmsg instance.
    @cvar fwdMsg: The fwdMsg object
    @type fwdMsg: fwdMsg
    @cvar semid: Semaphore identifier for process synchronisation
    @type semid: Int
    @cvar name: Identifier for the shared memory region
    @type name: String
    """
    def __init__(self,fm,semid,name):
        """Initialise the structure.
        @param fm: The fwdMsg object
        @type fm: fwdMsg
        @param semid: The semaphore identifier
        @type semid: Integer
        @param name: SHM identifier
        @type name: String
        """
        self.fwdMsg=fm
        self.semid=semid
        self.name=name
    def __del__(self):
        print "Destructing fwdMsgStruct (semid)",self.semid
        cmod.shmem.semdel(self.semid)
        
class SHMInfo:
    """A class holding information about the shared memory arrays, and used
    for creating and removing them.

    Class variables (important to simulation programmer):
     - name - string, the identifier for the data SHM array, should commence with /
     - dims - tuple, the dimensions of the data SHM array
     - dtype - char, the datatype of the data SHM array
     - nChild - int, the number of children reading the shared memory array
     - useFwdMsg - int, flag whether fwdMsg passing is to be used

    Class variables (not important to simulation programmer):
     - semid - int, the semaphore identifier, used for process synchronisation
     - fwdMsgDict - dict of fwdMsgStruct, holding information about fwdmsg shared memory arrays

    @cvar name: Identifier for the data SHM array, should commence with /
    @type name: String
    @cvar dims: Dimensions of the data SHM array
    @type dims: Tuple of Int
    @cvar dtype: Datatype for data SHM array
    @type dtype: Char
    @cvar nChild: Number of children to read the SHM array
    @type nChild: Int
    @cvar useFwdMsg: Flag, whether fwdMsg passing is used
    @type useFwdMsg: Int
    @cvar semid: Semaphore identifier used for process synchronisation
    @type semid: Int
    @cvar fwdMsgDict: Details about fwdMsg shared memory arrays.
    @type fwdMsgDict: Dict of fwdMsgStruct
    """
    def __init__(self,name,dims,dtype,nChild,useFwdMsg=1):
        """Initialise the SHMInfo.
        @param name: Identifier for the data SHM array, should commence with /
        @type name: String
        @param dims: Dimensions of the data SHM array
        @type dims: Tuple of Int
        @param dtype: Datatype for data SHM array
        @type dtype: Char
        @param nChild: Number of children to read the SHM array
        @type nChild: Int
        @param useFwdMsg: Flag, whether fwdMsg passing is used
        @type useFwdMsg: Int
        """
        self.data=None
        self.name=name
        self.dims=dims
        self.dtype=dtype
        self.nChild=nChild
        self.useFwdMsg=useFwdMsg
        self.semid=None
        self.openShmFile()
        self.fwdMsgDict={}
        if self.useFwdMsg:
            self.openFwdMsgArray()
        self.initDataLock()
    def __del__(self):
        print "Destructing SHMInfo object",self.semid
        cmod.shmem.semdel(self.semid)
        self.closeShmFile()
        if self.useFwdMsg:
            self.closeFwdMsgArray()
    def closeShmFile(self):
        """Close the data SHM array"""
        cmod.shmem.unmap(self.data)
        self.data=None
        cmod.shmem.unlink(self.name)
    def openShmFile(self):
        """Open the data SHM array (closing first if already open)"""
        if self.data!=None:
            self.closeShmFile()
        self.data=cmod.shmem.open(self.dims,self.dtype,self.name,1)
        self.data[:]=99#initialise memory to something stupid.
    def closeFwdMsgArray(self):
        """Close the shared memory regions used for the fwdMsg (one per
        child process)."""
        if self.fwdMsgDict!=None:
            for key in self.fwdMsgDict.keys():
                f=self.fwdMsgDict[key].fwdMsg
                cmod.shmem.unmap(f.arr)
                cmod.shmem.unlink(key)
        self.fwdMsgDict={}
    def openFwdMsgArray(self):
        """Open the shared memory regions used for the fwdMsg (one per
        child process)."""
        self.closeFwdMsgArray()
        for i in range(self.nChild):
            name=self.name+"%04d"%i
            f=base.fwdMsg.fwdMsg()
            f.arr=cmod.shmem.open((28,),"1",name,1)
            semid=cmod.shmem.newlock(name=name,setone=1)
            self.fwdMsgDict[name]=fwdMsgStruct(f,semid,name)

    def initDataLock(self):
        """Initialise the read/write lock for the SHM array.  The initial state
        depends on the value of useFwdMsg.
        Should be called only once.
        """
        self.semid=cmod.shmem.newsemid(name="/dev/shm"+self.name,existingFile=1,nSems=self.nChild*2+2)#+1 is for the none type flag and +2 is for the initialised flag (ie so remoteSHMGet instances don't continue until this flag has been set)...
        print "remoteSHMSend: Initialised semid for data lock (semid=%d)"%self.semid
        if self.useFwdMsg:#start without the write lock
            for i in range(self.nChild):#set all bitCnts to 1.
                cmod.shmem.initSemaphore(self.semid,i,1)
            for i in range(self.nChild,self.nChild*2):#set all bits to zero
                cmod.shmem.initSemaphore(self.semid,i,0)
        else:#start with the write lock.
            for i in range(self.nChild):#set all bitCnts to 0.
                cmod.shmem.initSemaphore(self.semid,i,0)
            for i in range(self.nChild,self.nChild*2):#set all bits to zero
                cmod.shmem.initSemaphore(self.semid,i,1)

            
        self.freeNoneFlag()
        cmod.shmem.initSemaphore(self.semid,self.nChild*2+1,self.nChild)
        #cmod.shmem.initSemaphore(self.semid,self.nChild*2,0)#not a None type.
    def getWriteLock(self):
        """waits until the data write lock can be obtained and all readers have read the data (for bitcnt==0)"""
        a=range(self.nChild)
        b=[0]*self.nChild
        cmod.shmem.semop(self.semid,a,b)#wait for lower half to be set to zero
        #shmem.semop(self.semid,range(self.nChild),[0]*self.nChild)
    def freeWriteLock(self):
        """Frees the SHM data write lock (set all bitcnt to 1)."""
        cmod.shmem.semop(self.semid,range(self.nChild*2),[1]*self.nChild+[-1]*self.nChild)#set lower half to one, and decrement upper half.
    def setNoneFlag(self):#data is None
        """Data from predecessor (parent) is None, sets flag accordingly"""
        cmod.shmem.initSemaphore(self.semid,self.nChild*2,1)
    def freeNoneFlag(self):#data is not None (ie is proper data)
        """Data from predecessor (parent) is Numeric array, sets flag accordingly"""
        cmod.shmem.initSemaphore(self.semid,self.nChild*2,0)
