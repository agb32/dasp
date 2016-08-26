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

import cmod.shmem,os,types,time,base.aobase
import numpy
import sys
#see also, shmGet2 below...
debugConnections=1
class shmGet:
    """
    Class to implement the receiving end of shared memory communications.
    Gets data from SHM array and return, blocking until the data is ready to
    be returned.  Returns the actual shared memory array and has a readlock
    on it.
    This should be called by an object whos science parent (predecessor) is
    communicating via SHM.

    Get a normal lock (not rwlock) on fwdMsg.  Use bitcnt to get locks on data.

    Class variables (important to simulation programmer):
     - shmParent - Object, storing information about the shared memory.
     - objID - string, identifier
     - debug - anything, printed with debug information
    Class variables (not important to simulation programmer):
     - generate - whether to do this iteration or not
     - dataValid - whether valid this iteration or not
     - newDataWaiting - whether to expect new data this iteration.
    @cvar shmParent: Details of the SHM memory
    @type shmParent: SHMParent object
    @cvar objID: Object identifier
    @type objID: String
    @cvar debug: If not None, will allow printing of debug messages
    @type debug: Anything
    @cvar generate: Whether to generate
    @type generate: Int
    @cvar dataValid: Whether data is valid
    @type dataValid: Int
    @cvar newDataWaiting: Whether expecting new data
    @type newDataWaiting: Int
    @cvar makeCopy: Whether to make copy of data
    @type makeCopy: Int
    """
    def __init__(self,shmParent,makeCopy=0,args={},debug=None):
        """Initialise the class.
        @param shmParent: SHM memory description
        @type shmParent: SHMParent object
        @param makeCopy: Whether to make copy of data
        @type makeCopy: Int
        @param args: Additional arguments, at present only idstr
        @type args: String
        @param debug: If not None, will allow printing of debug messages
        @type debug: Anything
        @param generate: Whether to generate
        @type generate: Int
        @param dataValid: Whether data is valid
        @type dataValid: Int
        @param newDataWaiting: Whether expecting new data
        @type newDataWaiting: Int
        """
        print "please use shmGet.shmGet2 or shmGet.newShmGet instead."
        self.objID=self.__module__.split(".")[-1]
        if args.has_key("idstr"):
            self.objID=self.objID+"_"+args["idstr"]
        self.debug=debug
        self.shmParent=shmParent
        #self.makeCopy=makeCopy#make a copy of the data before returning it?
        self.generate=1
        self.dataValid=0
        self.newDataWaiting=1
        self.makeCopy=makeCopy
        #self.firstIter=1
        #the data write lock earlier.
        
    def generateNext(self,msg=None):
        """called by the child when it wants more data.
        Get the data and return it.\n

        The read lock is freed once this finishes reading it.  However, other
        modules may still be reading it then, so care must be taken.  This is taken
        care of by the writer having to wait until the reader starts reading again, before
        it can write.  Might be slightly wasteful way of doing it, but means
        safety.

        The use of thread locks is important here, to ensure data isn't
        overwritten while it is being read:\n
        We get a read lock on the data, once new data has arrived in the array.
        This prevents the remoteSHMSend object from writing to the array until
        we release the read lock.
        If making a local copy of the data, this is made.

        @param msg: The fwdMsg to pass to the science parent (or None)
        @type msg: None or fwdMsg object
        """

##         if self.useFwdMsg:
##             if type(self.shmParent.fwdMsgArr)==types.NoneType:
##                 self.shmParent.openFwdMsgArr()
##             cmod.shmem.getlock(self.shmParent.fwdMsgLock)
##             msg.toArray(self.shmParent.fwdMsgArr)#write msg to the SHM array.
##             cmod.shmem.freelock(self.shmParent.fwdMsgLock)

##         if (self.makeCopy==0):#not making a copy, so need to free the shm data.
##             if self.firstIter==0 or self.useFwdMsg==1:
##                 print "shmGet: freeing readlock"
##                 self.shmParent.freeReadLock()#free the lock from previous reads.  Set bitcnt to zero
##                 print "shmGet: readlock freed"
##             else:
##                 self.firstIter=0
##         elif self.makeCopy==1 and self.useFwdMsg==1:#using the fwdmsg, which has now been written, so free the lock...
##             print "shmGet: freeing readlock"
##             self.shmParent.freeReadLock()#free the lock from previous reads.  Set bitcnt to zero
##             print "shmGet: readlock freed"

        
        
        #this should now allow the writelock to be obtained by the writer...
        if self.debug!=None: print "shmGet: In generateNext (debug=%s)"%str(self.debug)
        self.shmParent.setGenerateFlag(self.generate)
        if self.debug!=None: print "shmGet: getting readlock (debug=%s)"%str(self.debug)
        self.shmParent.getReadLock()#blocks if data has already been read - until new data is ready... waits until bitcnt==1.
        if self.debug!=None: print "shmGet: readlock obtained (debug=%s)"%str(self.debug)
        if self.shmParent.NoneFlagIsSet():
            self.outputData=None
            self.dataValid=0
        else:
            self.outputData=self.shmParent.data
            if self.makeCopy:
                self.outputData=self.outputData.copy()
            self.dataValid=1
        self.shmParent.freeReadLock()

        if self.debug!=None: print "shmGet: Done generateNext, dataValid=%d (debug=%s)"%(self.dataValid,str(self.debug))
    def setGenerate(self,val):
        """Set the generate flag.
        @param val: Value to set generate to
        @type val: Int"""
        self.generate=val

class shmParent:
    """An object which creates the shared memory array into which the data will
    be placed.
    nReaders should be the number of processes intending
    to READ the data (if a process is only writing, is isn't counted as a
    reader).  ie the number of children.

    Class variables (important to simulation programmer):
     - name - string, identifier for the shared memory, should commence with /.
     - dims - tuple, dimensions for the shared memory
     - dtype - char, datatype for shared memory array, e.g. Numeric.Int32
     - readCountBit - int, unique identifier, between 0 and less than nReaders, used to determine which semaphore can be used.  Other SHMParent objects reading the same shared memory must have a different readCountBit.
     - nReaders - int, number of instances of SHMParent reading the data.
     - useFwdMsg - int, whether to use fwdMsg passing or not.
    Class variables (not important to simulation programmer):
     - data - array, the shared memory data.
     - fwdMsgFName - string, name of the fwdMsg shared memory array for this
     - semid - int, semaphore ID
     - fwdMsgLock - shmem Lock object, used to lock the fwdMsg array when writing.
    @cvar name: identifier for the shared memory, should commence with /.
    @type name: String
    @cvar dims: dimensions for the shared memory
    @type dims: Tuple
    @cvar dtype: datatype for shared memory array, e.g. Numeric.Int32
    @type dtype: Char
    @cvar readCountBit: unique identifier, between 0 and less than nReaders, used to determine which semaphore can be used.  Other SHMParent objects reading the same shared memory must have a different readCountBit.
    @type readCountBit: Int
    @cvar nReaders: number of instances of SHMParent reading the data.
    @type nReaders: Int
    @cvar data: Shared memory data array
    @type data: Array
    @cvar semid: Semaphore identifier used for process synchronisation
    @type semid: Int
    """
    def __init__(self,name,dims,dtype,readCountBit,nReaders):
        """
        Initialise the shared memory object.
        @param name: identifier for the shared memory, should commence with /.
        @type name: String
        @param dims: dimensions for the shared memory
        @type dims: Tuple
        @param dtype: datatype for shared memory array, e.g. Numeric.Int32
        @type dtype: Char
        @param readCountBit: unique identifier, between 0 and less than nReaders, used to determine which semaphore can be used.  Other SHMParent objects reading the same shared memory must have a different readCountBit.
        @type readCountBit: Int
        @param nReaders: number of instances of SHMParent reading the data.
        @type nReaders: Int
        """
        self.data=None
        self.name=name+"_"+os.environ["USER"]
        name=self.name
        #self.fwdMsgFName=name+"%04d"%readCountBit
        self.dims=dims
        if type(dtype)==type(""):#eg "f" for float32
            self.dtype=dtype
        else:#probably numpy.float32, nmpy.int32 etc.  Convert these to char...
            self.dtype=None
            if dtype in numpy.typeDict.values():
                for key in numpy.typeDict:
                    if type(key)==type("") and len(key)==1:
                        if numpy.typeDict[key]==dtype:
                            self.dtype=key
            if self.dtype==None:
                raise Exception("shmSend: dtype %s not recognised"%str(dtype))
        self.readCountBit=readCountBit
        self.nReaders=nReaders#number of readers of this shm.
        self.openShmFile()
        #self.shmLock=shmRWBLock(name,nReaders,readCountBit,0)
        self.semid=cmod.shmem.newsemid(name="/dev/shm"+self.name,existingFile=1,nSems=self.nReaders*2+4)
        print "shmGet: Waiting for initialisation of semid"
        cmod.shmem.semop(self.semid,nReaders*2+3,-1)#wait for the initialisation to take place (shmSend.py in initDataLock())
        #self.fwdMsgArr=None
        #self.fwdMsgLock=None
        print "shmGet: Parent object Initialisation complete"
    def __del__(self):
        """Don't think this deletion is needed"""
        print "Destructing shmGet - not deleting any semaphores"
##         print "Destructing SHMParent object (self.semid)",self.semid
##         cmod.shmem.semdel(self.semid)
##         print "Destructing SHMParent object (self.fwdMsgLock)",self.fwdMsgLock
##         cmod.shmem.semdel(self.fwdMsgLock)
    def closeShmFile(self):
        """Close the shared memory array"""
        cmod.shmem.unmap(self.data)
        self.data=None
        #self.name=None
        #if self.created==1:
        #    shmem.unlink(self.name)
        #self.created=0
    def openShmFile(self):
        """Open a shared memory array"""
        if self.data!=None:
            self.closeShmFile()
        #if self.create:
        #    self.created=1
        while type(self.data)==types.NoneType:
            try:
                self.data=cmod.shmem.open(self.dims,self.dtype,self.name,0)
            except:
                print "shmGet couldn't open SHM array %s %s %s.  Waiting 2s and retrying"%(self.name,str(self.dims),str(self.dtype))
                self.data=None
                time.sleep(2)
        print "shmGet: Opened shm file %s okay"%self.name
    def freeReadLock(self):
        """Free the readlock on the shared memory data."""
        #cmod.shmem.semop(self.semid,range(self.nReaders*2),[-1]*self.nReaders+[1]*self.nReaders)
        cmod.shmem.semop(self.semid,[self.readCountBit,self.nReaders+self.readCountBit],[-1,1])
    def getReadLock(self):
        """Waits until new data has been written to the shared memory array, and until the writing of this data is complete."""
        #wait for all bitcnt==1.
        #cmod.shmem.semop(self.semid,range(self.nReaders,self.nReaders*2),[0]*self.nReaders)
        cmod.shmem.semop(self.semid,self.nReaders+self.readCountBit,0)
        
##     def openFwdMsgArr(self):
##         """Opens a fwdMsg shared memory array"""
##         if self.fwdMsgArr!=None:
##             cmod.shmem.unmap(self.fwdMsgArr)
##             self.fwdMsgArr=None
##         try:
##             self.fwdMsgArr=cmod.shmem.open((28,),"1",self.fwdMsgFName,0)
##         except:
##             print "Couldn't open fwdMsgArr %s in shmGet.  Retrying in 2s."%self.fwdMsgFName
##             time.sleep(2)
##             self.fwdMsgArr=cmod.shmem.open((28,),"1",self.fwdMsgFName,0)
            
    def NoneFlagIsSet(self):
        """Checks whether there is data in the shared memory array, or whether
        the science successor (parent) returned None.
        @return: Flag representing 1 for None, or 0 for array data
        @rtype: Int"""
        return cmod.shmem.getSemValue(self.semid,self.nReaders*2)
    def setGenerateFlag(self,val):
        """Decrement one semaphore of the set (to show that we have done this
        part, and then either increment or leave alone another semaphore,
        depending on the value of val.
        @param val: Whether to generate data or not
        @type val: Int (0 or 1)
        """        
        if val:
            cmod.shmem.semop(self.semid,self.nReaders*2+2,1)
        cmod.shmem.semop(self.semid,self.nReaders*2+1,-1)


class shmGet2(base.aobase.aobase):
    """
    Class to implement the receiving end of shared memory communications.
    Gets data from SHM array and return, blocking until the data is ready to
    be returned.  Returns the actual shared memory array and has a readlock
    on it.
    This should be called by an object whos science parent (predecessor) is
    communicating via SHM.
    This is actually implemented using semephores:
    Get                   Send
    Init l2=1             Init l1=1
    Then:
                          Wait on l1
    Dec l1
    Wait on l2
                          Write data
                          l1=1
                          Dec l2
    l2=1
    """
    def __init__(self,name,dims,dtype,makeCopy=0,args={},debug=None,idstr=None):
        if debugConnections==1:
            if debug==None:
                debug="shmGet %s"%name
        base.aobase.aobase.__init__(self,None,None,args=args,debug=debug,idstr=idstr)
        self.makeCopy=makeCopy
        self.name=name+"_"+os.environ["USER"]
        self.dims=dims
        if type(dtype)==type(""):#eg "f" for float32
            self.dtype=dtype
        else:#probably numpy.float32, nmpy.int32 etc.  Convert these to char...
            self.dtype=None
            if dtype in numpy.typeDict.values():
                for key in numpy.typeDict:
                    if type(key)==type("") and len(key)==1:
                        if numpy.typeDict[key]==dtype:
                            self.dtype=key
            if self.dtype==None:
                raise Exception("shmSend: dtype %s not recognised"%str(dtype))
        self.data=None
        self.openShmFile()
        self.semid=cmod.shmem.newsemid(name="/dev/shm"+self.name,existingFile=1,nSems=3)
        time.sleep(0.1)
        os.unlink("/dev/shm"+self.name)
        print "shmGet: waiting for initialisation of semid %d for %s"%(self.semid,self.name)
        cmod.shmem.semop(self.semid,2,-1)
        print "shmGet: semaphore %d initialised"%self.semid
        if self.makeCopy:
            self.outputData=self.data.copy()
        else:
            self.outputData=self.data
    def openShmFile(self):
        if self.data!=None:
            self.closeShmFile()
        tdelay=1
        cnt=0
        while type(self.data)==type(None):
            try:
                self.data=cmod.shmem.open(self.dims,self.dtype,self.name,0,0)
            except:
                if cnt>=tdelay:
                    print "shmGet2 couldn't open SHM array %s %s %s.  Waiting 2s and retrying (%ds till next message)\r"%(self.name,str(self.dims),str(self.dtype),tdelay*2*2),
                    cnt=0
                    tdelay*=2
                sys.stdout.flush()
                self.data=None
                time.sleep(2)
                cnt+=1
        print "shmGet2: opened shm file %s okay"%self.name
    def closeShmFile(self):
        """Close the shared memory array"""
        cmod.shmem.unmap(self.data)
        self.data=None
    def NoneFlagIsSet(self):
        """Checks whether there is data in the shared memory array, or whether
        the science successor (parent) returned None.
        @return: Flag representing 1 for None, or 0 for array data
        @rtype: Int"""
        return cmod.shmem.getSemValue(self.semid,2)

    def setGenerate(self,val):
        """Set the generate flag.
        @param val: Value to set generate to
        @type val: Int"""
        self.generate=val
    def setParentGenerateFlag(self):
        if self.generate:
            cmod.shmem.initSemaphore(self.semid,2,1)
        else:
            cmod.shmem.initSemaphore(self.semid,2,0)

    def generateNext(self,msg=None):
        self.setParentGenerateFlag()#can the parent generate?
        cmod.shmem.semop(self.semid,0,1)#allow shmSend to decrement...
        cmod.shmem.semop(self.semid,1,-1)#wait until we can decrement.
        #now check data is valid or None...
        if self.NoneFlagIsSet():
            self.outputData=None
            self.dataValid=0
        else:
            self.outputData=self.data
            if self.makeCopy:
                self.outputData=self.outputData.copy()
            self.dataValid=1
        
def newSHMGet(name,dims,dtype,idstr=None):
    sg=shmGet2(name,dims,dtype,idstr=None)
    return sg
