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

import cmod.shmem,base.fwdMsg,types,time
print "remoteSHMGet is depreciated.  Please use shmGet instead"
class remoteSHMGet:
    """
    DEPRECIATED
    Class to implement the receiving end of shared memory communications.
    Gets data from SHM array and return, blocking until the data is ready to
    be returned.  Returns the actual shared memory array and has a readlock
    on it.
    This should be called by an object whos science parent (predecessor) is
    communicating via SHM.

    Get a normal lock (not rwlock) on fwdMsg.  Use bitcnt to get locks on data.

    Class variables (important to simulation programmer):
     - shmParent - Object, storing information about the shared memory.
     - makeCopy - Int flag, determining whether a copy of the data is made before returning it.
     - useFwdMsg - Int flag, determining whether FwdMsg passing is used.
    Class variables (not important to simulation programmer):
     - firstIter - Int flag, whether currently in first iteration or not.
    @cvar shmParent: Details of the SHM memory
    @type shmParent: SHMParent object
    @cvar makeCopy: Flag, whether to make a copy of the data before returning
    @type makeCopy: Int
    @cvar useFwdMsg: Flag, whether to use fwdMsg passing
    @type useFwdMsg: Int
    @cvar firstIter: Whether SHM has been used or not
    @type firstIter: Int
    """
    def __init__(self,shmParent,makeCopy=1,args={}):
        """Initialise the class.
        @param shmParent: SHM memory description
        @type shmParent: SHMParent object
        @param makeCopy: Whether to make a copy of the data
        @type makeCopy: Int
        """
        self.objID=self.__module__.split(".")[-1]
        if args.has_key("idstr"):
            self.objID=self.objID+"_"+args["idstr"]
        self.shmParent=shmParent
        self.makeCopy=makeCopy#make a copy of the data before returning it?
        self.useFwdMsg=shmParent.useFwdMsg#if this is zero, it can be possible to free
        self.firstIter=1
        #the data write lock earlier.
        
    def next(self,childName,msg=None):
        """called by the child when it wants more data.
        Get the data and return it.\n

        The use of thread locks is important here, to ensure data isn't
        overwritten while it is being read:\n
        We get a read lock on the data, once new data has arrived in the array.
        This prevents the remoteSHMSend object from writing to the array until
        we release the read lock.

        If making a local copy of the data, this is made.  If this is the
        case, and we're not using fwdMsg passing, the readlock is freed,
        since the SHMSend object is now free to write into the array.

        However, if we're not making a local copy of the data, or are using
        fwdMsg passing, we cannot free the readlock yet.

        The data is now returned to the science object that requested it.

        When the science object next requests data:
        If we haven't yet freed the readlock, we do so (ie if we're not making
        a local copy of the data, or are using fwdMsg passing).

        Thats basically it, with a few subtleties.  If you need to change this
        code, make sure you understand it fully, and test all configurations
        as you go along.

        @param msg: The fwdMsg to pass to the science parent (or None)
        @type msg: None or fwdMsg object
        @return: The data
        @rtype: None or Numeric array
        """
        print "remoteSHMGet: In next()"
        if type(msg)==types.NoneType:
            if self.useFwdMsg:
                msg=base.fwdMsg.fwdMsg()
        else:
            if self.useFwdMsg==0:
                print "ERROR:  NOT USING FWDMSG, BUT IT IS USED",msg
                raise Exception("ERROR:  NOT USING FWDMSG, BUT IT IS USED")
        if self.useFwdMsg:
            if type(self.shmParent.fwdMsgArr)==types.NoneType:
                self.shmParent.openFwdMsgArr()
            cmod.shmem.getlock(self.shmParent.fwdMsgLock)
            msg.toArray(self.shmParent.fwdMsgArr)#write msg to the SHM array.
            cmod.shmem.freelock(self.shmParent.fwdMsgLock)

        if (self.makeCopy==0):#not making a copy, so need to free the shm data.
            if self.firstIter==0 or self.useFwdMsg==1:
                print "remoteSHMGet: freeing readlock"
                self.shmParent.freeReadLock()#free the lock from previous reads.  Set bitcnt to zero
                print "remoteSHMGet: readlock freed"
            else:
                self.firstIter=0
        elif self.makeCopy==1 and self.useFwdMsg==1:#using the fwdmsg, which has now been written, so free the lock...
            print "remoteSHMGet: freeing readlock"
            self.shmParent.freeReadLock()#free the lock from previous reads.  Set bitcnt to zero
            print "remoteSHMGet: readlock freed"
            
        #this should now allow the writelock to be obtained by the writer...
        print "remoteSHMGet: getting readlock"
        self.shmParent.getReadLock()#blocks if data has already been read - until new data is ready... waits until bitcnt==1.
        print "remoteSHMGet: readlock obtained"
        if self.shmParent.NoneFlagIsSet():
            data=None
            if self.makeCopy==1 and self.useFwdMsg==0:
                print "remoteSHMGet: freeing readlock"
                self.shmParent.freeReadLock()
                print "remoteSHMGet: readlock freed"
        else:
            if self.makeCopy:
                data=self.shmParent.data.copy()
                if self.useFwdMsg==0:#Only free datalock if not expecting msg.
                    print "remoteSHMGet: freeing readlock"
                    self.shmParent.freeReadLock()
                    print "remoteSHMGet: readlock freed"
            else:#return the shm array (keep the read lock on it!)
                data=self.shmParent.data
        return data


class SHMParent:
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
    @cvar useFwdMsg: whether to use fwdmsg passing
    @type useFwdMsg: Int
    @cvar nReaders: number of instances of SHMParent reading the data.
    @type nReaders: Int
    @cvar data: Shared memory data array
    @type data: Array
    @cvar fwdMsgFName: fwdMsg shared memory identifier
    @type fwdMsgFName: String
    @cvar semid: Semaphore identifier used for process synchronisation
    @type semid: Int
    @cvar fwdMsgLock: Lock used for process synchronisation
    @type fwdMsgLock: shmem lock object
    """
    def __init__(self,name,dims,dtype,readCountBit,nReaders,useFwdMsg=1):
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
        @param useFwdMsg: whether to use fwdmsg passing
        @type useFwdMsg: Int
        """
        self.data=None
        self.useFwdMsg=useFwdMsg
        self.name=name
        self.fwdMsgFName=name+"%04d"%readCountBit
        self.dims=dims
        self.dtype=dtype
        self.readCountBit=readCountBit
        self.nReaders=nReaders#number of readers of this shm.
        self.openShmFile()
        #self.shmLock=shmRWBLock(name,nReaders,readCountBit,0)
        self.semid=cmod.shmem.newsemid(name="/dev/shm"+name,existingFile=1,nSems=self.nReaders*2+2)
        print "remoteSHMGet: Waiting for initialisation of semid"
        cmod.shmem.semop(self.semid,nReaders*2+1,-1)#wait for the initialisation to take place (remoteSHMSend.py in initDataLock())
        self.fwdMsgArr=None
        self.fwdMsgLock=None
        if self.useFwdMsg:
            self.openFwdMsgArr()
            self.fwdMsgLock=cmod.shmem.newlock(name=self.fwdMsgFName)#simple shm lock (semid)
        print "remoteSHMGet: Initialisation complete"
    def __del__(self):
        """Don't think this deletion is needed"""
        print "Destructing remoteSHMGet - not deleting any semaphores"
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
                print "remoteSHMGet couldn't open SHM array %s.  Waiting 2s and retrying"%self.name
                self.data=None
                time.sleep(2)
        print "remoteSHMGet: Opened shm file %s okay"%self.name
    def freeReadLock(self):
        """Free the readlock on the shared memory data."""
        #cmod.shmem.semop(self.semid,range(self.nReaders*2),[-1]*self.nReaders+[1]*self.nReaders)
        cmod.shmem.semop(self.semid,[self.readCountBit,self.nReaders+self.readCountBit],[-1,1])
    def getReadLock(self):
        """Waits until new data has been written to the shared memory array, and until the writing of this data is complete."""
        #wait for all bitcnt==1.
        #cmod.shmem.semop(self.semid,range(self.nReaders,self.nReaders*2),[0]*self.nReaders)
        cmod.shmem.semop(self.semid,self.nReaders+self.readCountBit,0)
        
    def openFwdMsgArr(self):
        """Opens a fwdMsg shared memory array"""
        if self.fwdMsgArr!=None:
            cmod.shmem.unmap(self.fwdMsgArr)
            self.fwdMsgArr=None
        try:
            self.fwdMsgArr=cmod.shmem.open((28,),"1",self.fwdMsgFName,0)
        except:
            print "Couldn't open fwdMsgArr %s in remoteSHMGet.  Retrying in 2s."%self.fwdMsgFName
            time.sleep(2)
            self.fwdMsgArr=cmod.shmem.open((28,),"1",self.fwdMsgFName,0)
            
    def NoneFlagIsSet(self):
        """Checks whether there is data in the shared memory array, or whether
        the science successor (parent) returned None.
        @return: Flag representing 1 for None, or 0 for array data
        @rtype: Int"""
        return cmod.shmem.getSemValue(self.semid,self.nReaders*2)
