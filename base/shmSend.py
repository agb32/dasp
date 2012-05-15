#$Id: shmSend.py,v 1.13 2009/06/24 09:45:33 ali Exp $
import types,os,numpy
import cmod.shmem
import base.aobase
#also see shmSend2 below... I think shmSend is probably depreciated - I think they do the same thing anyway...
debugConnections=0
class shmSend:
    """
    Class to implement the sending end of shared memory communications.
    If you plan to alter this, then it needs some more thought when it
    comes to sending shm data.

    It may be most efficient to create this object, and then pass the shared memory array
    to the parent object to use for its outputData.  This would avoid copying the data twice.
    However, if you do this, you should then use the makeCopy=1 flag in shmGet, as otherwise
    you could be reading and writing the data at the same time.  So you have to copy somewhere,
    but where is up to you.

    Class variables (important to simulation programmer):
     - parent - object, the predecessor (parent) science object
     - shmInfo - SHMInfo instance, providing information about the shared memory connection
    @cvar parent: the predecessor (parent) science object
    @type parent: object
    @cvar shmInfo: information about shared memory connection
    @type shmInfo: SHMInfo instance
    @cvar objID: Object identifier
    @type objID: String
    @cvar debug: Message printed if debug messages are to be printed
    @type debug: Anything
    """
    def __init__(self,parent,shmInfo,args={},debug=None):
        """Initialise the SHM connection object.
        @param parent: predecessor (parent) object for retrieving science data
        @type parent: Instance
        @param shmInfo: shared memory information
        @type shmInfo: SHMInfo
        @param args: Additional arguments, currently supporting idstr.
        @type args: Dict
        @param debug: Message to send with debug messages
        @type debug: None or anything.
        """
        #Lock object for fwdMsg can be a standard lock (not rwlock).
        #Lock for data needs to be an array of bit counts.
        print "shmSend.shmSend is depreciated, please use shmSend.shmSend2, or newShmSend()"
        self.objID=self.__module__.split(".")[-1]
        if args.has_key("idstr"):
            self.objID=self.objID+"_"+args["idstr"]
        self.debug=debug
        self.parent=parent
        self.shmInfo=shmInfo
        #self.firstIter=1#this is because we start owning the write lock...

    def generateNext(self,msg=None):
        """Obtain data from the parent object, and place in the shared memory
        array.\n
        First wait for a sync call, and then until all children have finished reading
        the data and get the write lock.  Note, the write lock can only be obtained once all
        children have requested and freed the read lock.
        @param msg: None
        @type msg: None
        """
        if self.debug!=None: print "shmSend: in GenerateNext (debug=%s)"%str(self.debug)
        if self.parent.dataValid==1:#store the data
            #self.outputDataList.append(self.parent.outputData.copy())
            self.outputData=self.parent.outputData
            self.dataValid=1
        else:
            if self.debug!=None: print "shmSend: parent data not valid (debug=%s)"%str(self.debug)
            self.dataValid=0
            #self.outputDataList.append("")

        if self.debug!=None:
            print "shmSend: blocking on sem A (debug=%s)"%str(self.debug)
        self.shmInfo.waitForSync()
        if self.debug!=None: print "shmSend: Getting write lock (debug=%s)"%str(self.debug)
        self.shmInfo.getWriteLock()
        if self.debug!=None: print "shmSend: Writelock obtained (debug=%s)"%str(self.debug)
        #self.mpiComm.receive(self.syncmsg,self.mpiChild.rank,self.mpiChild.tag)
        nGen=cmod.shmem.getSemValue(self.shmInfo.semid,self.shmInfo.nChild*2+2)
        if nGen==self.shmInfo.nChild:#all generate=1
            self.parent.setGenerate(1)
            #data=self.outputDataList.pop(0)
            self.copyToSHM(self.outputData)
            self.shmInfo.freeNoneFlag()
        elif nGen==0:#all generate=0
            self.parent.setGenerate(0)
            self.shmInfo.setNoneFlag()
        else:
            print "shmSend: ERROR - readers differ in choice of generate or not.  %d should equal %d or 0 (debug=%s)"%(nGen,self.shmInfo.nChild,str(self.debug))
            raise Exception("shmSend error - readers differ in choice of generate")
        self.shmInfo.resetGenerateInfo()
        if self.debug!=None:
            print "shmSend: Freeing writelock (debug=%s)"%str(self.debug)
        self.shmInfo.freeWriteLock()



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
            self.shmInfo.data[:,]=data#copy data to shm region.
        else:#data is already in the shared memory array.
            pass
    
        
class shmInfo:
    """A class holding information about the shared memory arrays, and used
    for creating and removing them.

    Class variables (important to simulation programmer):
     - name - string, the identifier for the data SHM array, should commence with /
     - dims - tuple, the dimensions of the data SHM array
     - dtype - char, the datatype of the data SHM array
     - nChild - int, the number of children reading the shared memory array

    Class variables (not important to simulation programmer):
     - semid - int, the semaphore identifier, used for process synchronisation

    @cvar name: Identifier for the data SHM array, should commence with /
    @type name: String
    @cvar dims: Dimensions of the data SHM array
    @type dims: Tuple of Int
    @cvar dtype: Datatype for data SHM array
    @type dtype: Char
    @cvar nChild: Number of children to read the SHM array
    @type nChild: Int
    @cvar semid: Semaphore identifier used for process synchronisation
    @type semid: Int
    """
    def __init__(self,name,dims,dtype,nChild):
        """Initialise the SHMInfo.
        @param name: Identifier for the data SHM array, should commence with /
        @type name: String
        @param dims: Dimensions of the data SHM array
        @type dims: Tuple of Int
        @param dtype: Datatype for data SHM array
        @type dtype: Char
        @param nChild: Number of children to read the SHM array
        @type nChild: Int
        """
        self.data=None
        self.name=name+"_"+os.environ["USER"]
        name=self.name
        self.dims=dims
        self.dtype=dtype
        self.nChild=nChild
        self.semid=None
        self.openShmFile()
        self.initDataLock()
    def __del__(self):
        print "Destructing shmInfo object",self.semid
        #cmod.shmem.semdel(self.semid)
        #self.closeShmFile()
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
        self.data[:,]=0#initialise memory (prob not needed)
    def initDataLock(self):
        """Initialise the read/write lock for the SHM array.  The initial state
        is such that the write lock can be obtained, but not the read locks.
        Should be called only once.
        """
        self.semid=cmod.shmem.newsemid(name="/dev/shm"+self.name,existingFile=1,nSems=self.nChild*2+4)#+1 is for the none type flag and +4 is for the initialised flag (ie so shmGet instances don't continue until this flag has been set), +2 is for the countdown flag (block writer) and +3 is for the generate flag.
        print "shmSend: Initialised semid for data lock (semid=%d)"%self.semid
        for i in range(self.nChild):#set all bitCnts to 0.
            cmod.shmem.initSemaphore(self.semid,i,0)
        for i in range(self.nChild,self.nChild*2):#set all bits to zero
            cmod.shmem.initSemaphore(self.semid,i,1)
        self.freeNoneFlag()
        cmod.shmem.initSemaphore(self.semid,self.nChild*2+1,self.nChild)
        cmod.shmem.initSemaphore(self.semid,self.nChild*2+2,0)
        cmod.shmem.initSemaphore(self.semid,self.nChild*2+3,self.nChild)
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
    def waitForSync(self):
        """Wait for reader synchronisation (ready to read)"""
        cmod.shmem.semop(self.semid,self.nChild*2+1,0)
    def resetGenerateInfo(self):
        """Reset reader synchronisation for next time"""
        cmod.shmem.initSemaphore(self.semid,self.nChild*2+1,self.nChild)
        cmod.shmem.initSemaphore(self.semid,self.nChild*2+2,0)

class shmSend2(base.aobase.aobase):
    def __init__(self,parent,name,dims,dtype,args={},debug=None,idstr=None):
        if debugConnections==1:
            if debug==None:
                debug="shmSend %d"%name
        base.aobase.aobase.__init__(self,parent,None,args=args,debug=debug,idstr=idstr)
        self.name=name+"_"+os.environ["USER"]
        self.data=None
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
        self.semid=None
        self.openShmFile()
        self.initDataLock()

    def __del__(self):
        #print "Destructing shmSend object",self.semid
        try:
            cmod.shmem.semdel(self.semid)
        except:#This seems to raise an error, but still deletes the sem.
            #print "Error during shmSend.__del__"
            pass
        #self.closeShmFile()

    def closeShmFile(self):
        """Close the data SHM array"""
        cmod.shmem.unmap(self.data)
        self.data=None
        cmod.shmem.unlink(self.name)
    def openShmFile(self):
        """Open the data SHM array (closing first if already open)"""
        if self.data!=None:
            self.closeShmFile()
        #print "Opening shmFile %s with %s %s %s"%(self.name,str(type(self.name)),str(type(self.dims)),str(type(self.dtype)))
        self.data=cmod.shmem.open(self.dims,self.dtype,self.name,1)
        self.data[:,]=0#initialise memory (prob not needed)
    def initDataLock(self):
        """Initialise the semaphores for the SHM array.  The initial state
        is such that the write lock cant be obtained.
        Should be called only once.
        """
        self.semid=cmod.shmem.newsemid(name="/dev/shm"+self.name,existingFile=1,nSems=3)
        print "shmSend: Initialised semid for data lock (semid=%d)"%self.semid
        cmod.shmem.initSemaphore(self.semid,0,0)
        cmod.shmem.initSemaphore(self.semid,1,0)
        cmod.shmem.initSemaphore(self.semid,2,1)#tell the shmGet we're ready...

    def setNoneFlag(self,val):
        cmod.shmem.initSemaphore(self.semid,2,val)

    def copyToSHM(self,data):
        """Copy data into the shared memory, ready for passing to clients.
        @param data: The data to be shared
        @type data: None or Numeric array
        """
        if data is not self.data:#may already be the same objects...
            if type(self.data)==types.NoneType:
                self.dims=data.shape
                self.dtype=data.typecode()
                print "Possible error - opening SHM file"
                self.openShmFile()#open the shm file
            self.data[:,]=data#copy data to shm region.
        else:#data is already in the shared memory array.
            pass

    def generateNext(self,msg=None):
        if self.parent.dataValid==1:#store the data
            self.outputData=self.parent.outputData
            self.dataValid=1
        else:
            if self.debug!=None: print "shmSend: parent data not valid (debug=%s)"%str(self.debug)
            self.dataValid=0
            #self.outputDataList.append("")

        if self.debug!=None: print "shmSend: blocking on sem 0 (debug=%s)"%str(self.debug)

        cmod.shmem.semop(self.semid,0,-1)#block until can write
        self.generate=cmod.shmem.getSemValue(self.semid,2)
        if self.generate:
            self.parent.setGenerate(1)
            self.copyToSHM(self.outputData)
            self.setNoneFlag(0)
        else:
            self.parent.setGenerate(0)
            self.setNoneFlag(1)
        cmod.shmem.semop(self.semid,1,1)#unblock the reader...


def newSHMSend(parent,name,dims,dtype,idstr=None):
    ss=shmSend2(parent,name,dims,dtype,idstr=idstr)
    return ss
