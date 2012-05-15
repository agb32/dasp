import cmod.shmem,thread
class rwShmLock:
    """Class for implementing a read/write lock on shared memory.  This allows
    many readers exclusive with one writer.

    Class variables (important to simulation programmer):
     - name - string, shared memory identifier
     - create - int, whether to create the semaphore set
     - timeout - None or Float, whether to use a timeout while waiting for lock
    Class variables (important to simulation programmer):
     - semid - int, the semaphore ID
     - lock - int, whether has the lock
    @cvar name: Shared memory identifier
    @type name: String
    @cvar create: Whether to create the semaphore set
    @type create: Int
    @cvar semid: Semaphore identifier, used for locking
    @type semid: Int
    @cvar lock: Whether has lock
    @type lock: Int
    @cvar timeout: Time waiting for lock
    @type timeout: None or Float
    """
    def __init__(self,name,create=0,realfileflag=0,timeout=None):
        """Creates a read/write lock object based on file name (which must
        exist, and typically will be a shmem object.  If realfileflag==0,
        name will start with a / and actually be present in /dev/shm/, i.e.
        /test represents file /dev/shm/test
        timeout should be float (or int).
        Such a lock can have many processes reading it, but only one process
        writing it, iff there are no readers
        @param name: The shared memory identifier (in which case starts with /) or other filename
        @type name: String
        @param create: Create the semaphore set or not
        @type create: Int
        @param realfileflag: Whether using a shared memory object or a real file
        @type realfileflag: Int
        @param timeout: Timeout for obtaining the lock
        @type timeout: None or Float
        
        """
        self.name=name
        self.create=create
        self.semid=shmem.newlock(name,realfileflag=realfileflag,setone=create)
        self.lock=0#+1 if got read lock, -1 if got write lock
        if timeout==None:
            self.timeout=None
        else:
            self.timeout=float(timeout)
    def __del__(self):
        """Destructor - cleanup the semaphore set"""
        print "Destructing rwShmLock"
        cmod.shmem.semdel(self.semid)
    def getReadLock(self,timeout=None):
        """Obtain a readlock, optionally waiting only for a set time before failing if not obtained within this time.
        @param timeout: The timeout to wait for
        @type timeout: None or Float
        @return: Error indicator flag (1 means no error, lock obtained)
        @rtype: Int
        """
        rtval=1
        print thread.get_ident(),"rwShmLock: getting read lock",self.name
        if self.lock!=0:
            print thread.get_ident(),"Already got lock (+1=read, -1=write):",self.lock
            raise "ALREADY GOT LOCK"
        if timeout==None:
            if self.timeout==None:
                timeout=-1.0#a negative timeout means no timeout
            else:
                timeout=self.timeout
        else:
            timeout=float(timeout)
        print thread.get_ident(),"timeout=",timeout
        try:
            shmem.getreadlock(self.semid,timeout)#raises an error if fails
        except:
            print thread.get_ident(),"Failed to get shmreadlock"
            rtval=0
        if rtval==1:
            self.lock=1
        return rtval
    def getWriteLock(self,timeout=None):
        """Obtain a writelock, optionally waiting only for a set time before failing if not obtained within this time.
        @param timeout: The timeout to wait for
        @type timeout: None or Float
        @return: Error indicator flag (1 means no error, lock obtained)
        @rtype: Int
        """
        rtval=1
        print thread.get_ident(),"rwShmLock: getting write lock",self.name
        if self.lock!=0:
            print thread.get_ident(),"Already got lock (+1=read, -1=write):",self.lock
            raise "ALREADY GOT LOCK"#this is an error, because we really don't want this lock to be reenterant...
        if timeout==None:
            if self.timeout==None:
                timeout=-1.0#-ve timeout means no timeout
            else:
                timeout=self.timeout
        else:
            timeout=float(timeout)
        try:
            shmem.getlock(self.semid,timeout)#raises an error if fails.
        except:
            print thread.get_ident(),"Failed to get shmwritelock"
            rtval=0
        if rtval==1:#got lock successfully...
            self.lock=-1
        return rtval
    def freeLock(self):
        """Free the lock (whether read or write).  Allows other processes to obtain it.
        @return: Return error indicator (1 if freed okay).
        @rtype: Int
        """
        rtval=1
        if self.lock==1:#got a read lock
            print thread.get_ident(),"rwShmLock: Freeing shm read lock",self.name
            shmem.freereadlock(self.semid)
            self.lock=0
        elif self.lock==-1:#got a write lock
            print thread.get_ident(),"rwShmLock: freeing shm write lock",self.name
            shmem.freelock(self.semid)
            self.lock=0
        else:
            print thread.get_ident(),"Cannot free a shmlock we don't own (probably non-critical)"
            #this may happen if the lock is created while we're in a simloop, ie if a new shm connection is created because a client connects...
            rtval=0
            #dont raise an exception - may not be an error - ie if the lock was
            #introduced while we were elsewhere doing a calc which then tried to free the lock!
            #raise "CANNOT FREE LOCK: NO LOCK TO FREE"
        return rtval
    def promote(self,timeout=None):
        """Not yet implemented.
        Promote a readlock to a write lock"""
        #promote from read to write
        print thread.get_ident(),"promote not yet implemented"
    def demote(self,timeout=None):
        """Not yet implemented.
        Demote from write to read lock
        """
        #demote from write to read lock
        print thread.get_ident(),"demote not yet implemented"
