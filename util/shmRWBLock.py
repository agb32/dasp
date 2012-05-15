#$Id: shmRWBLock.py,v 1.5 2005/11/17 13:40:56 ali Exp $
import cmod.shmem,thread,sys
class shmRWBLock:
    """Class to create a read/write lock object with checking to ensure the
    writer can only get the lock once all readers have obtained and released it
    since the last write.

    Class variables (important to simulation programmer):
     - name - string, the name to create semaphore set on
     - nReaders - int, number of readers
     - myBit - int or None, the bit corresponding to this reader
     - create - int, flag whether to create the semaphore set
     - timeout - float or None, the timeout when acquiring a lock
    Class variables (not important to simulation programmer):
     - semid - int, semaphore set identifier
     - lock - int, specifying type of lock currently held
    
    @cvar name: the name to create semaphore set on
    @type name: String
    @cvar nReaders: number of readers
    @type nReaders: Int
    @cvar myBit: the bit corresponding to this reader
    @type myBit: None or Int
    @cvar create: flag whether to create the semaphore set
    @type create: Int
    @cvar timeout: the timeout when acquiring a lock
    @type timeout: Float or None
    @cvar semid: semaphore set identifier
    @type semid: Int
    @cvar lock: specifying type of lock currently held
    @type lock: Int
    """
    def __init__(self,name,nReaders,myBit,create=0,realfileflag=0,timeout=None):
        """Creates a read/write lock object based on file name (which must
        exist, and typically will be a shmem object.  If realfileflag==0,
        name will start with a / and actually be present in /dev/shm/, i.e.
        /test represents file /dev/shm/test
        timeout should be float.
        Such a lock can have many processes reading it, but only one process
        writing it, iff there are no readers
        This additionally implements a readCount bit set, which means that a
        bit (semaphore) is set for each reader (known in advance) when the
        writelock is obtained and then this bit for this reader is unset when
        the reader obtains the readlock.
        The writelock can only be obtained when all bits are zero, ie all
        the readers have read the data...
        myBit can be None if this is called by the process that writes but
        doesn't read the shm (probably when create==1).
        @param name: the name to create semaphore set on
        @type name: String
        @param nReaders: number of readers
        @type nReaders: Int
        @param myBit: the bit corresponding to this reader
        @type myBit: None or Int
        @param create: flag whether to create the semaphore set
        @type create: Int
        @param timeout: the timeout when acquiring a lock
        @type timeout: Float or None
        @param realfileflag: Whether a real file, or shared memory file
        @type realfileflag: Int
        """
        self.name=name
        self.nReaders=nReaders#number of readers expected (must be known)
        self.myBit=myBit#the bit to set for this process.
        if myBit>=nReaders:
            print "shmRWBlock error -",myBit,nReaders
            raise "MYBIT/NREADERS SHMRWBLOCK ERROR"
        self.create=create
        self.semid=shmem.newlock(name,realfileflag=realfileflag,setone=create,nReaders=nReaders)
        self.lock=0#+1 if got read lock, -1 if got write lock
        if timeout==None:
            self.timeout=None
        else:
            self.timeout=float(timeout)
    def __del__(self):
        """Destructor"""
        print "Destructing shmRWBlock"
        cmod.shmem.semdel(self.semid)
        
    def getReadLock(self,timeout=None):
        """Obtain a readlock.
        @param timeout: the timeout when acquiring a lock
        @type timeout: Float or None
        @return: 1 if successful, zero otherwise
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
            shmem.getreadlock(self.semid,timeout,self.myBit)#raises an error if fails or blocks if already have read the data.
        except:
            print thread.get_ident(),"Failed to get shmreadlock"
            rtval=0
        if rtval==1:
            self.lock=1
            print "readlock obtained"
        return rtval
    def getWriteLock(self,timeout=None):
        """Obtain the writelock.
        @param timeout: the timeout when acquiring a lock
        @type timeout: Float or None
        @return: 1 if successful, zero otherwise
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
            shmem.getlock(self.semid,timeout,self.nReaders)#raises an error if fails.
        except:
            print thread.get_ident(),"Failed to get shmwritelock"
            rtval=0
        if rtval==1:#got lock successfully...
            self.lock=-1
            print "write lock obtained"
        return rtval
    def freeLock(self):
        """Free the lock
        @return: 1 if successful, 0 otherwise
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
        """Not implemented"""
        #promote from read to write
        print thread.get_ident(),"promote not yet implemented"
    def demote(self,timeout=None):
        """Not implemented"""
        #demote from write to read lock
        print thread.get_ident(),"demote not yet implemented"

if __name__=="__main__":
    print "Testing shmRWBLock"
    if len(sys.argv)==1:
        create=1
        myBit=None
    else:
        myBit=int(sys.argv[1])
        create=0
    s=shmRWBLock("/home/ali/tmp2.ps",2,myBit,create,1)
    if not create:
        s.getReadLock()
        s.freeLock()
        s.getReadLock()
        s.freeLock()
        s.getReadLock()
        s.freeLock()
    else:
        s.getWriteLock()
        s.freeLock()
        s.getWriteLock()
        s.freeLock()
        s.getWriteLock()
        s.freeLock()
                        
