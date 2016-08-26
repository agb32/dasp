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
import cmod.shmem,thread
class rwShmLock:
    def __init__(self,name,create=0,realfileflag=0,timeout=None):
        """Creates a read/write lock object based on file name (which must
        exist, and typically will be a shmem object.  If realfileflag==0,
        name will start with a / and actually be present in /dev/shm/, i.e.
        /test represents file /dev/shm/test
        timeout should be float.
        Such a lock can have many processes reading it, but only one process
        writing it, iff there are no readers"""
        self.name=name
        self.create=create
        self.semid=cmod.shmem.newlock(name,realfileflag=realfileflag,setone=create)
        self.lock=0#+1 if got read lock, -1 if got write lock
        if timeout==None:
            self.timeout=None
        else:
            self.timeout=float(timeout)
    def __del__(self):
        print "Destructing rwShmLock"
        cmod.shmem.semdel(self.semid)
        
    def getReadLock(self,timeout=None):
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
            cmod.shmem.getreadlock(self.semid,timeout)#raises an error if fails
        except:
            print thread.get_ident(),"Failed to get shmreadlock"
            rtval=0
        if rtval==1:
            self.lock=1
        return rtval
    def getWriteLock(self,timeout=None):
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
            cmod.shmem.getlock(self.semid,timeout)#raises an error if fails.
        except:
            print thread.get_ident(),"Failed to get shmwritelock"
            rtval=0
        if rtval==1:#got lock successfully...
            self.lock=-1
        return rtval
    def freeLock(self):
        rtval=1
        if self.lock==1:#got a read lock
            print thread.get_ident(),"rwShmLock: Freeing shm read lock",self.name
            cmod.shmem.freereadlock(self.semid)
            self.lock=0
        elif self.lock==-1:#got a write lock
            print thread.get_ident(),"rwShmLock: freeing shm write lock",self.name
            cmod.shmem.freelock(self.semid)
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
        #promote from read to write
        print thread.get_ident(),"promote not yet implemented"
    def demote(self,timeout=None):
        #demote from write to read lock
        print thread.get_ident(),"demote not yet implemented"
