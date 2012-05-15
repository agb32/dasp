"""Simple reader-writer locks in Python
Many readers can hold the lock XOR one and only one writer.
Created elsewhere, slightly altered by agb.
"""
import threading,thread

#version = """$Id: rwlock.py,v 1.3 2005/11/17 13:40:56 ali Exp $"""

class RWLock:
  """
  A simple reader-writer lock Several readers can hold the lock
  simultaneously, XOR one writer. Write locks have priority over reads to
  prevent write starvation.

  Class variables (important to simulation programmer):
   - timeout - float, timeout.
  Class variables (important to simulation programmer):
   - rwlock - int, the current lock type.
   - writers_waiting - int, number of writers waiting.
   - monitor - threading.Lock object for synchronisation
   - readers_ok - threading.Condition object for synchronisation
   - writers_ok - threading.Condition object for synchronisation
  @cvar timeout: Timeout when obtaining a lock
  @type timeout: None or Float
  @cvar rwlock: Type of lock currently obtained
  @type rwlock: Int
  @cvar writers_waiting: Number of writers waiting to get lock
  @type writers_waiting: Int
  @cvar monitor: A lock object for synchronisation
  @type monitor: threading.Lock
  @cvar readers_ok: for synchronisation
  @type readers_ok: threading.Condition instance
  @cvar writers_ok: for synchronisation
  @type writers_ok: threading.Condition instance
  """
  def __init__(self,timeout=None):
    """Create a lock.  The timeouts dont really work.  Dont use them!
    @param timeout: A (not working) timeout
    @type timeout: Float
    """
    self.rwlock = 0
    self.writers_waiting = 0
    self.monitor = threading.Lock()
    self.readers_ok = threading.Condition(self.monitor)
    self.writers_ok = threading.Condition(self.monitor)
    self.timeout=timeout
  def acquire_read(self,timeout=None):
    """Acquire a read lock. Several threads can hold this typeof lock.
    It is exclusive with write locks.
    @param timeout: A (not working) timeout
    @type timeout: Float
    """
    #print thread.get_ident(),"rwlock: acquiring read"
    self.monitor.acquire()
    while self.rwlock < 0 or self.writers_waiting:
      if timeout==None:
        timeout=self.timeout
      self.readers_ok.wait(timeout)
    self.rwlock += 1
    self.monitor.release()
  def test_acquire_read(self):
    """Try to acquire read lock, but fail if cannot do so immediately.
    @return: Flag, whether have acquired the lock
    @rtype: Int
    """
    rtval=0
    self.monitor.acquire()
    if self.rwlock<0 or self.writers_waiting:
      #failed to acquire
      rtval= -1
    else:
      self.rwlock+=1
    self.monitor.release()
    #print thread.get_ident(),"rwlock: trying to aqcuire read (0==success)",rtval
    return rtval
  def test_acquire_write(self):
    """Try to acquire write lock, but fail if cannot do so immediately.
    @return: Flag, whether have acquired the lock
    @rtype: Int
    """
    rtval=0
    self.monitor.acquire()
    if self.rwlock!=0:
      #failed to acquire
      rtval=-1
    else:
      self.rwlock=-1
    self.monitor.release()
    #print thread.get_ident(),"rwlock: trying to aquire write (0==success)",rtval
    return rtval
  def acquire_write(self,timeout=None):
    """Acquire a write lock. Only one thread can hold this lock, and
    only when no read locks are also held.
    @param timeout: A (not working) timeout
    @type timeout: Float
    """
    #print thread.get_ident(),"rwlock: acquiring write"
    self.monitor.acquire()
    while self.rwlock != 0:
      self.writers_waiting += 1
      if timeout==None:
        timeout=self.timeout
      self.writers_ok.wait(timeout)
      self.writers_waiting -= 1
    self.rwlock = -1
    self.monitor.release()
  def promote(self,timeout=None):
    """Promote an already-acquired read lock to a write lock
    WARNING: it is very easy to deadlock with this method
    @param timeout: A (not working) timeout
    @type timeout: Float
    """
    self.monitor.acquire()
    self.rwlock -= 1
    while self.rwlock != 0:
      self.writers_waiting += 1
      if timeout==None:
        timeout=self.timeout
      self.writers_ok.wait(timeout)
      self.writers_waiting -= 1
    self.rwlock = -1
    self.monitor.release()
  def demote(self):
    """Demote an already-acquired write lock to a read lock
    """
    self.monitor.acquire()
    self.rwlock = 1
    self.readers_ok.notifyAll()
    self.monitor.release()
  def release(self):
    """Release a lock, whether read or write."""
    self.monitor.acquire()
    if self.rwlock < 0:#release the write lock.
      #print thread.get_ident(),"rwlock: releasing write lock"
      self.rwlock = 0
    elif self.rwlock==0:#called in error - doesn't hold lock
      print thread.get_ident(),"ERROR: CANNOT RELEASE RWLOCK THAT HASNT BEEN OBTAINED"
      raise "ERROR: cannot release rwlock"
    else:#release one read lock (a programming error may cause this to happen more than once).
      #print thread.get_ident(),"rwlock: releasing read lock"
      self.rwlock -= 1
    wake_writers = self.writers_waiting and self.rwlock == 0
    wake_readers = self.writers_waiting == 0
    self.monitor.release()
    if wake_writers:
      self.writers_ok.acquire()
      self.writers_ok.notify()
      self.writers_ok.release()
    elif wake_readers:
      self.readers_ok.acquire()
      self.readers_ok.notifyAll()
      self.readers_ok.release()

if __name__ == '__main__':
  import time
  rwl = RWLock()
  class Reader(threading.Thread):
    def run(self):
      print thread.get_ident(),self, 'start'
      rwl.acquire_read()
      print thread.get_ident(),self, 'acquired'
      time.sleep(5)    
      print self, 'stop'
      rwl.release()
  class Writer(threading.Thread):
    def run(self):
      print self, 'start'
      ac=-1
      while ac==-1:#try to acquire several times, returning if cant immediately
        
        ac=rwl.test_acquire_write()
        #rwl.acquire_write()
        print self, 'acquired?',ac
        time.sleep(10)    
      if ac==0:
        print self, 'stop'
        rwl.release()
  class ReaderWriter(threading.Thread):
    def run(self):
      print self, 'start'
      rwl.acquire_read()
      print self, 'acquired'
      time.sleep(5)    
      rwl.promote()
      print self, 'promoted'
      time.sleep(5)    
      print self, 'stop'
      rwl.release()
  class WriterReader(threading.Thread):
    def run(self):
      print self, 'start'
      rwl.acquire_write()
      print self, 'acquired'
      time.sleep(10)    
      print self, 'demoted'
      rwl.demote()
      time.sleep(10)    
      print self, 'stop'
      rwl.release()
  #Reader().start()
  #time.sleep(1)
  #Reader().start()
  #time.sleep(1)
  #ReaderWriter().start()
  #time.sleep(1)
  #WriterReader().start()
  #time.sleep(1)
  #Reader().start()
  Reader().start()
  time.sleep(1)
  Reader().start()
  time.sleep(1)
  Writer().start()
  time.sleep(1)
  Reader().start()
