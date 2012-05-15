#test threading on shm locks...
import rwShmLock
import time,thread,threading,sys
if len(sys.argv)>1:
    cr=0
else:
    cr=1
shmlock=rwShmLock.rwShmLock("/home/ali/tmp.ps",cr,1)


def rd():
    print thread.get_ident(),"Getting readlock"
    shmlock.getReadLock()
    for i in range(10):
        print thread.get_ident(),"Readlock sleeping"
        time.sleep(1)

def wr():
    print thread.get_ident(),"Getting writelock"
    shmlock.getWriteLock()
    for i in range(10):
        print thread.get_ident(),"Writelock sleeping"
        time.sleep(1)

def wt():
    for i in range(10):
        print thread.get_ident(),"Other thread waiting..."
        time.sleep(1)
if len(sys.argv)>1:
    thread.start_new_thread(wt,())
    wr()
else:
    thread.start_new_thread(wt,())
    rd()
shmlock.freeLock()
