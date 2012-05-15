# Tests semaphore locking ala Splitter.py
#

import threading, shmem, thread

def threadMain(ident):
   print "%6d : Entering threadMain" % (ident)
   if shmem.acquireAndWait(semid) == 0:
      shmem.trigEvent(semid)
   else:
      print "%6d : Waited for semaphore lock"
   print "%6d : Leaving threadMain" % (ident)

##

print "main ** beginning"

semid=shmem.newsemid(nSems=2)
shmem.initSemaphore(semid, 0, 1)
shmem.initSemaphore(semid, 1, 1)
tlist=[]

for x in [55,66]:
   print x
   t=threading.Thread(target=threadMain,args=([x]))
   tlist.append(t)
   t.start()

# wait for all threads to finish
for t in tlist:
   t.join()

print "main ** ending"

