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

