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
