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
import threading
class MPIEngine:
    """DEPRECIATED Class to synchronise calls to MPI methods.  Since MPI is not thread
    safe, we need to call MPI calls in a known order, hence the use of this
    class.
    Class variables (important to simulation user):
     - ncomms - int, the number of different MPI access points required.
     - debug - debug message to be printed.
    Class variables (not important to simulation user):
     - curOrder - int, the current access point which can use MPI calls.
     - events - list, the list of threading.Event objects
    @cvar ncomms: Number of MPI access points
    @type ncomms: Int
    @cvar curOrder: Current MPI access point
    @type curOrder:
    @cvar events: List of event objects
    @type events: List
    @cvar debug: Debug string (None if no string)
    @type debug: None or user defined
    """
    
    def __init__(self,ncomms,debug=None):
        """Initialise the MPIEngine object.
        @param ncomms: Number of different MPI access points
        @type ncomms: Int
        @param debug: Debug string (None if no string)
        @type debug: None or user defined
        """
        self.ncomms=ncomms*2
        self.debug=debug
        self.curOrder=0
        self.iter=0
        self.lock=threading.Lock()
        self.iterevents={}
        self.events=[]
        for i in range(self.ncomms):
            self.events.append(threading.Event())
        self.events[0].set()
        
    def myturn(self,order,iter):
        """Blocks thread until its turn for MPI communication
        @param order:  The access point requesting access.
        @type order: Int
        """
        if self.iterevents.has_key(self.iter):
            self.iterevents[self.iter].set()
        if iter!=self.iter:
            if self.debug!=None:
                print "MPIEngine: order %g called for iteration %d, but current is %d (debug=%s)"%(order/2.0,iter,self.iter,str(self.debug))
            #wait until correct iteration.
            if not self.iterevents.has_key(iter):
                if self.debug!=None:
                    print "MPIEngine: Creating event for iteration %d (order=%g, debug=%s)"%(iter,order/2.0,str(self.debug))
                self.iterevents[iter]=threading.Event()
            self.iterevents[iter].wait()
        #now sure its the correct iteration.
            
        while self.curOrder!=order or self.iter!=iter:
            if self.debug!=None:
                print "MPIEngine: Order %g waiting, current iter %d, waiting for iter %d (debug=%s)"%(order/2.0,self.iter,iter,str(self.debug))
            self.events[order].wait()
        #now its our turn...
        if self.debug!=None:
            print "MPIEngine: Order %g/%g (iter %d/%d) acquiring lock (debug=%s)"%(order/2.0,self.curOrder/2.0,iter,self.iter,str(self.debug))
        self.lock.acquire()
        if self.debug!=None:
            print "MPIEngine: Order %g (iter %d) acquired lock, got turn (debug=%s)"%(order/2.0,self.iter,str(self.debug))
        self.events[order].clear()
    
    def endturn(self):
        """End of turn - change curOrder and set an event ready for the next
        access point to awaken.
        """
        if self.debug!=None:
            print "MPIEngine: Order %g (iter %d) released lock, turn ending, setting event (debug=%s)"%(self.curOrder/2.0,self.iter,str(self.debug))
        self.curOrder+=1
        if self.curOrder>=self.ncomms:
            self.curOrder=0
            self.iter+=1
            if self.iterevents.has_key(self.iter):#delete the key
                del(self.iterevents[self.iter])
        self.events[self.curOrder].set()
    
        self.lock.release()

if __name__=="__main__":
    import time
    nthreads=4
    m=MPIEngine(nthreads)
    def mainloop(order):
        i=0
        while i<10:
            m.myturn(order)
            if order==1:#ok to grab it more than once...
                m.myturn(order)
            print "It's my turn! (I am %d)"%order
            m.endturn()
            if order==1:
                time.sleep(1)
            i+=1
            
    for i in range(nthreads):
        t=threading.Thread(target=mainloop,args=(i,))
        t.start()
