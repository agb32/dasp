#$Id: remoteMPISend.py,v 1.13 2006/05/03 08:25:27 ali Exp $
import Scientific.MPI,threading
import base.fwdMsg
from types import NoneType
print "remoteMPISend is depreciated.  Please use mpiSend instead"
class mpiSend:
    """
    DEPRECIATED
    Send data over an MPI connection to the child
    mpiChild is a MPIChild object::
    
    P
    |
    o  - remoteMPISend instance.
    |  - mpi connection.
    o  - remoteMPIGet instance.
    |
    C

    Note that MPI communications within a process may (will) lead to
    errors.  So, for example, don't let a single process send data via
    MPI to itself (it would be stupid to do this anyway, due to
    efficiency considerations).
    
    Class variables (important to simulation programmer):
     - parent - object, the simulation object generating data to be passed over the connection.
     - mpiChild - an MPIChild object specifying the MPI rank of the receiving process and the tag to be sent with the data
     - mpiComm - a Scientific.MPI.world object
     - mpiEngine - a MPIEngine object used for synchronising when this instance can make its MPI calls.  MPI isn't thread safe, so all calls must be made one after the other...
     - mpiOrder - the place in the queue where this mpi call will be carried out.
     - useFwdMsg - int flag, whether to use fwdMsg passing or not.
    Class variables (not important to simulation programmer):
     - msg - None or a FwdMsg instance
    @cvar parent: Predecessor (parent) object from which data is obtained.
    @type parent: Science object
    @cvar mpiChild: Structure detailing MPI connection.
    @type mpiChild: MPIChild object
    @cvar useFwdMsg: Flag, whether fwdMsg passing is to be used
    @type useFwdMsg: Int
    @cvar mpiComm:  MPI communicator object
    @type mpiComm: Scientific.MPI.world
    @cvar mpiEngine: An MPIEngine object for synchronising MPI calls
    @type mpiEngine: Instance
    @cvar mpiOrder: The order in which this object should do MPI calls
    @type mpiOrder: Int
    @cvar debug: Flag, whether to print debug message (if None, won't print)
    @type debug: None or user defined.
    """
    def __init__(self,parent,mpiChildObj,mpiComm,mpiEngine,mpiOrder,useFwdMsg=0,args={},debug=None):
        """Initialise the MPI connection.
        @param parent: Science object from which data is obtained
        @type  parent: Object
        @param mpiChildObj: Details of MPI connection
        @type  mpiChildObj: MPIChild
        @param mpiComm: Duplicate of Scientific.MPI.world
        @type mpiComm: Object
        @param mpiEngine: An MPIEngine object for synchronising MPI calls
        @type mpiEngine: Instance
        @param mpiOrder: The order in which this object should do MPI calls
        @type mpiOrder: Int
        @param mpiOrder: The order in which this object should do MPI synchronisation calls
        @type mpiOrder: Int
        @param useFwdMsg: Flag, whether fwdMsg passing is used
        @type  useFwdMsg: Int
        @param debug: Flag, whether to print debug message (if None, won't print)
        @type debug: None or user defined.
        """
        print "remoteMPISend: Initialising"
        self.objID=self.__module__.split(".")[-1]
        if args.has_key("idstr"):
            self.objID=self.objID+"_"+args["idstr"]
        self.iter=0
        self.lock=threading.Lock()
        self.parent=parent
        self.mpiChild=mpiChildObj
        self.useFwdMsg=useFwdMsg
        self.mpiEngine=mpiEngine
        self.realmpiOrder=mpiOrder
        if self.useFwdMsg:
            self.mpiOrder=2*self.realmpiOrder+1
            self.mpiSyncOrder=2*self.realmpiOrder
        else:
            self.mpiOrder=2*self.realmpiOrder
            self.mpiSyncOrder=2*self.realmpiOrder+1
        #print "remoteMPISend: creating mpi.world.duplicate()"
        self.mpiComm=mpiComm#Scientific.MPI.world.duplicate()
        #print "remoteMPISend: done MPI.world.duplicate() (rank %d)"%self.mpiComm.rank
        if useFwdMsg:
            self.msg=base.fwdMsg.fwdMsg()
        else:
            self.msg=None
            self.tmpmsg=base.fwdMsg.fwdMsg()
        self.debug=debug
        print "remoteMPISend: Initialised"
    def next(self,childName,msg=None):
        """Obtain the next data from the parent object, and pass over MPI
        connection to successor object.
        @param msg: fwdMsg or None
        @type msg: Instance or None
        @return: The data
        @rtype: None or Numeric array
        """
        self.lock.acquire()#must happen atomically
        thisiter=self.iter
        self.iter+=1
        self.lock.release()
        if self.useFwdMsg:#wait for the message to arrive.
            self.mpiEngine.myturn(self.mpiSyncOrder,thisiter)
            self.mpiComm.receive(self.msg.arr,self.mpiChild.rank,self.mpiChild.tag)
            self.mpiEngine.endturn()
            self.msg.fromArray()#populate msg.i, msg.f, msg.strng
        if self.debug!=None: print "remoteMPISend: Calling parent.next() (debug=%s)"%str(self.debug)
        data=self.parent.next(self.objID,self.msg)
        if type(data)==NoneType:
            tdata=""
        else:
            tdata=data
        self.mpiEngine.myturn(self.mpiOrder,thisiter)#wait til its our turn!
        if self.debug!=None: print "remoteMPISend: Sending data (debug=%s)"%str(self.debug)
        self.mpiComm.send(tdata,self.mpiChild.rank,self.mpiChild.tag)#blocking
        if self.debug!=None: print "remoteMPISend: Data sent (debug=%s)"%str(self.debug)
        self.mpiEngine.endturn()
        if not self.useFwdMsg:#synchronise MPI calls...
            self.mpiEngine.myturn(self.mpiSyncOrder,thisiter)
            self.mpiComm.receive(self.tmpmsg.arr,self.mpiChild.rank,self.mpiChild.tag)
            self.mpiEngine.endturn()
        return data

class MPIChild:
    """Structure class holding information about connecting object via MPI.
    @cvar rank: MPI rank of the connecting object (to which data is sent).
    @type rank: Int
    @cvar tag: Tag sent with data to connecting object
    @type tag: Int
    """
    def __init__(self,rank,tag):
        """Rank is the rank of the child process.  tag is a tag used
        to identify which signals we want"""
        self.rank=rank
        self.tag=tag
