#$Id: remoteMPIGet.py,v 1.14 2006/05/03 08:25:27 ali Exp $
import Scientific.MPI,types
import base.fwdMsg
import Numeric,threading
print "remoteMPIGet is depreciated.  Please use mpiGet instead"
class remoteMPIGet:
    """
    DEPRECIATED
    A class for receiving MPI communication data.  This class should be
    connected to an remoteMPISend object.
    This should be called by an object whos true parent is in
    MPI space::
    
     parent (generating data for child)
       |
     remoteMPISend
       |   - mpi connection.
     remoteMPIGet
       |
     child (calls mpiget.next())
    
    Class variables (important to simulation programmer):
     - mpiParent - object, an MPIParent object, giving details about the MPI parent that messages are expected from (rank, data tags etc).
     - useFwdMsg - int, whether the use of fwdMsg passing is allowed

    Class variables (not important to simulation programmer):
     - mpiComm - MPI communicator (Scientific.MPI.world.duplicate())
    @cvar mpiParentObj: MPIParent object holding information about the MPI
    connection from which data is expected.
    @type mpiParentObj: MPIParent object
    @cvar mpiComm: MPI communicator
    @type mpiComm: Scientific.MPI.world object
    @cvar mpiEngine: An MPIEngine object for synchronising MPI calls
    @type mpiEngine: Instance
    @cvar mpiOrder: The order in which this object should do MPI calls
    @type mpiOrder: Int
    @cvar mpiOrder: The order in which this object should do MPI synchronisation calls
    @type mpiOrder: Int
    @cvar useFwdMsg: Flag determining whether forward message passing should be
    expected and used
    @type useFwdMsg: Int
    @cvar debug: Flag, whether to print debug message (if None, won't print)
    @type debug: None or user defined.
    """
    def __init__(self,mpiParentObj,mpiComm,mpiEngine,mpiOrder,useFwdMsg=0,args={},debug=None):
        """Initialise the MPI connection details.
        @param mpiParentObj: MPIParent object with details about connection
        @type mpiParentObj: MPIParent object
        @param mpiComm: Duplicate of Scientific.MPI.world
        @type mpiComm: Object
        @param mpiEngine: An MPIEngine object for synchronising MPI calls
        @type mpiEngine: Instance
        @param mpiOrder: The order in which this object should do MPI calls
        @type mpiOrder: Int
        @param mpiOrder: The order in which this object should do MPI synchronisation calls
        @type mpiOrder: Int
        @param useFwdMsg: Flag, whether to use fwdMsg passing
        @type useFwdMsg: Int
        @param debug: Flag, whether to print debug message (if None, won't print)
        @type debug: None or user defined.
        """
        self.objID=self.__module__.split(".")[-1]
        if args.has_key("idstr"):
            self.objID=self.objID+"_"+args["idstr"]
        self.iter=0
        self.lock=threading.Lock()
        self.mpiParent=mpiParentObj
        self.mpiComm=mpiComm#Scientific.MPI.world.duplicate()
        self.useFwdMsg=useFwdMsg
        self.mpiEngine=mpiEngine
        self.realmpiOrder=mpiOrder
        if self.useFwdMsg:
            self.mpiOrder=2*self.realmpiOrder+1
            self.mpiSyncOrder=2*self.realmpiOrder
        else:
            self.mpiOrder=2*self.realmpiOrder
            self.mpiSyncOrder=2*self.realmpiOrder+1

        self.tmpmsg=base.fwdMsg.fwdMsg()
        self.debug=debug
    def next(self,childName,msg=None):
        """called by the child when it wants more data.
        Get the data from the MPI connection and return it.
        msg is an instance of fwdMsg.fwdMsg().

        @param msg: The message to pass to the parent (predecessor) object
        over the MPI connection (or None).
        @type msg: None or fwdMsg object
        @return: The data from the parent (predecessor) object
        @rtype: None or Numeric array
        """
        self.lock.acquire()#must happen atomically
        thisiter=self.iter
        self.iter+=1
        self.lock.release()
        if not self.useFwdMsg:
            if type(msg)!=types.NoneType:
                print "Error - fwdmsg received, when MPI not using FWDMSG"
                raise Exception("Error - fwdmsg received but mpi connection not using fwdmsg")
        else:#using fwdmsg.
            if type(msg)==types.NoneType:
                msg=base.fwdMsg.fwdMsg()#empty message...
            msg.toArray()
            self.mpiEngine.myturn(self.mpiSyncOrder,thisiter)#ok to grab this more than once - won't block if currently is our turn...
            self.mpiComm.send(msg.arr,self.mpiParent.sourceRank,self.mpiParent.tag)
            self.mpiEngine.endturn()
        self.mpiEngine.myturn(self.mpiOrder,thisiter)
        if self.debug!=None: print "remoteMPIGet: Receiving MPI message from rank %d, tag %d (debug=%s)"%(self.mpiParent.sourceRank,self.mpiParent.tag,str(self.debug))
        data,rank,tag,nelements=self.mpiComm.receive(self.mpiParent.data,self.mpiParent.sourceRank,self.mpiParent.tag)
        self.mpiEngine.endturn()
        if self.debug!=None: print "remoteMPIGet: Received MPI data (debug=%s)"%str(self.debug)
        if not self.useFwdMsg:#synchronise MPI calls...
            self.mpiEngine.myturn(self.mpiSyncOrder,thisiter)
            self.mpiComm.send(self.tmpmsg.arr,self.mpiParent.sourceRank,self.mpiParent.tag)
            self.mpiEngine.endturn()
        if nelements==0:#nelements is the number of elements of the array (not bytes) returned.
            data=None
        return data

class MPIParent:
    """An object which creates the array to hold the data
    received by the mpi connection and also specifies the MPI process from
    which data is expected.
    @cvar sourceRank: The MPI rank of the source object
    @type sourceRank: Int
    @cvar tag: The tag expected, associated with the data.
    @type tag: Int
    @cvar data: The data array for holding the MPI data received.
    @type data: Numeric array
    """
    def __init__(self,dims,dtype,sourceRank,tag):
        """Initialise the MPIParent object.
        @param dims: The dimensions of the Numeric array expected from the parent (predecessor) object
        @type dims: Tuple of Int
        @param dtype: The data type of the Numeric array
        @type dtype: Char
        @param sourceRank: The MPI rank of the source process
        @type sourceRank: Int
        @param tag: The tag associated with the MPI data
        @type tag: Int
        """
        self.sourceRank=sourceRank
        self.tag=tag
        if dims==None:#create the array every time (more memory/cpu intensive)
            self.data=dtype
        else:#always pass into this array.
            self.data=Numeric.zeros(dims,dtype)
