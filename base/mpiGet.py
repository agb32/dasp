#$Id: mpiGet.py,v 1.11 2009/02/18 10:28:41 ali Exp $
#import Scientific.MPI
import numpy,time
import base.aobase

debugConnections=0
class mpiGet(base.aobase.aobase):
    """
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
    @cvar objID: Object identifier
    @type objID: String
    @cvar generate: Whether to generate this iteration
    @type generate: Int
    @cvar dataValid: Whether the data is valid
    @type dataValid: Int
    @cvar syncarr: Array for synchronisation
    @type syncarr: Array
    @cvar debug: Flag, whether to print debug message (if None, won't print)
    @type debug: None or user defined.
    """
    def __init__(self,mpiParentObj,mpiComm,args={},debug=None,idstr=None):
        """Initialise the MPI connection details.
        @param mpiParentObj: MPIParent object with details about connection
        @type mpiParentObj: MPIParent object
        @param mpiComm: Duplicate of Scientific.MPI.world
        @type mpiComm: Object
        @param args: Dictionary of arguments, at present can contain keys: idstr
        @type args: Dict
        @param debug: Flag, whether to print debug message (if None, won't print)
        @type debug: None or user defined.
        """
        if debugConnections==1:
            if debug==None:
                debug="mpiGet %d"%mpiComm.rank
        base.aobase.aobase.__init__(self,None,None,args=args,debug=debug,idstr=idstr)
        self.rank=mpiComm.rank
        self.generateNextTime=0.
        self.mpiParent=mpiParentObj
        self.outputData=self.mpiParent.data

        self.syncarr=numpy.ones((1,),numpy.uint8)
        self.mpiComm=mpiComm#Scientific.MPI.world.duplicate()

    def generateNext(self,msg=None):
        """called by the child when it wants more data.
        Get the data from the MPI connection and return it.
        msg is an instance of fwdMsg.fwdMsg().

        @param msg: The message to pass to the parent (predecessor) object
        over the MPI connection (or None).
        @type msg: None or fwdMsg object
        """
        t1=time.time()
        if self.debug!=None: print "mpiGet: generate=%d, Receiving MPI message from rank %d, tag %d (debug=%s)"%(self.generate,self.mpiParent.sourceRank,self.mpiParent.tag,str(self.debug))
        if self.debug!=None: print "mpiGet: Sending sync response (generate=%s) (debug=%s)"%(str(self.generate),str(self.debug))
        #try:
        self.mpiComm.send(self.syncarr,self.mpiParent.sourceRank,self.mpiParent.tag)
        #except Scientific.MPI.core.MPIError:
        #    print "MPI seems not to be configured for numpy - try Numeric instead."
        #    self.syncarr=Numeric.array(self.syncarr,copy=0)
        #    self.mpiComm.send(self.syncarr,self.mpiParent.sourceRank,self.mpiParent.tag)
        #    self.mpiParent.data=Numeric.array(self.mpiParent.data,copy=0)
        if self.debug!=None: print "mpiGet: Waiting to receive (debug=%s)"%str(self.debug)
        data,rank,tag,nelements=self.mpiComm.receive(self.mpiParent.data,self.mpiParent.sourceRank,self.mpiParent.tag)
        if self.debug!=None: print "mpiGet: Received MPI data (debug=%s)"%str(self.debug)
        #print "mpiGet: data=",data
        if nelements==0:#nelements is the number of elements of the array (not bytes) returned.
            self.dataValid=0
            sh="no data received"
        else:
            self.outputData=numpy.array(data,copy=0)
            self.dataValid=1
            sh=self.outputData.shape
        if self.debug!=None:
            print "mpiGet: Done generateNext, dataValid=%d shape=%s (debug=%s)"%(self.dataValid,str(sh),str(self.debug))
        self.generateNextTime=time.time()-t1

    def setGenerate(self,val):
        """set value of generate class variable
        @param val: Value to set generate to
        @type val: Int"""
        if self.debug!=None:
            print "mpiGet: setGenerate(%d) called (debug=%s)"%(val,str(self.debug))
        self.generate=val
        if val:
            self.syncarr[0]=1
        else:
            self.syncarr[0]=0

class mpiParent:
    """An object which creates the array to hold the data
    received by the mpi connection and also specifies the MPI process from
    which data is expected.
    @cvar sourceRank: The MPI rank of the source object
    @type sourceRank: Int
    @cvar tag: The tag expected, associated with the data.
    @type tag: Int
    @cvar data: The data array for holding the MPI data received.
    @type data: numpy array
    """
    def __init__(self,dims,dtype,sourceRank,tag):
        """Initialise the MPIParent object.
        @param dims: The dimensions of the numpy array expected from the parent (predecessor) object
        @type dims: Tuple of Int
        @param dtype: The data type of the numpy array
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
            self.data=numpy.zeros(dims,dtype)

def newMPIGet(dims,dtype,sourceRank,tag,mpiComm,args={},debug=None,idstr=None):
    """This function generates the mpiGet object (automatically creating the mpiParent object), and returns it to the user."""
    p=mpiParent(dims,dtype,sourceRank,tag)
    return mpiGet(p,mpiComm,args,debug,idstr=idstr)
