#$Id: mpiSend.py,v 1.13 2009/02/18 10:28:41 ali Exp $
#import Scientific.MPI#,threading
import numpy
import types,time
import base.aobase

debugConnections=0
class mpiSend(base.aobase.aobase):
    """Send data over an MPI connection to the child
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
     - objID - string, the object identifier
     - debug - if not None will print debug messages
    Class variables (not important to simulation programmer):
     - generate - whether it will generate this iteration (usually set/unset by successor (child) object).
     - dataValid - whether data is valid.
     - newDataWaiting - whether new data is waiting this iteration.
     - syncmsg - synchronisation message
    @cvar parent: Predecessor (parent) object from which data is obtained.
    @type parent: Science object
    @cvar mpiChild: Structure detailing MPI connection.
    @type mpiChild: MPIChild object
    @cvar mpiComm:  MPI communicator object
    @type mpiComm: Scientific.MPI.world
    @cvar generate: Whether to generate this iteration
    @type generate: Int
    @cvar dataValid: Whether data is valid
    @type dataValid: Int
    @cvar newDataWaiting: Whether to collect data from parent (predecessor) object this iteration.
    @type newDataWaiting: Int
    @cvar objID: Object identifier
    @type objID: String
    @cvar debug: Flag, whether to print debug message (if None, won't print)
    @type debug: None or user defined.
    @cvar syncmsg: Synchronisation message
    @type syncmsg: Array
    """
    #def __init__(self,parent,mpiChildObj,mpiComm,mpiEngine,mpiOrder,useFwdMsg=0,args={},debug=None):
    def __init__(self,parent,mpiChildObj,mpiComm,args={},debug=None,idstr=None):
        """Initialise the MPI connection.
        @param parent: Science object from which data is obtained
        @type  parent: Object
        @param mpiChildObj: Details of MPI connection
        @type  mpiChildObj: MPIChild
        @param mpiComm: Duplicate of Scientific.MPI.world
        @type mpiComm: Object
        @param args: Optional argument dictionary, can have keys of idstr
        @type args: Dict
        @param debug: Flag, whether to print debug message (if None, won't print)
        @type debug: None or user defined.
        """
        print "remoteMPISend: Initialising"
        if debugConnections==1:
            if debug==None:
                debug="mpiSend %d"%mpiComm.rank
        base.aobase.aobase.__init__(self,parent,None,args=args,debug=debug,idstr=idstr)
        self.generateNextTime=0.
        #self.iter=0
        #self.lock=threading.Lock()
        self.mpiChild=mpiChildObj
        #self.useFwdMsg=useFwdMsg
        #self.mpiEngine=mpiEngine
        #self.realmpiOrder=mpiOrder
        #if self.useFwdMsg:
        #    self.mpiOrder=2*self.realmpiOrder+1
        #    self.mpiSyncOrder=2*self.realmpiOrder
        #else:
        #    self.mpiOrder=2*self.realmpiOrder
        #    self.mpiSyncOrder=2*self.realmpiOrder+1
        #print "remoteMPISend: creating mpi.world.duplicate()"
        self.mpiComm=mpiComm#Scientific.MPI.world.duplicate()
        self.syncmsg=numpy.zeros((1,),numpy.uint8)
            
        print "remoteMPISend: Initialised"
    def generateNext(self,msg=None):
        """Obtain the next data from the parent object, and pass over MPI
        connection to successor object.
        @param msg: Does nothing
        @type msg: None
        """
        t1=time.time()
        #self.lock.acquire()#must happen atomically
        #thisiter=self.iter
        #self.iter+=1
        #self.lock.release()
        #if self.useFwdMsg:#wait for the message to arrive.
        #    self.mpiEngine.myturn(self.mpiSyncOrder,thisiter)
        #    self.mpiComm.receive(self.msg.arr,self.mpiChild.rank,self.mpiChild.tag)
        #    self.mpiEngine.endturn()
        #    self.msg.fromArray()#populate msg.i, msg.f, msg.strng
        if self.debug!=None: print "mpiSend: in GenerateNext (debug=%s)"%str(self.debug)
        if self.parent.dataValid==1:#store the data
            #self.outputDataList.append(self.parent.outputData.copy())

            self.outputData=self.parent.outputData
            self.dataValid=1
        else:
            if self.debug!=None: print "MPISend: parent data not valid (debug=%s)"%str(self.debug)
            self.dataValid=0
            #self.outputDataList.append("")

        for i in range(len(self.mpiChild.rank)):
            rank=self.mpiChild.rank[i]
            tag=self.mpiChild.tag[i]
            if self.debug!=None: print "mpiSend: Receiving synchronisation rank %d tag %s (debug=%s)"%(rank,str(tag),str(self.debug))
            self.mpiComm.receive(self.syncmsg,rank,tag)
            if self.syncmsg[0]==1:#generate=1
                self.parent.setGenerate(1)
                #data=self.outputDataList.pop(0)
                data=self.outputData
                sh=data.shape
            else:#self.outputData stored for next time... assuming next time, parent.dataValid will be zero.
                self.parent.setGenerate(0)
                data=""
                sh="data=''"
            if self.debug!=None: print "mpiSend: (got generate=%s) Sending data rank %d tag %s shape %s (debug=%s)"%(str(self.syncmsg[0]),rank,str(tag),str(sh),str(self.debug))
            self.mpiComm.send(data,rank,tag)#blocking
            if self.debug!=None: print "mpiSend: data sent rank %d tag %s (debug=%s)"%(rank,str(tag),str(self.debug))
        self.generateNextTime=time.time()-t1
class mpiChild:
    """Structure class holding information about connecting object via MPI.
    @cvar rank: MPI rank of the connecting object (to which data is sent).
    @type rank: Int or List
    @cvar tag: Tag sent with data to connecting object
    @type tag: Int or List
    """
    def __init__(self,rank,tag):
        """Rank is the rank of the child process.  tag is a tag used
        to identify which signals we want.  Rank and tag can be single integers of lists.  If one is a list and the other an integer, the integer is made into a list of the same length."""
        if type(rank)==types.IntType:
            if type(tag)==types.IntType:
                self.rank=[rank]
                self.tag=[tag]
            elif type(tag)==types.ListType:
                self.rank=[rank]*len(tag)
                self.tag=tag
            else:
                raise Exception("MPISend tag must be list or int")
        elif type(rank)==types.ListType:
            if type(tag)==types.IntType:
                self.rank=rank
                self.tag=[tag]*len(rank)
            elif type(tag)==types.ListType:
                if len(rank)!=len(tag):
                    raise Exception("MPISend tag and rank lists must be same length")
                self.rank=rank
                self.tag=tag
            else:
                raise Exception("MPISend tag must be list or int")
        else:
            raise Exception("MPISend rank must be list or int")

def newMPISend(parent,rank,tag,mpiComm,args={},debug=None,idstr=None):
    """This function generates the mpiSend object (automatically creating the mpiChild object), and returns it to the user."""
    c=mpiChild(rank,tag)
    return mpiSend(parent,c,mpiComm,args,debug,idstr=idstr)
