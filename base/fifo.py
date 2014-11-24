import base.aobase
import numpy
class fifo(base.aobase.aobase):
    """
    A class for delaying an output by an arbitrary number of iterations (to facilitate AO loop delay/latency).  The number of cycles delay can be specified as the idstr, or in param file as fifoDelay.
    @cvar parent: the parent
    @type parent: object
    @cvar config: config object
    @type config: object
    @cvar args: dict of arguements
    @type args: dict
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Initialise the splitOutput module. 
        @param parent: Parent object which generates the data
        @type parent: Parent object
        @param config: Configuration object
        @type config: Object
        @param args: Dictionary of arguments, at present can contain keys: idstr, code, shape, dtype, makecontiguous.   Code must be a string of python code that extracts part of an array from parentData, and writes it to outputData, e.g. "outputData=parentData[0:10]".  This can be as complex as you like.  shape is the shape of the output array, and dtype is the data type of it.  If shape doesn't agree with what code produces, an error will be raised.  dtype can be "fromparent".
        @type args: Dict
        @param forGUISetup: whether for GUI or not
        @type forGUISetup: Int
        @param debug: Flag, whether to print debug message (if None, won't print)
        @type debug: None or user defined.
        """
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)

        try:
            delay=int(idstr)
        except:
            delay=1
        delay=self.config.getVal("fifoDelay",default=delay)
        self.delay=delay
        self.initialised=0
        self.outputData=None
        self.dataList=[None]*self.delay
        self.dataValidList=[1]*self.delay
        self.forGUISetup=forGUISetup
        self.spareArray=None
        self.args=args
        self.finalInitialisation()#this may fail...

    def finalInitialisation(self):
        """once we've got all the info we need, initialise here..."""
        if self.initialised==0:
            try:
                self.shape=self.parent.outputData.shape
                self.dtype=self.parent.outputData.dtype
                self.initialised=1
            except:
                self.shape=None
                self.dtype=None
                print "base.fifo object not yet initialised"

            if self.initialised:
                if self.forGUISetup:
                    self.outputData=(self.shape,self.dtype)
                else:
                    self.outputData=numpy.zeros(self.shape,self.dtype)
            else:
                if self.forGUISetup:
                    #Needs to know shape and dtype so this can be returned...
                    if self.args.has_key("shape") and self.args.has_key("dtype"):
                        self.outputData=self.args["shape"],self.args["dtype"]

    def newParent(self,parent,idstr=None):
        self.parent=parent
        self.initialised=0
        self.finalInitialisation()

    def generateNext(self,msg=None):
        """called by the child when it wants more data.
        @param msg: The message to pass to the parent (predecessor) object
        @type msg: None or fwdMsg object
        """
        if self.debug!=None:
            print "fifo: generateNext (debug=%s)"%str(self.debug)
        if self.generate==1:
            if self.newDataWaiting:
                self.dataValidList.append(self.parent.dataValid)
                if self.parent.dataValid==1:
                    if self.spareArray!=None and self.spareArray.shape==self.parent.outputData.shape and self.spareArray.dtype==self.parent.outputData.dtype:
                        #do this to avoid an unnecessary malloc...
                        self.spareArray[:]=self.parent.outputData
                        self.dataList.append(self.spareArray)
                        self.spareArray=None
                    else:
                        self.dataList.append(self.parent.outputData.copy())
                else:
                    print "fifo: waiting for data but not valid (debug=%s)"%self.debug
                self.dataValid=self.dataValidList.pop(0)
                if self.dataValid:
                    self.outputData=self.dataList.pop(0)
                    if self.outputData==None:#make the array
                        if self.spareArray!=None and self.spareArray.shape==self.shape and self.spareArray.dtype==self.dtype:
                            self.outputData=self.spareArray
                            self.outputData[:]=0
                        else:
                            self.outputData=numpy.zeros((self.shape),self.dtype)
                    self.spareArray=self.outputData#the outputData memory can be reused next time...
                else:
                    self.spareArray=None
        else:
            self.dataValid=0


if __name__=="__main__":
    class dummy:
        dataValid=1
        outputData=numpy.zeros((2,2),"f")
    d=dummy()
    f=fifo(d,None,idstr="2")
    for i in range(1,10):
        d.outputData[:]=i
        #print i,f.outputData
        f.generateNext()
        print i,f.outputData
