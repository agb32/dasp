#$Id: splitOutput.py,v 1.15 2011/05/25 05:50:45 ali Exp $
import types
import base.fwdMsg
#import Numeric
#import numpy
import cmod.utils
import base.aobase
class splitOutput(base.aobase.aobase):
    """
    A class for splitting the output from a science module up.  This class just returns part of the parent output array.
    @cvar parent: the parent
    @type parent: object
    @cvar config: config object
    @type config: object
    @cvar args: dict of arguements
    @type args: dict
    @cvar startOffset: Offset into array to read freom
    @type startOffset: Int
    @cvar shape: Shape of the output array
    @type shape: Tuple of Ints
    @cvar dtype: Datatype of output array
    @type dtype: Char
    @cvar debug: Flag, whether to print debug message (if None, won't print)
    @type debug: None or user defined.
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

        self.initialised=0
        self.code=config.getVal("splitOutputCode",default=args.get("code"))
        #Code must be python code that extracts part of an array from parentData, and writes it to outputData, e.g. "outputData=parentData[0:10]".  This can be as complex as you like.
        self.outputData=None
        self.compiledCode=compile(self.code,"<string>","exec")
        self.shape=None
        self.dtype=None
        self.feedback=0#not a feedback object by default...
        self.forGUISetup=forGUISetup
        if args.has_key("outputData"):
            self.outputData=args["outputData"]
            self.shape=self.outputData.shape
            self.dtype=self.outputData.dtype.char#typecode()

        if args.has_key("shape"):
            self.shape=args["shape"]#the shape of the final array.
        elif self.shape==None:
            self.shape="fromparent"
        #print "dtyp1e/shape:%s %s"%(str(self.shape),str(self.dtype))
        if args.has_key("dtype"):
            self.dtype=args["dtype"]#the datatype of the final array.
        elif self.dtype==None and not hasattr(self.dtype,"char"):#note numpy.flaot32 etc evaluate to None for some reason...
            self.dtype="fromparent"
        self.makecontiguous=0
        if args.has_key("makecontiguous"):
            self.makecontiguous=args["makecontiguous"]
        if args.has_key("feedback"):
            self.feedback=args["feedback"]
        if self.feedback:
            self.dataValid=1
        #print "dtype/shape:%s %s"%(str(self.shape),str(self.dtype))
        try: #attempt to initialise...
            self.finalInitialisation()
        except: #probably means the parent will be supplied later...
            print "Cannot yet initialise splitOutput object."
        if self.forGUISetup:
            if type(self.dtype)!=type(""):
                dtype=self.dtype.char
            else:
                dtype=self.dtype
            if len(dtype)!=1:
                print "WARNING - splitOutput dtype='%s'"%dtype
            self.outputData=[self.shape,dtype]
    def finalInitialisation(self):
        """once we've got all the info we need, initialise here..."""
        if self.initialised==0:
            dtype=None
            if self.shape=="fromparent":#attempt to get shape from parent (could be dangerous if parent not yet created properly).
                d={"parentData":self.parent.outputData,"outputData":self.outputData,"parent":self.parent}
                exec self.compiledCode in d
                self.shape=d["outputData"].shape
                print "splitOutput got shape of %s from parent (code=%s)"%(self.shape,self.code)
                if self.dtype=="fromparent":
                    dtype=d["outputData"].dtype#typecode()
            #print self.dtype
            if type(self.dtype)==type("") and self.dtype=="fromparent":
                if dtype==None:
                    self.dtype=self.parent.outputData.dtype.char#typecode()
                else:
                    self.dtype=dtype
            self.ds=None
            if type(self.outputData)==type(None):
                class dummyshape:
                    shape=self.shape
                    dtype=self.dtype
                    #def dtype(self):
                    #    return self.dtype
                self.outputData=dummyshape()#just so other things can get the shape.  Note, we don't want to create the array here, cos it might be large/unused etc.
                self.ds=self.outputData
            #self.outputData=Numeric.zeros((1,),self.dtype)
            #self.itemsize=self.outputData.itemsize()
            #self.parentMem=self.parent.outputData.itemsize()*reduce(lambda x,y:x*y,self.parent.outputData.shape)
            #self.mymem=reduce(lambda x,y:x*y,self.shape)*self.outputData.itemsize()
            #self.endOffset=self.startOffset+self.mymem
            if self.forGUISetup==1:
                self.outputData=[self.shape,self.dtype]
            else:
                pass

            self.initialised=1

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
            print "splitOutput: generateNext (debug=%s)"%str(self.debug)
        if self.ds is self.outputData:#output data is only dummy, so we set to None here, which means it will be created.
            self.outputData=None
        if self.generate==1:
            if self.newDataWaiting:
                if self.parent.dataValid==1:
                    self.dataValid=1
                    d={"parentData":self.parent.outputData,"outputData":self.outputData,"parent":self.parent}
                    exec self.compiledCode in d
                    outputData=d["outputData"]
                    if outputData.shape!=self.shape:
                        print "ERROR: splitOutput - outputdata after exec is not same as that given %s %s"%(str(outputData.shape),str(self.shape))
                        raise Exception("Wrong shape in splitOutput")
                    if outputData.dtype!=self.dtype:
                        outputData=outputData.astype(self.dtype)
                    if self.makecontiguous and not outputData.iscontiguous():
                        outputData=outputData.copy()
                    self.outputData=outputData
                    
##                     if !self.parent.outputData.iscontiguous():
##                         outputData=Numeric.array(self.parent.outputData)
##                     else:
##                         outputData=self.parent.outputData
##                     #first take memory from the parent offset.
##                     outputData=cmod.utils.arrayFromArray(outputData,((self.parentMem-self.startOffset)/self.itemsize,),self.dtype,self.startOffset)

                                                             
##                     self.outputData=cmod.utils.arrayFromArray(self.parent.outputData[self.startOffset:self.endOffset],self.shape,self.dtype)
                else:
                    print "splitOutput: waiting for data but not valid (debug=%s)"%debug
                    self.dataValid=0
        else:
            self.dataValid=0

