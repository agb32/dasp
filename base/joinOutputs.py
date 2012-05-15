#$Id: joinOutputs.py,v 1.8 2011/12/15 07:26:25 ali Exp $
import types
import base.fwdMsg
#import Numeric
import numpy
import cmod.utils
import base.aobase
class joinOutputs(base.aobase.aobase):
    """
    A class for joining the outputs from a science module previously split up.  This class returns the full parent output array.  This can be useful for parallelising a science module - eg DM phase is split into 4 before being send to 4 separate wfscent, each of which compute centroids for this part of the DM phase, and then the outputs are rejoined back together here.
    @cvar parent: the parent
    @type parent: object
    @cvar config: config object
    @type config: object
    @cvar args: dict of arguements which can include shape, dtype (shape and datatype of the joined up array), and codeDict, a dictionary of strings (one for each parent object), which are used to place the parent output into the output here.
    @type args: dict
    @cvar debug: Flag, whether to print debug message (if None, won't print)
    @type debug: None or user defined.
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Initialise the splitOutput module. 
        @param parent: Parent object which generates the data
        @type parent: Parent object
        @param config: Configuration object
        @type config: Object
        @param args: Dictionary of arguments, at present can contain keys: idstr
        @type args: Dict
        @param forGUISetup: whether for GUI or not
        @type forGUISetup: Int
        @param debug: Flag, whether to print debug message (if None, won't print)
        @type debug: None or user defined.
        """
        if type(parent)!=types.DictType:
            parent={1:parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        if args.has_key("shape"):
            self.shape=args["shape"]
        else:
            self.shape=self.config.getVal("joinOutputsShape")
        self.dtype=args["dtype"]
        if type(self.shape)==type(""):
            self.shape=eval(self.shape)
        self.nParents=len(self.parent.keys())
        self.zeroParent=1
        if args.has_key("zero"):
            self.zeroParent=args["zero"]

##         if args.has_key("shape"):#the shape of the output array...
##             self.shape=args["shape"]
##         else:
##             self.shape=self.config.getVal("joinOutputShape")
##         if args.has_key("dtype"):
##             self.dtype=args["dtype"]
##         else:
##             self.dtype=self.config.getVal("joinOutputDtype")
##         if args.has_key("parentInfo"):
##             self.parentInfo=args["parentInfo"]
##         else:
##             self.parentInfo=self.config.getVal("joinOutputParentInfo")
        #parent info should be a dictionary with keys equal to the keys for self.parent.  The values of these is then a tuple of 2 tuples giving details about how to copy the parent outputData into the outputData here.  Each of these tuples contains tuples of (start, end, step) for the output data here, and the parent outputdata.  So, for example, if you had (((0,6,2),(0,10,1)),((0,3,1),(10,30,2))) this would be the same as self.outputData[0:6:2,0:10:1]=parent.outputData[0:3:1,10:30:2]
            
##         self.itemsize=self.outputData.itemsize()
##         self.arrsize=self.itemsize
##         for i in self.shape:
##             self.arrsize*=i
        #self.endOffset=self.startOffset+self.arrsize
        if forGUISetup==1:
            self.outputData=[self.shape,self.dtype]
        else:
            self.outputData=numpy.zeros(self.shape,self.dtype)
            self.codeDict=args["codeDict"]
            # codeDict should have a entry for each parent.  This is then the python that gets exec'd to copy the parent data into the output data array.  It should use parentData and outputData for this purpose.  eg to put the data into half of the outputData array, you could use "outputData[:outputData.shape[0]]=parentData", or more simply, "outputData[:50]=parentData" etc.  Note, the output array will exist, and must not be overwritten - should be written into only.
            self.compiledCodeDict={}
            for key in self.codeDict.keys():
                self.compiledCodeDict[key]=compile(self.codeDict[key],"<string>","exec")


    def newParent(self,parent):
        raise Exception("Please give parent upon initialisation")

    def generateNext(self,msg=None):
        """called by the child when it wants more data.
        @param msg: The message to pass to the parent (predecessor) object
        @type msg: None or fwdMsg object
        """
        if self.debug!=None:
            print "splitOutput: generateNext (debug=%s)"%str(self.debug)
        if self.generate==1:
            if self.newDataWaiting:
                nin=0
                if self.zeroParent:
                    self.outputData[:]=0
                for key in self.parent.keys():
                    if self.parent[key].dataValid==1:
                        d={"outputData":self.outputData,"parentData":self.parent[key].outputData,"parent":self.parent[key],"numpy":numpy}
                        try:
                            exec self.compiledCodeDict[key] in d
                        except:
                            print "ERROR: joinOutputs - copying outputs for key %s: %s"%(str(key),self.codeDict[key])
                            print "output, parent shape: %s %s"%(str(self.outputData.shape),str(self.parent[key].outputData.shape))
                            raise
                        nin+=1
                if nin>0:
                    self.dataValid=1
                    if nin!=self.nParents:
                        print "joinOutputs: Received unexpected number of datas from parents (%d out of %d)"%(nin,self.nParents)
                else:
                    self.dataValid=0
        else:
            self.dataValid=0

