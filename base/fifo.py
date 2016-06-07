from __future__ import print_function
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
                print("INFORMATION::**fifo:{:s}**: delay = {:d}".format(
                         str(self.idstr), self.delay )
                    )
            except:
                self.shape=None
                self.dtype=None
                print("INFORMATION::**fifoDelay:{:s}**: not yet initialized".format(
                         str(self.idstr) )
                    )

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
                    else:
                        self.outputData=self.config.getVal("fifoShape"),self.config.getVal("fifoDtype")

    def newParent(self,parent,idstr=None):
        self.parent=parent
        self.initialised=0
        self.finalInitialisation()

    def generateNext(self,msg=None):
        """called by the child when it wants more data.
        @param msg: The message to pass to the parent (predecessor) object
        @type msg: None or fwdMsg object
        """
        if self.debug is not None:
            print("INFORMATION::**fifo:{:s}**: generateNext (debug={:s})".format(
                     str(self.idstr), str(self.debug) )
                )
        if self.generate==1:
            if self.newDataWaiting:
                self.dataValidList.append(self.parent.dataValid)
                if self.parent.dataValid==1:
                    if (self.spareArray is not None) and self.spareArray.shape==self.parent.outputData.shape and self.spareArray.dtype==self.parent.outputData.dtype:
                        #do this to avoid an unnecessary malloc...
                        self.spareArray[:]=self.parent.outputData
                        self.dataList.append(self.spareArray)
                        self.spareArray=None
                    else:
                        self.dataList.append(self.parent.outputData.copy())
                elif self.debug is not None:
                    print(("INFORMATION::**fifo:{:s}**: waiting for data but not valid "+
                           "(debug={:s})").format( str(self.idstr), self.debug )
                        )
                self.dataValid=self.dataValidList.pop(0)
                if self.dataValid:
                    self.outputData=self.dataList.pop(0)
                    if self.outputData is None:#make the array
                        if (self.spareArray is not None) and self.spareArray.shape==self.shape and self.spareArray.dtype==self.dtype:
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
   def test1():
       print("TEST1::")
       class dummy:
           dataValid=1
           outputData=numpy.zeros((2,2),"f")
       d=dummy()
       f=fifo(d,None,idstr="2")
       for i in range(1,10):
           d.outputData[:]=i
           #print i,f.outputData
           f.generateNext()
           print("\t"),
           print(i,f.outputData)
       #
       print("::ENDS")

   def test2():
       print("TEST2::")

#import numpy
#import base.fifo
       import base.readConfig
       class parent:
           dataValid=0
           outputData=numpy.arange(4)
       
       c=base.readConfig.AOXml([])
       p=parent()
       f=base.fifo.fifo(p,c,idstr="2")
       
       
       print("Running with a fifo delay of 2, and parent valid every 3 frames:")
       for i in range(12):#run 12 iterations
           p.dataValid=int((i%3)==2)
           if p.dataValid:
               p.outputData[:]=i
           f.doNextIter()
           print("\t")
           print(p.dataValid,p.outputData,f.dataValid,f.outputData)
       print("\n\n\nRunning with a fifo delay of 2 and parent valid every other frame:")
       f=base.fifo.fifo(p,c,idstr="2")
       for i in range(10):#run 10 iterations
           p.dataValid=i%2
           if p.dataValid:
               p.outputData[:]=i
           f.doNextIter()
           print("\t")
           print(p.dataValid,p.outputData,f.dataValid,f.outputData)
       print("\n\n\nRunning with a fifo delay of 1 and parent valid every 3 frames:")
       f=base.fifo.fifo(p,c,idstr="1")
       for i in range(10):#run 10 iterations
           p.dataValid=int((i%3)==2)
           if p.dataValid:
               p.outputData[:]=i
           f.doNextIter()
           print("\t")
           print(p.dataValid,p.outputData,f.dataValid,f.outputData)
       
       print("\n\n\nRunning with a fifo delay of 1 and parent valid every 2 frames (ie every other):")
       f=base.fifo.fifo(p,c,idstr="1")
       for i in range(10):#run 10 iterations
           p.dataValid=i%2
           if p.dataValid:
               p.outputData[:]=i
           f.doNextIter()
           print("\t")
           print(p.dataValid,p.outputData,f.dataValid,f.outputData)
       #
       print("::ENDS")

   test1()
   test2()
