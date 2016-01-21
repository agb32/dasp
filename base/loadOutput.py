#$Id: loadOutput.py,v 1.5 2011/08/18 05:12:19 ali Exp $
import types,os
import numpy
import cmod.utils
import base.aobase
import util.FITS
class loadOutput(base.aobase.aobase):
    """
    A class for loading a fits file into the output at each iteration.
    Once the fits file has been exhausted, it will print a warning and start from the beginning again.  Typically,
    this is used on a file that has been saveoutput'd before.
    
    @cvar parent: the parent
    @type parent: object
    @cvar config: config object
    @type config: object
    @cvar args: dict of arguements
    @type args: dict
    @cvar shape: Shape of the output array
    @type shape: Tuple of Ints
    @cvar dtype: Datatype of output array
    @type dtype: Char
    @cvar debug: Flag, whether to print debug message (if None, won't print)
    @type debug: None or user defined.
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        """Initialise the loadOutput module. 
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
        if parent!=None:
            print "WARNING: base.loadOutput - parent will be ignored"
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)

        self.initialised=0
        self.outputData=None
        self.forGUISetup=forGUISetup
        if args.has_key("filename"):
            self.filename=args["filename"]
        else:
            self.filename=self.config.getVal("loadOutputFilename")
        self.header=util.FITS.ReadHeader(self.filename)
        self.dataOffset=len(self.header["raw"])*80
        if self.dataOffset%2880!=0:
            raise Exception("base.loadOutput - dataOffset not multiple of 2880 for file %s"%filename)
        self.header=self.header["parsed"]
        self.shape=[]
        for i in range(int(self.header["NAXIS"])-1):
            self.shape.insert(0,int(self.header["NAXIS%d"%(i+1)]))
        self.shape=tuple(self.shape)
        self.nd=int(self.header["NAXIS%d"%(int(self.header["NAXIS"]))])#this is the number of iterations for which we have data...
        self.bitpix=int(self.header["BITPIX"])
        if self.bitpix==-64:
            self.dtype="d"
        elif self.bitpix==-32:
            self.dtype="f"
        elif self.bitpix==32:
            self.dtype="i"
        elif self.bitpix==16:
            self.dtype=="s"
        elif self.bitpix==8:
            self.dtype="1"
        else:
            raise Exception("base.loadOutput - bitpix not known (please update)")
        self.dataSize=reduce(lambda x,y:x*y,self.shape)*numpy.fabs(self.bitpix)/8
        self.dataValid=0#never changes... no output...
        if forGUISetup:
            self.outputData=[self.shape,self.dtype]
        else:
            self.outputData=numpy.zeros(self.shape,self.dtype)
        #self.ff=open(self.filename,"r")
        #self.ff.seek(self.dataOffset)
        self.map=numpy.memmap(self.filename,dtype=self.dtype,mode="r")
        self.dataOffset/=self.map.itemsize
        self.littleEndianData=int(eval(self.header.get("LITTLE_E",0)))
        self.doByteSwap=(self.littleEndianData!=numpy.little_endian)
        self.dataCnt=0
   
    def generateNext(self,msg=None):
        """called by the child when it wants more data.
        @param msg: The message to pass to the parent (predecessor) object
        @type msg: None or fwdMsg object
        """
        if self.debug!=None:
            print "splitOutput: generateNext (debug=%s)"%str(self.debug)
        if self.generate==1:
            self.dataValid=1
            data=self.map[self.dataOffset+self.dataCnt*self.dataSize/self.map.itemsize:self.dataOffset+(self.dataCnt+1)*self.dataSize/self.map.itemsize]
            self.outputData.flat[:]=data
            #data=self.ff.read(self.dataSize)
            #self.outputData.flat[:]=numpy.fromstring(data,dtype=self.dtype)
            if self.doByteSwap:
                self.outputData[:]=self.outputData.byteswap()
            self.dataCnt+=1
            if self.dataCnt>=self.nd:
                self.dataCnt=0
                print "WARNING: loadOutput has reached end of data - wrapping round."
                #self.ff.seek(self.dataOffset)
        else:
            self.dataValid=0

if __name__=="__main__":
    class dumconfig:
        def __init__(self):
            pass
        def getVal(self,val,default=None,raiseerror=1):
            return "tmp.fits"
        def setSearchOrder(self,so):
            pass
    config=dumconfig()
    so=loadOutput(None,config)
    for i in range(28):
        so.generateNext()
        print so.outputData
