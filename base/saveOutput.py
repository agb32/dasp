#$Id: saveOutput.py,v 1.10 2011/01/20 08:08:36 ali Exp $
import types,os
import numpy
import cmod.utils
import base.aobase
import util.FITS
class saveOutput(base.aobase.aobase):
    """
    A class for saving the output to a FITS file.  Should you wish to save only part of the output, or
    the output of several different objects, use splitOutput or joinOutput first.
    
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
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.ended=0
        self.initialised=0
        self.finalised=0
        self.outputData=None
        self.forGUISetup=forGUISetup
        if args.has_key("filename"):
            self.filename=args["filename"]
        else:
            if len(self.idstr)==0 or self.idstr[0]==None:
                ids=""
            else:
                ids="_%s"%self.idstr[0]
            self.filename=self.config.getVal("saveOutputFilename",default="saveOutput%sb%d.fits"%(ids,config.batchno))
        self.doByteSwap=self.config.getVal("doByteSwap",default=0)
        self.dataValid=0#never changes... no output...
        if forGUISetup:
            self.outputData=[None,None]

    def __del__(self):
        self.endSim()

    def endSim(self):
        if self.ended==0:
            self.ended=1
            if self.finalised==0:
                self.finalised=1
                print "saveOutput: Finalising FITS file %s"%self.filename
                import stat
                try:
                    s=os.stat(self.filename)
                    size=s[stat.ST_SIZE]
                except:
                    print "Unable to stat file %s"%self.filename
                    size=0
                #except:
                #    size=0
                tmp=size%2880
                if tmp>0:
                    try:
                        os.ftruncate(self.ff.fileno(),size+2880-tmp)
                    except:
                        print "saveOutput: unable to finalise FITS file"
                self.ff.close()

    def finalInitialisation(self):
        """once we've got all the info we need, initialise here..."""
        if self.initialised==0:
            self.initialised=1
            self.shape=self.parent.outputData.shape
            print self.idstr,self.parent
            self.dtype=self.parent.outputData.dtype
            if type(self.dtype)!=type(""):
                self.dtype=self.dtype.char
            if self.dtype=="f":
                self.bitpix=-32
            elif self.dtype=="d":
                self.bitpix=-64
            elif self.dtype=="i":
                self.bitpix=32
            elif self.dtype=="s":
                self.bitpix=16
            elif self.dtype=="1":
                self.bitpix=8
            else:
                raise Exception("base.saveOutput - bitpix not known (please update)")
            try:
                self.ff=open(self.filename,"wb+",buffering=4096)
            except:#there seemed to be a problem when the python code was called on a ext3 mount from a directory in a cifs mount... no idea why!  The solution is to copy files to where you run them from...
                print "Unable to open file with buffering - trying without"
                try:
                    self.ff=open(self.filename,"wb+")
                except:
                    print "Unable to open file in wb+ mode.  Trying w+."
                    try:
                        self.ff=open(self.filename,"w+")
                    except:
                        print "Unable to open file in w+ mode.  Trying w."
                        self.ff=open(self.filename,"w")
                
            #now write the header...
            self.nd=len(self.shape)+1
            self.dims=list(self.shape[:])
            self.dims.reverse()
            self.dims+=[0]
            util.FITS.WriteKey(self.ff,"SIMPLE","T")
            util.FITS.WriteKey(self.ff,"BITPIX",self.bitpix)
            util.FITS.WriteKey(self.ff,"NAXIS",str(self.nd))
            for i in range(self.nd):
                util.FITS.WriteKey(self.ff,"NAXIS%d"%(i+1),str(self.dims[i]))
            self.axisIncPos=self.ff.tell()-80
            util.FITS.WriteKey(self.ff,"EXTEND","T")
            le=numpy.little_endian and self.doByteSwap==0
            util.FITS.WriteKey(self.ff,"LITTLE_E",le,"whether file is little endian")
            if le:
                util.FITS.WriteKey(self.ff,"UNORDERD","T")
                util.FITS.WriteComment(self.ff,"Note, this file is be non-FITS complient, saved on a little endian machine.  All the data will be byteswapped.")
            util.FITS.EndHeader(self.ff)

    def newParent(self,parent,idstr=None):
        self.parent=parent
        self.initialised=0

    
    def generateNext(self,msg=None):
        """called by the child when it wants more data.
        @param msg: The message to pass to the parent (predecessor) object
        @type msg: None or fwdMsg object
        """
        if self.debug!=None:
            print "splitOutput: generateNext (debug=%s)"%str(self.debug)
        if self.generate==1:
            if self.newDataWaiting:
                if self.parent.dataValid==1:
                    outputData=self.parent.outputData
                    if outputData.shape!=self.shape:
                        print "ERROR: saveOutput - outputdata not same as that given %s %s"%(str(outputData.shape),str(self.shape))
                        raise Exception("Wrong shape in splitOutput")
                    if outputData.dtype.char!=self.dtype:
                        outputData=outputData.astype(self.dtype)
                        print "Warning: saveOutput - converting outputData to type %s"%self.dtype
                    if self.doByteSwap and numpy.little_endian:
                        outputData=outputData.byteswap()
                    self.ff.seek(0,2)#move to end of file.
                    self.ff.write(outputData)
                    self.ff.seek(self.axisIncPos)
                    key=self.ff.read(80)
                    self.ff.seek(self.axisIncPos)
                    newdim=int(key[10:])+1
                    util.FITS.WriteKey(self.ff,"NAXIS%d"%self.nd,str(newdim))
                else:
                    print "saveOutput: waiting for data but not valid (debug=%s)"%self.debug
                    self.dataValid=0
        else:
            self.dataValid=0

if __name__=="__main__":
    class dummy:
        outputData=numpy.zeros((5,10),"i")
        dataValid=1
    class dumconfig:
        def __init__(self):
            pass
        def getVal(self,val,default=None,raiseerror=1):
            return "tmp.fits"
        def setSearchOrder(self,so):
            pass
    parent=dummy()
    config=dumconfig()
    so=saveOutput(parent,config)
    so.finalInitialisation()
    for i in range(14):
        parent.outputData[:]=i
        so.generateNext()
