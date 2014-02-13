#$Id: Ctrl.py,v 1.77 2011/12/06 05:14:39 ali Exp $
import types,time,string,select

import numpy
import sys,thread,threading,os,socket
import getopt,re
import base.readConfig,util.SockConn,cmod.shmem,util.serialise
import Scientific.MPI
class Ctrl:
    """A control class for the AO simulation.  This is responsible for keeping
    the simulation running.

    Class variables (important to simulation programmer):
     - globals - The globals() dictionary from the top level python file.
     - config - The configuration file object
    Class variables (not important to simulation programmer):
     - All the rest.
    @cvar go: Set to zero to end the simulation
    @type go: Int
    @cvar debug: Debug message if debugging (None if not).
    @type debug: None or user defined
    @cvar niter: number of iterations to run for (-1 means forever)
    @type niter: Int
    @cvar thisiter: Current iteration number
    @type thisiter: Int
    @cvar frametime: Time to compute last iteration
    @type frametime: Float
    @cvar paused: Whether simulation paused or running
    @type paused: Int
    @cvar nextniters: Determines how many iterations to run for before pausing.  If None, will continue running until simulation ends (or is paused).
    @type nextniters: Int
    @cvar globals: Dict of global variables.  Usually A Ctrl instance will be created by Ctrl(globals()).
    @type globals: Dict
    @cvar indent: The indent text to use when printing the simulation variable tree
    @type indent: String
    @cvar batchno: Current batch number
    @type batchno: Int
    @cvar rank: MPI rank
    @type rank: Int
    @cvar paramfile: XML Parameter filename
    @type paramfile: String
    @cvar config: A AOXml object holding simulation configuration information
    @type config: AOXml instance
    @cvar sockConn: A SockConn object which is listening for external instructions (e.g. from a GUI).
    @type sockConn: SockConn instance
    """
    def __init__(self,globals=None,paramfile=[],debug=None):
        """Initialise the Ctrl object.
        @param globals: The globals() dictionary
        @type globals: Dict
        @param debug: Debug message or None
        @type debug: None or user defined.
        """

        ## INITIALISATION of variables:
        self.simInitTime=time.time()
        self.debug=debug
        self.go=1
        self.niter=-1
        self.thisiter=0
        self.frametime=0
        self.paused=0
        self.nextniters=None
        self.globals=globals
        self.indent="    "
        self.batchno=0
        self.mpiComm=Scientific.MPI.world.duplicate()
        self.listenSTDIN=1
        self.rank=self.mpiComm.rank
        self.paramfile=paramfile
        self.simID=""
        self.compList=[]
        self.compListPos=0
        self.compListNames=[]
        #broadcast a tag for SHM from rank 0 to all others...
        arr=numpy.array([time.time()],numpy.float64)
        self.mpiComm.broadcast(arr,0)
        self.shmtag=long(arr.view(numpy.int64)[0])
        self.initCmdList=[]#text commands that are run just before the main loop is started.
        self.slowMotion=0#whether to pause between each iteration (eg for viewing gfx)...
        nice=19
        self.paramString=None
        self.initParamString=None
        sys.stdout=myStdout(self.rank)

        ## OPTION parsing (options from the command line):
        optlist,arglist=getopt.gnu_getopt(sys.argv[1:],"hp",["batchno=","start-paused","param-file=","help","iterations=","id=","nonice","param=","nostdin","debug-connections","user=","init="])
        for o, a in optlist:
            if o=="--displaythread":
                sys.stdout.displayThread=1
            if o=="--batchno":
                self.batchno=int(a)
                print "Batch number:",self.batchno
            if o in ["-p","--start-paused"]:
                self.paused=1
            if o=="--param-file":
                self.paramfile+=a.split(",")
                print "Using %s"%self.paramfile
            if o=="--nonice":
                nice=0
            if o=="--nostdin":
                self.listenSTDIN=0
            if o=="--debug-connections":
                self.debugConnections=1
                if globals!=None:
                    print globals
                    globals["base"].mpiGet.debugConnections=1
                    globals["base"].mpiSend.debugConnections=1
                    globals["base"].shmGet.debugConnections=1
                    globals["base"].shmSend.debugConnections=1
            if o=="--iterations":
                self.niter=eval(a)
                if type(self.niter)!=types.IntType and type(self.niter)!=types.ListType:
                    print "Warning, must use integer or list for number of iterations"
                    self.niter=-1
                if type(self.niter)==types.ListType:
                    if len(self.niter)!=self.mpiComm.size:
                        print "Warning, list length must be equal to number of MPI processes"
                    if self.rank>=len(self.niter):
                        self.niter=-1
                    else:
                        self.niter=self.niter[self.rank]
            if o=="--id":
                self.simID=a
            if o=="--param":#optional commandline to change parameters.
                self.paramString=a
            if o=="--init":
                self.initParamString=a
            if o in ["-h","--help"]:
                print 'HELP:\nRun simulation with\n--batchno=xxx\n--start-paused\n--param-file=paramfile\n--iterations=niters\n--id=simulationID (string)\n--param="Text string, e.g. this.globals.nLayers=2, used to alter parameter file on the fly"\n--displaythread (to print thread identity with messages)\n--nostdin to stop listening to stdin\n-p Same as --start-paused\nfile.xml same as --param-file=paramfile\n--debug-connections to print mpiget/send messages\n--user=xxx for user specific options\n'
                sys.exit(0)

        ## ARGUMENT parsing (from the command line):
        for a in arglist:#could be a param file, or ?
            if a[-4:]==".xml":
                self.paramfile+=a.split(",")
                print "Using %s"%self.paramfile
            else:
                try:
                    i=int(a)
                    self.patchno=i
                    print "Batch number: %d"%self.batchno
                except:
                    print "Unrecognised option %s"%a

        ## DEFAULT parameter file name?
        if len(self.paramfile)==0:
            self.paramfile=["params.xml"]
        initDict=None
        if self.initParamString!=None:
            d={"numpy":numpy}
            initDict={}
            exec self.initParamString in d,initDict
        self.config=base.readConfig.AOXml(self.paramfile,batchno=self.batchno,initDict=initDict)
        if self.paramString!=None:
            tmpDict={"this":self.config.this}
            print "Changes to parameter file:\n%s"%self.paramString
            exec self.paramString in tmpDict
            if self.simID=="":
                if self.initParamString!=None:
                    self.simID=self.paramString+";"+self.initParamString
                else:
                    self.simID=self.paramString
        if self.initParamString!=None and self.simID=="":
            self.simID=self.initParamString

        self.config.this.simID=self.simID
        self.sockConn=util.SockConn.SockConn(self.config.getVal("testrunport",default=9000)+self.rank,globals=self.globals,startThread=1,listenSTDIN=self.listenSTDIN)
        os.nice(self.config.getVal("nice",default=nice))
        if self.config.getVal("connectPortDict",default=0,warn=0):
            try:
                t=threading.Thread(target=self.connectPortDict)
                t.daemon=True
                t.start()
                #self.connectPortDict()
            except:
                print "Unable to connect to portdict.py"
        if self.niter==-1:
            try:
                exptime=self.config.getVal("AOExpTime")
                tstep=self.config.getVal("tstep")
                if exptime>0:
                    self.niter=exptime/tstep
                else:
                    self.niter=-1
            except:
                self.niter=-1

####### END of __init__ #############################################################

    def __del__(self):
        
        print "Destroying Ctrl object at iteration %d: self.sockConn.endLoop, cmod.shmem.cleanUp"%self.thisiter
        msg=""
        if self.paramString!=None:
            msg="Commandline parameter changes used for this simulation: %s"%self.paramString
        if self.initParamString!=None:
            msg+=" And %s"%self.initParamString
        if msg!="":
            print msg
        try:
            self.sockConn.endLoop()
        except:
            pass
            #cmod.shmem.cleanUp()
    def connectPortDict(self):
        """This is run as a thread, connects to portdict, and sends info,
        and remains connected until simulation dies..."""
        try:
            file=self.globals["__file__"]
        except:
            file="Unknown"
        data=["add",self.sockConn.host,self.sockConn.port,self.rank,os.environ["USER"],self.simID,file,self.batchno]
        conn=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            conn.connect(("129.234.187.10",8999))
        except:
            print "Couldn't connect to portdict daemon - trying again"
            conn=None
        if conn==None:
            conn=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            try:
                conn.connect(("192.168.3.30",8999))
            except:
                print "Couldn't connect to portdict daemon on gbit network"
                conn=None
        if conn is not None:
            print "Connected to portdict"
            try:
                util.serialise.Send(data,conn)
            except:
                print "Couldn't send info to portdict.py"
        self.portdictsock=conn#keep it in scope!

    def setGlobals(self,globals):
        """Set the globals object.
        @param globals: The globals() dictionary
        @type globals: Dict
        """
        self.globals=globals
    def running(self):
        """Determine whether the simulation is running
        @return: Running flag
        @rtype: Int
        """
        return (self.thisiter<self.niter or self.niter<0) and self.go
    def addnextniters(self,niter):
        """Increase the number of iterations to run before pausing again.
        @param niter:  Number of iterations to increase by
        @type niter: Int
        """
        if self.nextniters==None:
            self.nextniters=niter
        else:
            self.nextniters+=niter
    def getTree(self,regexp=None,base=None,printit=0,ignoreList=[],indent=None):
        """Get the simulation object tree.
        @param regexp: Regular expression to search for
        @type regexp: String
        @param base: Object to start from (in future this may be a regexp)
        @type base: String
        @param printit: Whether to print to terminal
        @type printit: Int
        @param ignoreList: List of objects to be ignored (typically, use parent and config)
        @type ignoreList: List of String
        @param indent: The indentation string
        @type indent: String
        @return: The simulation object tree
        @rtype: String
        """
        saveIndent=self.indent
        if self.globals==None:
            txt="Set globals first"
        else:
            if indent!=None:
                self.indent=indent
            txt=self.parseTree(base,ignoreList=ignoreList)
            if regexp!=None:
                try:
                    regexp=re.compile(regexp)
                except:
                    regexp=None
            if regexp!=None:
                lines=txt.split("\n")
                oklines=[]
                for line in lines:
                    if regexp.match(line.strip())!=None:
                        oklines.append(line)
                txt=string.join(oklines,"\n")
                
            #print "-"*40
            #print txt
        if printit:
            print txt
        self.indent=saveIndent
        return txt
    def parseTree(self,base=None,treeIndent=0,ignoreList=[]):
        """Parse the simulation objects and create a tree.
        @param base: Starting point for parsing
        @type base: String
        @param treeIndent: Number of indents to indent by
        @type treeIndent: Int
        @param ignoreList: List of objects to ignore
        @type ignoreList: List of String
        @return: The tree
        @rtype: String
        """
        txt=""
        istr=self.indent*treeIndent
        if base==None:
            base=self.globals
            #if type(base)==types.DictType:#ie globals
            for key in base.keys():
                if key[0:2]!="__" and key not in ignoreList:
                    #print key,base[key]
                    t=type(base[key])
                    if t==types.InstanceType:
                        txt+=self.parseTree(key,treeIndent,ignoreList=ignoreList)
                    elif t==numpy.ndarray:
                        tmptxt=istr+"%s\t%s (%s, %s)\n"%(key,str(type(base[key])),str(base[key].shape),str(base[key].dtype))
                        txt+=tmptxt
                    else:
                        tmptxt=istr+"%s\t%s\n"%(key,str(type(base[key])))
                        txt+=tmptxt
                        #print tmptxt
        else:#a recursive call... key will be a string
            mydict={}
            execStr="t=type(%s) ; d=dir(%s)"%(base,base)
            #print "Executing %s at indent level %d"%(execStr,treeIndent)
            try:
                exec  execStr in self.globals,mydict
            except:
                tmptxt=istr+"%s\tUNKNOWN\n"%base
                txt+=tmptxt
                #print tmptxt
            else:
                t=mydict["t"]
                if t==types.InstanceType:#a class instance-ie object.
                    tmptxt=istr+"%s\t%s\n"%(base,str(t))
                    txt+=tmptxt
                    #print tmptxt
                    
                    for dirstr in mydict["d"]:
                        if dirstr[0:2]!="__" and dirstr not in ignoreList:#evaluate the object.
                            #print "Evaluating %s . %s"%(base,dirstr)
                            txt+=self.parseTree("%s.%s"%(base,dirstr),treeIndent+1,ignoreList=ignoreList)
                elif t==numpy.ndarray:
                    exec "s=%s.shape ; tc=%s.dtype"%(base,base) in self.globals,mydict
                    tmptxt=istr+"%s\t%s (%s, %s)\n"%(base,str(t),str(mydict["s"]),str(mydict["tc"]))
                    txt+=tmptxt
                elif t in [types.IntType,types.LongType,types.FloatType,types.BooleanType,types.ComplexType,types.NoneType]:
                    exec "v=str(%s)"%base in self.globals,mydict
                    tmptxt=istr+"%s\t%s %s\n"%(base,str(mydict["t"]),mydict["v"])
                    txt+=tmptxt
                    #print tmptxt
                elif t in [types.StringType,types.ListType,types.TupleType,types.DictType]:
                    exec "l=str(len(%s))"%base in self.globals,mydict
                    tmptxt=istr+"%s\t%s length: %s\n"%(base,str(mydict["t"]),mydict["l"])
                    txt+=tmptxt
                else:
                    tmptxt=istr+"%s\t%s\n"%(base,str(mydict["t"]))
                    txt+=tmptxt
                    #print tmptxt
        return txt
    #def end(self):
    #    """Force an end to the simulation"""
    #    tlist=threading.enumerate()
    #    print "Threads left: %s"%str(tlist)
    #    for t in tlist:
    #        if t!=threading.currentThread():
    #            t._Thread__stop()
    #    cmod.shmem.cleanUp()
    def initialCommand(self,cmd,action=None,ret="",freq=1,startiter=0,repeat=1):
        """Add commands to be computed initially when the simulation runs.
        This is useful e.g. for creating a poke matrix automatically etc.
        This should be called before ctrl.mainloop is called.
        @param cmd: The python command string to be executed
        @type cmd: String
        @param action: The action to take, ie when to execute the command.
        @type action: None or String
        @param ret: The return value, usually a zero length string
        @type ret: String
        @param freq: The frequency with which the command is computed.  If this is -1, the command is executed only once.
        @type freq: Int
        @param startiter: The first iteration number after which the command is executed.
        @type startiter: Int
        @param repeat: Whether the command is appended to the repeating list, or to the single execution list
        @type repeat: Int
        """
        data=util.ConnObj.ConnMsg(cmd,action,ret=ret)
        if repeat:
            self.sockConn.rptCmdList.append([data,None,freq,startiter])
        else:
            self.sockConn.cmdList.append([data,None])
    
    def mainloop(self,compList,sockConn=None,cleanShmem=1):
        """Run the simulation main loop.
        @param compList:  The list of simulation objects who's next method should be called each iteration.
        @type compList: List
        @param sockConn: SockConn object for communication
        @type sockConn: SockConn instance
        @param cleanShmem: Flag, clean up shared memory at end
        @type cleanShmem: Int
        """
        self.compList=compList
        self.compListNames=[]
        for module in compList:
            self.compListNames.append(module.objID)
            module.finalInitialisation()
        self.simStartTime=time.time()
        print "Took %g seconds to initialise"%(self.simStartTime-self.simInitTime)
        self.thisIterTiming=numpy.zeros((len(compList),),numpy.float64)
        self.meanTiming=numpy.zeros((len(compList),),numpy.float64)
        self.meanTiming2=numpy.zeros((len(compList),),numpy.float64)
        self.meanClock=numpy.zeros((len(compList),),numpy.float64)
        rangeLenCompList=range(len(self.compList))
        self.simctrlXML=self.createQueryObjs()
        self.simctrlXMLRestricted=self.createQueryObjs(addDir=0)
        print "Waiting at MPI barrier..."
        self.mpiComm.barrier()#wait til they're all ready to start - not essential, but nice...
        if sockConn==None:
            sockConn=self.sockConn
        #and now share the connection settings...
        print "Sharing connection setting"
        connParamsDict={}
        for i in range(self.mpiComm.size):
            connParamsDict[i]=numpy.zeros((5,),numpy.int32)
            if i==self.mpiComm.rank:
                connParamsDict[i][:4]=map(int,socket.gethostbyname(sockConn.host).split("."))
                connParamsDict[i][4]=sockConn.port
            self.mpiComm.broadcast(connParamsDict[i],i)
            # and put them into a dict...
            connParamsDict[i]=(string.join(map(str,connParamsDict[i][:4]),"."),int(connParamsDict[i][4]))
        #and store in config so that they're accessible to all...
        self.config.connectionParamsDict=connParamsDict
        
        mem=self.checkMem()
        print "Using approximately %gGB"%(mem/1024./1024./1024.)
        # if mem>4*1024*1024*1024:
        #     print "Warning: Using more than 4GB memory - inefficient on the cray (%gGB)"%(mem/1024./1024/1024)
        pausedMsg=0
        for cmd in self.initCmdList:
            exec cmd in self.globals,{}
        try:
            sockConn.doCmdList(self.thisiter)
        except:
            print "\n\n\n\nError in sockConn.doCmdList\n\n\n\n"
        print "Entering main loop"
        while self.running():
            if not self.paused:
                if self.debug!=None:
                    print "Ctrl: rank %d (debug=%s) doing iteration %d"%(self.rank,str(self.debug),self.thisiter)
                pausedMsg=0
                t=time.time()
                for i in rangeLenCompList:
                    self.compListPos=i
                    module=compList[i]
                    if self.running():
                        if self.debug!=None:
                            print "Ctrl: starting %s"%module.objID
                        t1=time.time()
                        c1=time.clock()
                        try:
                            module.doNextIter()#generateNext()
                        except:
                            print "Error in generate next for %dth module (iter %d)"%(self.compListPos,self.thisiter)
                            raise
                        c2=time.clock()-c1#get the CPU time... (resolution typically 0.01s).
                        t2=time.time()-t1
                        self.thisIterTiming[i]=t2
                        self.meanTiming[i]+=t2
                        self.meanTiming2[i]+=t2*t2
                        self.meanClock[i]+=c2
                        if self.debug!=None:
                            print "Ctrl: Time taken for %s was %g seconds"%(module.objID,t2)
                    else:
                        break
                #data=parent.next("stickman")
                self.frametime=time.time()-t
                self.thisiter+=1
                #print "Done %d iterations"%self.thisiter
                if self.nextniters!=None:
                    self.nextniters-=1
                    if self.nextniters==0:
                        self.nextniters=None
                        self.paused=1
                if self.slowMotion!=0:
                    time.sleep(self.slowMotion)
            else:
                if self.debug!=None:
                    print "Ctrl: rank %d (debug=%s) paused at iteration %d"%(self.rank,str(self.debug),self.thisiter)
                #print "Paused..."
                if pausedMsg==0:
                    print "Paused at iteration %d"%self.thisiter
                    pausedMsg=1
                #time.sleep(1)
                select.select(sockConn.selIn,[],sockConn.selIn,1.)
                time.sleep(0.001)#allow other thread to operate.
            try:
                sockConn.doCmdList(self.thisiter)
            except:
                print "\n\n\n\nError in sockConn.doCmdList\n\n\n\n"
##         tlist=threading.enumerate()
##         print "Threads left: %s"%str(tlist)
##         for t in tlist:
##             if t.getName()!="MainThread":
##                 t._Thread__stop()
        print "Finished - calling endSim for each module"
        #Sum = sum(self.meanTiming) # UB: to sum up the time spent at modules
        for i in rangeLenCompList:
            self.compListPos=i
            module=compList[i]
           # print "%s:  %.4f s, %.2f %%"%(module.objID, self.meanTiming[i],
           #                               self.meanTiming[i]/Sum*100) # UB 2012Jul23
            module.endSim()
        #print "Sum over modules: {0} s".format(Sum) # UB
        print "waiting at mpi barrier for all other processes to finish"
        self.mpiComm.barrier()#wait til they're all ready to finish - not essential, but nice...
        if cleanShmem:
            cmod.shmem.cleanUp()
        t=time.time()
        print "Total time %gs, running time %gs"%(t-self.simInitTime,t-self.simStartTime)
        time.sleep(1)#allow a bit of time before abort is called - to allow all semaphores to be cleaned up.

    def createQueryObjs(self,addHeader=1,addDir=1):
        """Create XML for all the science objects, such that it can be used
        by a GUI for querying..."""
        if addHeader:
            s='<?xml version="1.0"?>\n<simdata>\n'
        else:
            s=''
        gkey=self.globals.keys()
        indx=-1
        for sciobj in self.compList:
            indx+=1
            objname=""
            for key in gkey:
                if self.globals[key]==sciobj:
                    objname=key
                    break
            if objname=="":#index using the list...
                objname="ctrl.compList[%d]"%indx
            s+='<simobj name="%s">\n'%sciobj.objID#objname
            #print objname
            if hasattr(sciobj,"plottable"):
                s+=sciobj.plottable(objname)
            if addDir:#add a listing of all elements of the object...
                for a in dir(sciobj):
                    if a[:2]!="__":
                        obj=getattr(sciobj,a)
                        if type(obj) not in [types.MethodType,types.InstanceType]:
                            s+='<plot title="%s.%s" cmd="data=%s.%s" ret="data" texttype="1" wintype="mainwindow" when="cmd"/>\n'%(objname,a,objname,a)
                        if (type(obj)==numpy.ndarray and obj.dtype not in [numpy.complex64,numpy.complex128]):
                            s+='<plot title="Plot: %s.%s" cmd="data=%s.%s" ret="data" type="gist" window="1" palette="gray.gp" when="rpt"/>\n'%(objname,a,objname,a)
            s+='</simobj>\n'
        s+='<plot title="Poke all" ret="data" when="cmd" texttype="1" wintype="mainwindow">\n<cmd>\nfor obj in ctrl.compList:\n if obj.control.has_key("poke"):\n  obj.control["poke"]=1\n if obj.control.has_key("zero_dm"):\n  obj.control["zero_dm"]=1\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=1\n print obj.control\nprint "Starting poking"\ndata="Starting poking"</cmd>\n</plot>\n'
        s+='<plot title="Get ref cents" ret="data" when="cmd" texttype="1" wintype="mainwindow">\n<cmd>\nfor obj in ctrl.compList:\n if obj.control.has_key("takeRef"):\n  obj.control["takeRef"]=1\n if obj.control.has_key("zero_dm"):\n  obj.control["zero_dm"]=1\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=1\n print obj.control\nprint "Starting reference centroids"\ndata="Starting reference centroids"</cmd>\n</plot>\n'
        s+='<plot title="Close loop, zero science" ret="data" when="cmd" texttype="1" wintype="mainwindow">\n<cmd>\nfor obj in ctrl.compList:\n if obj.control.has_key("zero_science"):\n  obj.control["zero_science"]=1\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=0\n if obj.control.has_key("zero_dm"):\n  obj.control["zero_dm"]=1\n if obj.control.has_key("close_dm"):\n  obj.control["close_dm"]=1\n print obj.control\ndata="Zeroing science, closing loop etc"\nprint data</cmd>\n</plot>\n'
        s+='<plot title="Open loop, zero science" ret="data" when="cmd" texttype="1" wintype="mainwindow">\n<cmd>\nfor obj in ctrl.compList:\n if obj.control.has_key("zero_science"):\n  obj.control["zero_science"]=1\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=0\n if obj.control.has_key("zero_dm"):\n  obj.control["zero_dm"]=1\n if obj.control.has_key("close_dm"):\n  obj.control["close_dm"]=0\n print obj.control\ndata="Zeroing science, opening loop etc"\nprint data</cmd>\n</plot>\n'
        s+='<plot title="Reset science" ret="data" when="cmd" texttype="1" wintype="mainwindow">\n<cmd>\nfor obj in ctrl.compList:\n if obj.control.has_key("zero_science"):\n  obj.control["zero_science"]=1\n print obj.control\ndata="Zeroing science"\nprint data</cmd>\n</plot>\n'
        s+="""<plot title="Control status" ret="data" when="cmd" texttype="1" wintype="mainwindow">\n<cmd>\ndata=""\nfor obj in ctrl.compList:\n print '%s: %s'%(obj.objID,str(obj.control))\n data+=obj.objID+": "+str(obj.control)+"\\n"\n</cmd>\n</plot>\n"""
        s+='<plot title="Timings" ret="data" when="rpt" texttype="1" wintype="ownwindow" textreplace="1">\n<cmd>\ndata=(ctrl.compListNames,ctrl.thisIterTiming,ctrl.meanTiming,ctrl.meanTiming2,ctrl.thisiter,ctrl.frametime,ctrl.meanClock)\n</cmd>\nimport numpy\nres="\tmean\tstdev\tmean cpu clocks\tThis iter\\n"\nfor i in range(len(data[0])):\n  res+="%s:\t%g\t%g\t%g\t%g\\n"%(data[0][i],data[2][i]/data[4],numpy.sqrt(data[3][i]/data[4]-(data[2][i]/data[4])**2),data[6][i]/data[4],data[1][i])\ndata=res+"Frame %d (%g fps, %g spf)"%(data[4],1./data[5],data[5])\n</plot>\n<plot title="Iteration counter" cmd="data=(ctrl.thisiter,ctrl.frametime)" ret="data" when="rpt" texttype="1" wintype="ownwindow" textreplace="1">\ndata="Frame %d (%g fps, %g spf)"%(data[0],1.0/data[1],data[1])\n</plot>\n'
        if addHeader:
            s+="</simdata>\n"
        #self.simctrlXML=s#save the XML for the simctrl gui.
        return s
    def checkMem(self):
        try:
            mem=int(open("/proc/%d/stat"%os.getpid()).read().split()[22])
            #mem2=open("/proc/%d/statm"%os.getpid()).read().split()[0]*4096
        except:
            mem=-1
        return mem

    def doInitialPoke(self,startiter=0,reconNumber=None):
        """This is a function that can be called during simulation setup to make the simulation perform a poke at the start of iteration zero.
        """
        if reconNumber==None:
            cmd="""for obj in ctrl.compList:\n if obj.control.has_key("poke"):\n  obj.control["poke"]=1\n if obj.control.has_key("science_integrate"):\n  obj.control["science_integrate"]=0\n if obj.control.has_key("zero_dm"):\n  obj.control["zero_dm"]=1\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=1\n print obj.control\nprint "Starting poking"\n"""
        elif type(reconNumber)==type(0):
            cmd="""for obj in ctrl.compList:\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=1\n if obj.control.has_key("science_integrate"):\n  obj.control["science_integrate"]=0\nreconList[%s].control["poke"]=1\nreconList[%s].control["zero_dm"]=1\nprint "Starting poking DM %s"\n"""%(reconNumber,reconNumber,reconNumber)
        elif type(reconNumber)==type(""):
            #its an idstr...
            cmd="""for obj in ctrl.compList:\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=1\n if obj.control.has_key("science_integrate"):\n  obj.control["science_integrate"]=0\nfor r in reconList:\n if r.idstr[0]=='%s':\n  r.control["poke"]=1\n  r.control["zero_dm"]=1\nprint "Starting poking DM %s"\n"""%(reconNumber,reconNumber)
        self.initialCommand(cmd,freq=-1,startiter=startiter)
    def doInitialReferenceCentroids(self,startiter=0):
        """This function can be called during simulation setup to make the simulation take reference centroids.  Useful for cases with LGS spot elongation.
        """
        cmd="""for obj in ctrl.compList:\n if obj.control.has_key("takeRef"):\n  obj.control["takeRef"]=1\n if obj.control.has_key("zero_dm"):\n  obj.control["zero_dm"]=1\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=1\nprint "Taking reference centroids"\n"""
        self.initialCommand(cmd,freq=-1,startiter=startiter)
        
    def doInitialPokeThenRun(self,startiter=0):
        """This function can be called during sim setup to make the simulation perform a poke from iteration zero, followed by computation of the reconstructor, and then run using this.
        """
        cmd="""for obj in ctrl.compList:\n if obj.control.has_key("poke"):\n  obj.control["poke"]=1\n if obj.control.has_key("science_integrate"):\n  obj.control["science_integrate"]=0\n if obj.control.has_key("zero_dm"):\n  obj.control["zero_dm"]=1\n if obj.control.has_key("cal_source"):\n  obj.control["cal_source"]=1\n print obj.control\nprint "Starting poking"\n"""
        self.initialCommand(cmd,freq=-1,startiter=startiter)

        cmd="""maxmodes=0\nfor obj in ctrl.compList:\n if hasattr(obj,"nmodes") and obj.nmodes>maxmodes:\n  maxmodes=obj.nmodes\nctrl.initialCommand("ctrl.doSciRun()",freq=-1,startiter=maxmodes+10)\nprint "Will close loop after %d iterations"%(maxmodes+10)\n"""
        self.initCmdList.append(cmd)
    def doInitialSciRun(self,startiter=0):
        cmd="ctrl.doSciRun()"
        self.initialCommand(cmd,freq=-1,startiter=startiter)
        
    def doSciRun(self):
        """Can be called during sim setup, or by GUI to zero science/dms etc and close the loop.
        Used by doInitialPokeThenRun().
        """
        for obj in self.compList:
            if obj.control.has_key("zero_dm"):
                obj.control["zero_dm"]=1
            if obj.control.has_key("zero_science"):
                obj.control["zero_science"]=10
            if obj.control.has_key("science_integrate"):
                obj.control["science_integrate"]=1
            if obj.control.has_key("cal_source"):
                obj.control["cal_source"]=0
            if obj.control.has_key("close_dm"):
                obj.control["close_dm"]=1
            print obj.control
        print "Zeroing science,closing loop etc"
        
        
class myStdout:
    """A class for customising print messages"""
    def __init__(self,rank):
        self.rank=rank
        self.displayRank=1
        self.displayThread=0
    def write(self,txt):
        if len(txt.strip())>0:
            if self.displayRank:
                r="%d"%self.rank
            else:
                r=""
            if self.displayThread:
                t=" (%012d)"%thread.get_ident()
            else:
                t=""
            sys.__stdout__.write(">>>%s%s: %s\n"%(r,t,txt))
            sys.__stdout__.flush()
    def flush(self):
        sys.__stdout__.flush()
