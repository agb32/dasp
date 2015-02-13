#$Id: readConfig.py,v 1.45 2011/01/24 08:00:40 ali Exp $
import xml.parsers.expat,types,string
import time,os,sys,numpy
import socket,util.serialise
# try:#this is used for sharing live data between modules/processes.
#     from Scientific.MPI import world as MPIWorld
# except:
#     class dummy:
#         rank=0
#         size=1
#         def abort(self,rt):
#             sys.exit(0)
#     MPIWorld=dummy()

class AOXml:
    """A class for reading an AO XML config file.  If given, file (filename)
    is parsed on init.  If given, searchOrder defines the order in which
    objects are searched when looking for a variable.  searchOrder should be
    a list of module names, first will be searched first.  Additionally, a
    batch number can be specified, and then only variables that are in
    agreement with this are returned.
    Class variables (important to simulation programmer):
     - writeSchema - int flag, allows schema to be written to XML file.  Not yet implemented
     - filename - string, the file to parse
     - batchno - int, the batch number to use when parsing the XML file.
    Class variables (not important to simulation programmer):
     - p - object, the XML parser
     - this - object, the internal structure of the XML file
     - fileData - object, used for the parameter GUI, holding information for display
     - whichElement - int, element to be used for current batch
     - inSchema - int flag, whether currently in schema part of XML file
     - inside - list, of tags which tags currently inside parsing
     - batchInfo - list, of batch numbers of parent tags
     - curbatch - list, not used
     - invar - int flag, whether currently in a var tag
     - curmodule - string, name of module currently being evaluated
     - storedattrs - dictionary, of attributes from XML tag
     - storedtxt - string, of text from the XML value tag
     - searchOrder - list, the order of modules to search when getting a variable
     - batchReset - list, not used
     - ignoreModule - int flag, whether a module should be ignored (batch number) while parsing
    @cvar writeSchema: If set, allows schema to be written to XML file.  Not yet implemented
    @type writeSchema: int flag
    @cvar filename: Filename of XML file to be parsed
    @type filename: string
    @cvar batchno: Batch number to search for when parsing
    @type batchno: int
    @cvar p: XML parser object
    @type p: XML parser object
    @cvar this: Internal structure holding XML file info
    @type this: Object
    @cvar fileData: Structure holding XML info
    @type fileData: Object
    @cvar whichElement: element used for current batch
    @type whichElement: int
    @cvar inSchema: If set, currently in schema part of XML file
    @type inSchema: Int
    @cvar inside: List of tags containing tags currently inside parsing
    @type inside: List
    @cvar batchInfo: list of batch numbers of parent tags
    @type batchInfo: List
    @cvar curbatch: not used
    @type curbatch: List
    @cvar invar: flag, whether evaluating a var tag
    @type invar: Int
    @cvar curmodule: Name of module currently being evaluated
    @type curmodule: String
    @cvar storedattrs: Attributes from current XML tag
    @type storedattrs: Dict
    @cvar storedtxt: Text from current XML tag
    @type storedtxt: String
    @cvar searchOrder: Order to search modules for a variable
    @type searchOrder: List
    @cvar batchReset: Not used
    @type batchReset: List
    @cvar ignoreModule: Flag, whether module should be ignored (batch number mismatch)
    @type ignoreModule: Int
    """ 
    def __init__(self,file=None,batchno=0,writeSchema=0,ignoreError=0,initDict=None,mpiRankSizeAbort=(0,1,None)):
        """Initialise a AO XML parser object.
        @param file: Filename
        @type  file: String
        @param writeSchema: Whether to write schema when saving XML file
        @type  writeSchema: Int flag
        @param batchno: Batch number that we're interested in
        @type  batchno: Int
        mpiRankSizeAbort is the mpi rank, worldsize, and an abort function.
        """
        self.p=None
        self.writeSchema=writeSchema
        if type(file)!=type([]):
            file=[file]
        self.fileList=file
        self.filename=None
        self.batchno=batchno
        self.initDict=initDict
        self.mpiRankSizeAbort=mpiRankSizeAbort
        self.rank=self.mpiRankSizeAbort[0]
        self.reset()
        for f in file:
            self.p = xml.parsers.expat.ParserCreate()
            self.p.StartElementHandler = self.start_element
            self.p.EndElementHandler = self.end_element
            self.p.CharacterDataHandler = self.char_data
            self.p.returns_unicode=0

            self.open(f,ignoreError)
    def abort(self,rt=0):
        rank,size,abt=self.mpiRankSizeAbort
        if abt!=None:
            if size==1:
                abt()
            else:
                abt(rt)

    def open(self,file,ignoreError=0):
        """Open an XML file and parse it.
        @param file: Filename
        @type  file: String
        """
        if file!=None:
            print "\n%s\n\nReading param file %s\n\n%s\n"%("*"*70,file,"*"*70)
            setattr(self.this,"filename",file)
            if file[-4:]==".xml":
                #txt=open(file).read()
                txt=PreFormatXML(file).txt
                self.filename=file
                self.error=0
                self.parse(txt)
            elif file[-3:]==".py":
                print "USING PYTHON CONFIG FILE %s"%file
                txt=open(file).read()
                self.filename=file
                self.error=0
                self.loadPyConfig(txt)
            else:
                raise Exception("Unknown file type %s"%file)
            if self.error and ignoreError==0:
                self.reset()
            else:
                if hasattr(self.this,"globals"):
                    if not hasattr(self.this.globals,"ncpu"):
                        setattr(self.this.globals,"ncpu",self.this.ncpu)

    def loadPyConfig(self,txt):
        """Load a python configuration file.
        For more info, see example in test/scao/params.py
        """
        self.this.globals=This()
        d={"new":This,"this":self.this,"batchno":self.batchno,"batchNumber":self.batchno,"numpy":numpy,"filename":self.filename,"ncpu":self.getCpus()}
        try:
            exec txt in d
        except:
            self.error=1
            raise
        e={}
        exec "" in e
        for key in e.keys():
            if d.has_key(key):
                del(d[key])
        this=d["this"]
        for key in d.keys():
            if key not in ["this","new","batchNumber","batchno","filename","ncpu"]:
                setattr(self.this.globals,key,d[key])
        so=dir(self.this)
        rem=["__doc__","__init__","__module__","ncpu","numpy",'batchNumber', 'batchno', 'filename',"globals"]
        for r in rem:
            try:
                so.remove(r)
            except:
                pass
        so.sort(reverse=True)
        so.append("globals")
        self.searchOrder=so#["globals"]
    def reset(self):
        """Reset the object ready to parse a different file"""
        self.error=0
        self.p = xml.parsers.expat.ParserCreate()
        self.p.StartElementHandler = self.start_element
        self.p.EndElementHandler = self.end_element
        self.p.CharacterDataHandler = self.char_data
        self.p.returns_unicode=0
        self.this=This()
        self.used=This()#all the used variables...
        self.fileData=FileData()
        setattr(self.this,"batchNumber",self.batchno)
        setattr(self.this,"batchno",self.batchno)
        setattr(self.this,"filename",self.filename)#parameter file name
        setattr(self.this,"ncpu",self.getCpus())
        if self.initDict!=None:
            setattr(self.this,"globals",This())
            for key in self.initDict.keys():
                setattr(self.this.globals,key,self.initDict[key])
        #self.thisBatch=This()
        self.whichElement=0#used if not current batch...
        self.inSchema=0
        self.inside=[]
        self.batchInfo=[]
        self.curbatch=[]
        self.invar=0
        self.curmodule=None
        self.storedattrs=None
        self.storedtxt=None
        self.searchOrder=[]
        self.batchReset=[]
        self.ignoreModule=0
        self.postList=[]#list of data that have been posted by simulation processes... (live sharing of data)
        self.connectionParamsDict={}#dictionary of rank:(host,port)
        self.rank=self.mpiRankSizeAbort[0]#MPIWorld.rank
        self.rankSize=self.mpiRankSizeAbort[1]#MPIWorld.size
    def start_element(self,name, attrs):
        """Called by XML parser at opening of tag.
        @param name: Name of tag
        @type  name: String
        @param attrs: Attributes set in tag
        @type attrs: Dict
        """
        #print 'Start element:', name, attrs
        self.inside.append(name)
        if self.batchno==None:
            self.this.batchNumber=0
        if self.curmodule!=None and self.ignoreModule==1:#already in a module
            self.batchInfo.append(-1)
        elif self.invar==1 and self.ignoreVariable==1:
            self.batchInfo.append(-1)
        elif attrs.has_key("batchno"):
            tryexec=0
            try:
                b=eval(attrs["batchno"])
            except:
                tryexec=1
            if tryexec:
                d={}
                attrs["batchno"]=attrs["batchno"].replace("\\n","\n")
                exec(attrs["batchno"]) in d
                b=d["batchno"]
            if type(b)==types.IntType:
                b=[b]
            if type(b)!=types.ListType:
                print "ERROR:  batchno value must be evaluatable to a list"
                raise Exception("ERROR:  batchno value must be evaluatable to a list")
            if self.batchno==None or self.batchno in b:
                self.batchInfo.append(1)#use for this batch
            else:
                self.batchInfo.append(-1)#ignore for this batch
            if self.batchno==None:
                self.this.batchNumber=b[0]
        else:
            if self.curmodule!=None and self.batchModule==1:
                self.batchInfo.append(1)#only this batch
            elif self.invar==1 and self.batchVariable==1:
                self.batchInfo.append(1)#only this batch
            else:
                self.batchInfo.append(0)#use for all batches...
        self.storedattrs=None
        self.storedtxt=None
        if name=="module":
            if self.inSchema==0:
                self.ignoreModule=0
                self.ignoreVariable=0
                self.ignoreVar=0
                self.batchModule=0
                self.batchVariable=0
                self.batchVar=0
                self.fileData.modules.append(ModuleClass(attrs["name"],attrs))
                if self.batchInfo[-1]==-1:
                    self.ignoreModule=1
                else:
                    if hasattr(self.this,attrs["name"]):
                        if self.batchInfo[-1]==0:#not specific to this batch and already defined...
                            print "WARNING- redefining module %s.  Continuing anyway..."%attrs["name"]
                    else:
                        setattr(self.this,attrs["name"],This())
                self.searchOrder.insert(0,attrs["name"])
                self.curmodule=attrs["name"]
            else:
                if self.writeSchema:
                    self.prepareSchema(name,attrs)
        elif name=="variables":
            if self.curmodule==None:
                print "ERROR: <variables> must be within <module>"
                self.invar=0
                raise "ERROR: <variables> must be within <module>"
            else:
                self.invar=1
                self.ignoreVariable=0
                self.ignoreVar=0
                self.batchVariable=0
                self.batchVar=0
                if self.ignoreModule==1 or self.batchInfo[-1]==-1:
                    self.ignoreVariable=1
                elif self.batchInfo[-1]==1:
                    self.batchVariable=1
                
        elif name=="var":
            if self.invar==0:
                print "ERROR - <var> must be within <variables>"
                raise "ERROR - <var> must be within <variables>"
            else:
                self.storedattrs=attrs
                self.ignoreVar=0
                self.batchVar=0
                if self.ignoreVariable==1 or self.batchInfo[-1]==-1:
                    self.ignoreVar=1
                elif self.batchInfo[-1]==1:
                    self.batchVar=1
        elif name=="schema":
            self.inSchema=1
            if self.writeSchema:
                if attrs.has_key("filename"):
                    self.schemaFile=open(attrs["filename"],"w")
                else:
                    self.schemaFile=dudFile()
        elif name=="include":
            p = xml.parsers.expat.ParserCreate()
            p.StartElementHandler = self.start_element
            p.EndElementHandler = self.end_element
            p.CharacterDataHandler = self.char_data
            p.returns_unicode=0
            if attrs.has_key("file"):
                file=attrs["file"]
                if file[0]=="'":
                    file=eval(file)
                txt=open(file,"r").read()
                p.Parse(txt)
    def end_element(self,name):
        """Called by XML parser at end of element
        @param name: Name of tag to be closed
        @type name: String
        """
        #print 'End element:', name
        n=self.inside.pop()
        bInfo=self.batchInfo.pop()
        if n=="variables":
            self.invar=0
        elif n=="module":
            self.curmodule=None
        elif n=="var":#ending a var, so now create it...
            attrs=self.storedattrs
            if self.storedtxt==None:
                txt=None
            else:
                txt=self.storedtxt.strip()
            if not attrs.has_key("value"):
                if txt[0]=='\n':
                    txt=txt[1:]
                if txt[-1]=='\n':
                    txt=txt[:-1]
                attrs["value"]=txt#which may be none - which is an error.
                attrs["extendedTxt"]=1
            if not attrs.has_key("type"):#unknown type... eval it...
                attrs["type"]="eval"
            if self.ignoreVar==0 and attrs.has_key("name") and len(attrs["name"])>0:
                #v=self.makeVal(attrs)
                try:
                    v=self.makeVal(attrs)
                except:
                    print "ERROR - in XML file, while making value."
                    print "This occured for element %s near line %d of file %s."%(attrs["name"],self.p.CurrentLineNumber,self.filename)
                    #print sys.exc_info()
                    #import traceback
                    #traceback.print_tb(sys.exc_info()[2])
                    print sys.exc_type
                    print sys.exc_value
                    v=None
                    self.error=1
                    raise
                if hasattr(self.this,self.curmodule):
                    if hasattr(getattr(self.this,self.curmodule),attrs["name"]):
                        if bInfo==1:#overwrite prev value
                            if type(v)==numpy.ndarray:
                                strv="Array..."
                            else:
                                strv=str(v)
                            print "Overwriting previous value for variable %s (new value: %s)"%(str(attrs["name"]),strv)
                            setattr(getattr(self.this,self.curmodule),attrs["name"],v)
                        else:
                            print "WARNING - value %s defined previously.  Will not overwrite"%attrs["name"]
                    else:
                        setattr(getattr(self.this,self.curmodule),attrs["name"],v)
                else:
                    print "ERROR: doesn't have module",self.curmodule
##                 for module in self.fileData.modules:
##                     if module.name==self.curmodule:
##                         module.variableList.append(attrs)
##                         break
                self.fileData.modules[-1].variableList.append((attrs,v))    
        elif n=="schema":
            self.inSchema=0
            if self.writeSchema:
                self.schemaFile.close()
            
        if n!=name:
            print "ERROR - exiting element we weren't inside",n,name
            raise "XML ERROR - not inside element"
        if self.batchno==None:
            self.this.batchNumber=None

    def char_data(self,data):
        """Called while parsing a multi-line tag.  Stores the data internally.
        @param data: The data to be stored
        @type data: String
        """
        if self.storedtxt==None:
            self.storedtxt=data
        else:
            self.storedtxt+=data
        #print 'Character data:', data,len(data)#repr(data)
            
    def makeVal(self,attrs):
        """Calculate a value from an XML tag.
        @param attrs: XML tag attributes
        @type attrs: Dict
        @return: The value
        @rtype: XML file defined
        """
        #put current module contents in local space...
        if "_" in self.curmodule:
            mlist=["globals",reduce(lambda x,y:x+"_"+y,self.curmodule.split("_")[:-1]),self.curmodule]
        else:
            mlist=["globals",self.curmodule]
        glob={}
        for mod in mlist:
            #put all the variables into the locals dictionary...
            if hasattr(self.this,mod):
                obj=getattr(self.this,mod)
                vars=dir(obj)
                for var in vars:
                    if var[:2]!="__":
                        glob[var]=getattr(obj,var)
                glob[mod]=obj
        glob["this"]=self.this
        glob["numpy"]=numpy
        glob["module"]=self.curmodule
        #if len(attrs["value"])==0:
        #    raise KeyError
        v=attrs["value"]#okay to raise error if no key...
        if attrs.has_key("type"):
            t=attrs["type"]
            if t=="i":#eg value will be "5"
                v=int(v)
            elif t=="f":#eg value will be "5.9"
                v=float(v)
            elif t=="numpy":#eg value will be "numpy.array([1,2,3])"
                v=eval(v,glob)
            elif t=="copy":#eg value will be "this.globals.npup"
                v=eval(v,glob)
            elif t=="list":#eg value will be "[1,2,3]"
                v=eval(v,glob)
            elif t=="dict":#eg value will be "{1:2, 3:4}"
                v=eval(v,glob)
            elif t=="code":#eg if attrs["name"]="hi" value could be "hi=5.0**2"
                v=self.doexec(v,attrs["name"],glob)
            elif t=="eval":#eg could be anything that can be eval'ed...
                v=eval(v,glob)
            elif t=="string" or t=="s":
                pass
            else:
                print("ERROR:readConfig: Unrecognised type for: '**"+   
                     str(attrs["name"].strip())+"**'(='**"+str(t)+"**'); "+
                     "Assuming string")
        return v

    def doexec(self,strng,name,glob):
        """Function to exec a string safely with the "this" dictionary,
        returning the value obtained from the exec.
        @param strng: The string to be exec'd
        @type strng: String
        @param name: The value to be returned
        @type name: String
        @return: The value calculated by execing strng.
        @rtype: User defined.
        """
        dct=glob#{"this":self.this}
        if len(strng.strip())==0:
            print "ERROR - Zero length strng in doexec...",name
            raise Exception("ZERO LENGTH STRNG FOR EXEC")
        try:
            exec strng in dct
        except:
            print "ERROR - executing XML python code:"
            print strng
            print sys.exc_type
            print sys.exc_value
            raise
        return dct[name]

    def parse(self,txt,isfinal=1):
        """A soft wrapper for p.Parse()
        @param txt: The XML text to be parsed
        @type txt: String
        @param isfinal: Not used
        @type isfinal: Int
        """
        self.p.Parse(txt)
    def __repr__(self):
        exclude=["__doc__","__init__","__module__","__repr__","batchNumber","batchno"]
        s="Parameters for batch number %d\n"%self.this.batchNumber
        thisobj=dir(self.this)
        
        for e in exclude:
            if e in thisobj:
                thisobj.remove(e)
        for sobj in thisobj:
            s+=sobj+"\n"
            obj=getattr(self.this,sobj)
            a=dir(obj)
            for e in exclude:
                if e in a:
                    a.remove(e)
            for sobj2 in a:
                obj2=getattr(obj,sobj2)
                s+="  "+str(sobj2)+"="+str(obj2)+"\n"
        s+="\nMODULE TEXT\n-----------\n"
        for module in self.fileData.modules:
            s+="MODULE "+module.name+":"+str(module.args)+"\n"
            for var in module.variableList:
                s+=str(var[0])+"\n"
            s+="\n"
        return s

    def setSearchOrder(self,searchOrder):
        """The default search order is in reverse order from the config file
        layout.  This probably isn't usually desired (eg an atmos module
        doesn't want to pick up a variable means for a reconstructor module).
        @param searchOrder: A list of strings defining the search order when
        looking for a value
        @type searchOrder: List
        """
        self.searchOrder=searchOrder
    def getVal(self,varname,default=None,searchOrder=None,raiseerror=1,warn=1):
        """Return the value of a variable stored in the XML file.
        This would be used e.g. myvar=x.getVal("myvar") where x is the
        instance of AOXml (self).  A default value can be given.  If no default
        is given and nothing is found, an error is raised if raiseerror is set.
        If the serachOrder isn't specified, the default is used.  Otherwise,
        the searchOrder must be a list of strings.
        @param varname: The variable to obtain
        @type varname: String
        @param default: The value to return if variable isn't found (None means no default)
        @type default: User defined
        @param searchOrder: The order in which to search modules
        @type searchOrder: List of strings
        @param raiseerror: Flag determining whether an error should be raised if the variable isn't found and if the default is None.
        @type raiseerror: Int
        @return:  The value of variable in the XML file
        @rtype: User defined
        """
        if searchOrder==None:
            searchOrder=self.searchOrder
        found=0
        val=default
        for module in searchOrder:
            if hasattr(self.this,module):
                tmod=getattr(self.this,module)
                if hasattr(tmod,varname):
                    found=1
                    val=getattr(tmod,varname)
                    #now mark this object as used...
                    if not hasattr(self.used,module):
                        setattr(self.used,module,This())
                    umod=getattr(self.used,module)
                    setattr(umod,varname,1)
                    break
        if found==0:
            if val==None and raiseerror==1:
                print "ERROR: value not found %s"%str(varname)
                raise Exception("ERROR: value not found: %s %s"%(str(varname),str(searchOrder)))
            else:
                if warn:
                    print "INFORMATION: using default value of **%s** for **%s**, not found in: %s"%(str(val),str(varname),str(searchOrder))
        return val

    def setVal(self,varname,value,searchOrder=None,raiseerror=1):
        """Change a stored value. This does not affect the XML file, only the
        stored representation of it.
        @param varname:  The variable to be changed
        @type varname: String
        @param value: The new value of the variable to be changed
        @type value: User defined
        @param searchOrder: The order in which to search modules for the variable (None to use default)
        @type searchOrder: List of Strings
        @param raiseerror: Whether to raise an error if variable not found in searchPath
        @type raiseerror: Int
        """
        if searchOrder==None:
            searchOrder=self.searchOrder
        found=0
        for module in searchOrder:
            if hasattr(self.this,module):
                tmod=getattr(self.this,module)
                if hasattr(tmod,varname):
                    setattr(tmod,varname,value)
                    found=1
                    break
        if found==0 and raiseerror==1:
            print "ERROR: value not found",varname
            raise Exception("Error: value not found")

    def strUsed(self,indent="    "):
        """Get all used objects"""
        txt="#Warning, if an object is used within the config file, but not in simulation, it won't appear here - you should be aware of that.\n"
        objs=dir(self.used)
        for obj in objs:
            if obj[:2]!="__" and obj!="simID" and obj!="filename":
                txt+="%s\n"%obj
                vars=dir(getattr(self.used,obj))
                for var in vars:
                    if var[:2]!="__":
                        #the variable is used...
                        txt+="%s%s\n"%(indent,var)
        return txt

    def strUnused(self,indent="    "):
        """return all unused objects as a string"""
        txt="#Warning - just because an object appears unused, it may be used within the config file itself - you should check that before deleting it.\n"
        objs=dir(self.this)
        for obj in objs:
            if obj[:2]!="__" and obj!="simID" and obj!="filename":
                tobj=getattr(self.this,obj)
                if hasattr(self.used,obj):
                    uobj=getattr(self.used,obj)
                else:
                    uobj=This()
                otxt="%s\n"%obj
                vtxt=""
                for var in dir(tobj):
                    if var[:2]!="__":
                        if not hasattr(uobj,var):
                            #object is not used...
                            vtxt+="%s%s\n"%(indent,var)
                if len(vtxt)>0:
                    txt+=otxt+vtxt
        return txt
    def prepareSchema(self,name,attrs):
        """Prepare the Schema python script from part of the XML file and
        write it.
        @param name: The object name, not used
        @type name: String
        @param attrs: Dictionary of tag attributes
        @type attrs: Dict
        """
        obj=None
        params=""
        oname=None
        mod=None
        if attrs.has_key("mod"):
            self.schemaFile.write("import "+attrs["mod"]+"\n")
            obj=attrs["mod"]
            mod=attrs["mod"]
        if attrs.has_key("obj"):
            obj=attrs["obj"]
        if attrs.has_key("params"):
            params=attrs["params"]
        if attrs.has_key("name"):
            oname=attrs["name"]
        if oname!=None and obj!=None and mod!=None:
            self.schemaFile.write(oname+"="+mod+"."+obj+"("+params+")\n")

    def writeOut(self,filename=None):
        """Write the current settings to an xml file.  This can be used to
        duplicate a parameter XML file, or to save changes.
        @param filename: The filename (if None, use current filename)
        @type filename: None or String
        """
        if filename==None:
            filename=self.filename
        txt=self.writeXML()
        f=open(filename,"w")
        f.write(txt)
        f.close()

    def writeXML(self):
        """Prepare the XML string for writing to a file, using current internal
        module and variable information.
        @return: String for writing to XML file
        @rtype: String
        """
        s=""
        s+='<?xml version="1.0"?>\n<aosim>\n<author name="%s"/>\n<created date="%s"/>\n'%(os.environ["USER"],time.strftime("%y%m%d",time.localtime()))
        for module in self.fileData.modules:
            s+="\n<module"
            for key in module.args.keys():
                marg=string.replace(module.args[key],'"',"'")
                s+=' %s="%s"'%(key,marg)
            s+=">\n<variables>\n"
            for vars in module.variableList:
                var=vars[0]
                s+="<var"
                orderList=["name","type","value","comment"]
                for key in orderList:
                    if key in var.keys():
                        if key!="value" or (key=="value" and not var.has_key("extendedTxt")):
                            vvar=string.replace(str(var[key]),'"',"'")
                            s+=' %s="%s"'%(key,vvar)
                        
                for key in var.keys():
                    if key not in orderList:#!="value" and key!="extendedTxt":
                        vvar=string.replace(str(var[key]),'"',"'")
                        s+=' %s="%s"'%(key,vvar)
                if var.has_key("extendedTxt"):
                    if var["value"][0]=="\n":
                        s+=">%s"%var["value"]
                    else:
                        s+=">\n%s"%var["value"]
                    if var["value"][-1]=="\n":
                        s+="</var>\n"
                    else:
                        s+="\n</var>\n"
                else:
                    s+="/>\n"#' value="%s"/>\n'%string.replace(var["value"],'"',"'")
            s+="</variables>\n"
            s+="</module>\n"
        s+="</aosim>\n"
        return s
    def newVar(self,moduleName,attribs,moduleArgs=None,pos=-1):
        """Add a new variable to the current settings (used by GUI).
        pos determines where in the xml file the new variable is added -
        -1 means at end of current vars in the module, 0 means at beginning
        etc.
        Does not affect self.this but only the fileData object.
        @param moduleName: The name of the module inwhich to insert the variable
        @type moduleName: String
        @param attribs: Attributes for the XML tag
        @type attribs: Dict
        @param moduleArgs: Optional dictionary to specify which module to insert (allows specification of batch number etc).
        @type moduleArgs: None or Dict
        @param pos: The position at which to insert the new variable
        @type pos: Int
        """
        if moduleArgs==None:
            moduleArgs={"name":moduleName}
        done=0
        for module in self.fileData.modules:
            if module.name==moduleName and module.args==moduleArgs:
                done=1
                if pos>=0:
                    module.variableList.insert(pos,(attribs,"UPDATE"))
                else:
                    module.variableList.append((attribs,"UPDATE"))
                break
        if done==0:#add new module
            self.fileData.modules.append(ModuleClass(moduleName,moduleArgs))
            self.fileData.modules[-1].variableList.append((attribs,"UPDATE"))
    def deleteVar(self,moduleName,attribs,moduleArgs=None):
        """Detete a variable from the current settings.  This will affect the filedata and not self.this.
        @param moduleName: Name of module
        @type moduleName: String
        @param attribs: Attributes of the variable to be delected
        @type attribs: Dict
        @param moduleArgs: Specifics of module to search
        @type moduleArgs: None or Dict
        """
        if moduleArgs==None:
            moduleArgs={"name":moduleName}
        done=0
        for module in self.fileData.modules:
            if module.name==moduleName and module.args==moduleArgs:#module ok
                if attribs==None:#delete whole module
                    self.fileData.modules.remove(module)
                    done=1
                    break
                else:#delete var in this module
                    for var in module.variableList:
                        if var[0]==attribs:
                            module.variableList.remove(var)
                            done=1
                            break
                if done==1:
                    break
        if done==0:
            raise Exception("Module/var not found to delete",moduleName,attribs)
        

    def changeVar(self,moduleName,attribs,moduleArgs=None):
        """Change attributes for a given variable.
        @param moduleName: Name of module to search
        @type moduleName: String
        @param attribs: Attributes of variable to change
        @type attribs: Dict
        @param moduleArgs: Specifics of module to search
        @type moduleArgs: Dict or None"""
        if moduleArgs==None:
            moduleArgs={"name":moduleName}
        done=0
        for module in self.fileData.modules:
            if module.name==moduleName and module.args==moduleArgs:#module ok
                for var in module.variableList:
                    if var[0].has_key("name") and attribs.has_key("name"):
                        if var[0]["name"]==attribs["name"]:
                            done=1
                            for key in attribs.keys():
                                var[0][key]=attribs[key]
                            for key in var[0].keys():
                                if not attribs.has_key(key):
                                    del(var[0][key])
                            var[1]="UPDATE"
                            break
                    else:
                        raise Exception("No var name found in args/attribs")
                    if done==1:
                        break
                if done==1:
                    break
        if done==0:
            raise Exception("Variable doesn't exist - cannot change:"+moduleName)
                                

    def postAdd(self,namedata):#used internally...
        """Add a tuple to the self.postList"""
        #print "postadd %s"%namedata[0]
        remlist=[]
        for nd in self.postList:
            if nd[0]==namedata[0]:
                remlist.append(nd)
        for r in remlist:
            self.postList.remove(r)
        self.postList.insert(0,namedata)
    def postGet(self,name,default=None,raiseerror=1):
        """Get data from the self.postList.  This is data that has been shared by other DASP science modules"""
        #print "postget %s"%name
        for n,d in self.postList:
            if n==name:
                return d
        if default==None and raiseerror==1:
            raise Exception("config.postGet could not find name %s"%name)
        return default
    
    def post(self,name,data):
        """Share data variable name between all MPI processes...
        This should be used for one-off sharing of data, ie data that is computed only once, and typically will be called from the GUI... No module should call this automatically.  With the possible exception of infScrn posting the initial screens for infAtmos to get...
        """
        self.postAdd((name,data))
        #Now share to other MPI processes... (by socket...)
        for i in range(self.rankSize):
            if i!=self.rank:#not our rank, so send...
                if self.connectionParamsDict.has_key(i):#we can connect to this process...
                    host,port=self.connectionParamsDict[i]
                    conn=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    try:
                        conn.connect((host,port))
                        print "config.post: Connected to %s %d"%(host,port)
                    except:
                        print "util.readConfig - Couldn't connect"
                        conn=None
                    if conn!=None:
                        lst=["add",(name,data),None,None]
                        #print "serialise.Send"
                        util.serialise.Send(lst,conn)
                        #print "serialise.Sent..."
                        #try:
                        #    util.serialise.Send(lst,conn)
                        #except:
                        #    print "util.readConfig - Unable to post data %s - serialise failed"%name
                        conn.close()
    def getCpus(self):
        ncpu=None
        try:
            lines=open("/proc/cpuinfo").readlines()
            ncpu=0
            for line in lines:
                if ("processor" in line) and (":" in line):
                    ncpu+=1
            if ncpu==0:
                print "Warning - couldn't determine number of CPUs, assuming 1"
                ncpu=1
        except:
            pass
        if ncpu==None:
            try:
                cmd="sysctl -a hw | grep ncpu | tail -n 1"
                print "/proc/cpuinfo not found - may be OSX?  Trying commandline '%s'"%cmd
                import popen2
                line=popen2.popen2(cmd)[0].read().strip()
                ncpu=int(line.split("=")[1])
            except:
                ncpu=1
                print "Cannot detemine number of CPUs - assuming 1"
        return ncpu
    
class This:
    """A class to hold values obtained fromt the XML file as it is parsed"""
    def __init__(self):
        """Does nothing"""
        pass

class FileData:
    """A class to hold information (including values) obtained from the XML file as it is parsed, to be used e.g. for a GUI
    @cvar modules: List of currently parsed modules
    @type modules: List of ModuleClass objects"""
    def __init__(self):
        """Initialises the module list to empty"""
        self.modules=[]
class ModuleClass:
    """A class to hold information about modules and variables contained within.
    @cvar name: Name of module
    @type name: String
    @cvar args: Arguments to create module
    @type args: Dict
    @cvar variableList: List of variables within module
    @type variableList: List of Dict
    """
    def __init__(self,name,args):
        self.name=name
        self.args=args
        self.variableList=[]#a list of tuples of (dictionaries which contain the attributes of each variable, evaluated value).
    def copy(self):
        m=ModuleClass(self.name,self.args.copy())
        for v in self.variableList:
            m.variableList.append((v[0].copy(),v[1]))
        return m
class dudFile:
    """An empty class containing methods to mimic a file object, but do nothing
    with data read or written"""
    def __init__(self):
        """pass"""
        pass
    def write(self,data):
        """pass"""
        pass
    def read(self,data):
        """pass"""
        pass
    def close(self):
        """pass"""
        pass
    def flush(self):
        """pass"""
        pass


class PreFormatXML:
    """Used to expand modules with id tags into separate modules."""
    def __init__(self,file,batchno=0):
        self.p=None
        self.file=file
        self.filename=None
        self.batchno=batchno
        self.reset()
        self.p = xml.parsers.expat.ParserCreate()
        self.p.StartElementHandler = self.start_element
        self.p.EndElementHandler = self.end_element
        self.p.CharacterDataHandler = self.char_data
        self.p.returns_unicode=0
        self.open(file)

    def reset(self):
        pass
    def open(self,file):
        """Open an XML file and parse it.  Expand any module tags with id arguments, and then save the result in self.txt
        @param file: Filename
        @type  file: String
        """
        txt=open(file).read()
        self.filename=file
        self.error=0
        self.p.Parse(txt)
        if self.error:
            self.reset()

    def reset(self):
        """Reset the object ready to parse a different file"""
        self.error=0
        self.tagOpenList=[]
        self.txt=""
        self.p = xml.parsers.expat.ParserCreate()
        self.p.StartElementHandler = self.start_element
        self.p.EndElementHandler = self.end_element
        self.p.CharacterDataHandler = self.char_data
        self.p.returns_unicode=0
    def start_element(self,name, attrs):
        """Called by XML parser at opening of tag.
        @param name: Name of tag
        @type  name: String
        @param attrs: Attributes set in tag
        @type attrs: Dict
        """
        self.tagOpenList.append(Tag(name,attrs))
        #print name
    def end_element(self,name):
        """Called by XML parser at end of element
        @param name: Name of tag to be closed
        @type name: String
        """
        tag=self.tagOpenList.pop()
        if len(self.tagOpenList)>0:
            self.tagOpenList[-1].tags.append(tag)
        else:#have finished parsing...
            self.txt=self.finalise(tag)
        #print "/"+name
    def char_data(self,data):
        """Called while parsing a multi-line tag.  Stores the data internally.
        @param data: The data to be stored
        @type data: String
        """
        tag=self.tagOpenList[-1]
        tag.tags.append(data)
        #if tag.storedtxt==None:
        #    tag.storedtxt=data
        #else:
        #    tag.storedtxt+=data

    def finalise(self,tag):
        """Run through tag and rewrite the tags..."""
        if type(tag)==type(""):#character data only...
            txt=tag
        else:#get this tag, and everything under it...
            id=[]
            if tag.name=="module" and tag.attrs.has_key("id"):
                id=eval(tag.attrs["id"])
                if type(id)==type(()):
                    id=list(id)#convert tuple to list.
                if type(id)!=type([]):
                    id=[id]
            if len(id)==0:
                id=[""]
            
            intxt=""
            for t in tag.tags:
                intxt+=self.finalise(t)
            txt=""
            for i in id:
                attrs=tag.attrs.copy()
                if len(i)>0:#no specific id for module name
                    attrs["name"]+="_"+i
                txt+="<%s"%tag.name
                for attr in attrs.keys():
                    txt+=' %s="%s"'%(attr,attrs[attr])
                if len(tag.tags)==0:
                    txt+="/>"
                else:
                    txt+=">"+intxt+"</%s>"%tag.name
                if len(id)>1:
                    txt+="\n"
                    
            #txt="<"+tag.name
            #for attr in tag.attrs.keys():
            #    txt+=' %s="%s"'%(attr,tag.attrs[attr])
            #if len(tag.tags)==0:
            #    txt+="/>"
            #else:
            #    txt+=">"
            #    for t in tag.tags:
            #        txt+=self.finalise(t)
            #    txt+="</%s>"%tag.name
        #print "finalising: ",txt
        return txt

class Tag:
    def __init__(self,name,attrs):
        self.name=name
        self.attrs=attrs
        #self.storedtxt=None
        self.tags=[]



if __name__=="__main__":
    """Code for testing purposes"""
    x=AOXml(batchno=9)
    import sys
    if len(sys.argv)==1:
        txt1="""<?xml version="1.0"?>
        <aosim>
          <author name="Alastair Basden"/>
          <created date="050322"/>
          <schema>
          <comment value="Hi there"/>
          </schema>
          <module name="globals">
            <variables>
              <var name="l0" type="f" value="30.0"/>
              <var name="atmosLinearInterp" type="i" value="2"/>
              <var name="npup" type="i" value="64"/>
              <var name="nwfs" type="i" >28</var>
            </variables>
          </module>
          <module name="atmos" batchno="range(3,10)">
            <variables>
              <var name="atmosLinearInterp" type="i" value="1"/>
              <var name="tstep" type="eval" value="0.005*this.batchNumber"/>
              <var name="arr" type="Numeric" value="Numeric.zeros((this.atmos.atmosLinearInterp,10),'f')"/>
              <var name="ntel" type="copy" value="this.globals.npup"/>
              <var name="ntel2" type="code" value="ntel2=this.globals.npup*this.globals.l0**2"/>
              <var name="wierd" type="abcxyz" value="type is abcxyz which isn't known, so is assumed to be a string"/>
            </variables>
          </module>
          <module name="atmos" batchno="2">
            <variables>
              <var name="atmosLinearInterp" type="i" value="122"/>
              <var name="ntel" type="copy" value="this.globals.npup*2"/>
            </variables>
          </module>
          <module name="atmos">
            <variables>
              <var name="atmosLinearInterp" type="i" value="1122"/>
              <var name="ntel3" type="copy" value="this.globals.npup*4"/>
            </variables>
          </module>
        </aosim>
        """
    else:
        try:
            txt1=open(sys.argv[1]).read()
        except:
            print "Couldn't open requested file - using default file"
            txt1=open("/home/ali/py/mpi/params.xml").read()
    x.parse(txt1)

    print x
    print "testing getVal:"
    print "ntel=",x.getVal("ntel")
    print "atmosLinearInterp=",x.getVal("atmosLinearInterp")
    print "atmosLinearInterp=",x.getVal("atmosLinearInterp",searchOrder=["globals","atmos"])
    print "atmosLinearInterp=",x.getVal("atmosLinearInterp",99.,searchOrder=["badmod"])
    print "tstep=",x.getVal("tstep"),"defined as 0.005*batchNumber for some batches"
    writeFile=0
    if writeFile:
        print "testing writeout to /tmp/tst.xml"
        x.writeOut("/tmp/tst.xml")
