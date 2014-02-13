"""Code to parse a sim setup xml file, and create objects from it.  These can then be used to create a python simulation script"""
import xml.parsers.expat


class parseSimXml:
    def __init__(self,txt=""):
        self.p = xml.parsers.expat.ParserCreate()
        self.p.StartElementHandler = self.startElement
        self.p.EndElementHandler = self.endElement
        self.p.CharacterDataHandler = self.charData
        self.p.returns_unicode=0
        self.initialise()
        self.parse(txt)
    def initialise(self):
        self.objList=[]
        self.intagList=[]
        self.storedTxt=None
        self.mpiTag=0
        self.shmCnt=0
        self.globalPreCode=""
        self.globalPostCode=""
    def parse(self,txt):
        if len(txt.strip())>0:
            self.p.Parse(txt)
        self.sortSharing()
    def sortSharing(self):
        #now sort out the shareTags stuff...
        for obj in self.objList:
            if obj.type=="simObj":
                for stag in obj.sharedTo:
                    sobj=self.getObjWithTag(stag)
                    if obj.shareTags==None:
                        if sobj.shareTags==None:
                            obj.shareTags=[]
                            sobj.shareTags=obj.shareTags
                        else:
                            obj.shareTags=sobj.shareTags
                    else:
                        if sobj.shareTags==None:
                            sobj.shareTags=obj.shareTags
                        else:
                            if obj.shareTags is not sobj.shareTags:
                                tmp=sobj.shareTags
                                for t in tmp:
                                    if t not in obj.shareTags:
                                        obj.shareTags.append(t)
                                    self.getObjWithTag(t).shareTags=obj.shareTags
                                sobj.shareTags=obj.shareTags
                    if obj.tag not in obj.shareTags:
                        obj.shareTags.append(obj.tag)
                    if stag not in obj.shareTags:
                        obj.shareTags.append(stag)
                    
    def startElement(self,name,attrs):
        self.intagList.append(name)
        #if len(self.intagList)>2 and self.intagList[-2]=="simSetup" and self.intagList[-3]=="aosim":
        if "aosim" in self.intagList and "simSetup" in self.intagList:
            self.storedTxt=None
            if name=="simulationObject":
                self.simObjAttrs=attrs
                self.objList.append(simObj(dictionary=attrs))
            elif name=="groupShare":
                self.simObjAttrs=attrs
                self.objList.append(groupShare(dictionary=attrs))
    def charData(self,data):
        if self.storedTxt==None:
            self.storedTxt=data
        else:
            self.storedTxt+=data
    def endElement(self,name):
        n=self.intagList.pop()
        if n!=name:
            print "ERROR Parse error - muddled up somehow %s %s"%(n,name)
            raise Exception("XML parse error")
        if name=="simulationObject":
            o=self.objList[-1]
            o.finalise()
        elif name=="groupShare":
            o=self.objList[-1]
            o.finalise()
        elif len(self.intagList)>0:
            if self.intagList[-1]=="simulationObject":
                o=self.objList[-1]
                if name=="precode":
                    o.precodetxt=self.storedTxt.strip()
                elif name=="postcode":
                    o.postcodetxt=self.storedTxt.strip()
                elif name=="lines":
                    o.linestxt=self.storedTxt.strip()
                elif name=="endlines":
                    o.endlinestxt=self.storedTxt.strip()
                elif name=="parentNames":
                    o.parentNamesTxt=self.storedTxt.strip()
                elif name=="sharedTo":
                    o.sharedToTxt=self.storedTxt.strip()
            elif self.intagList[-1]=="simSetup":
                if name=="precode":
                    self.globalPreCode=self.storedTxt.strip()
                elif name=="postcode":
                    self.globalPostCode=self.storedTxt.strip()

    def expandGroups(self):
        """Here, we expand any grouped objects into their physical components."""
        groupList=[]
        objList=[]
        for o in self.objList:#separate the real sci objects from the group objects
            if o.type=="groupShare":
                groupList.append(o)
            else:
                objList.append(o)
        groupedObjList=[]#a tiered list of all (unexpanded) sci objects that are grouped
        expandedObjList=[]#a tiered list of all expanded sci objects that are grouped.  
        allGroupedObjects=[]#a flat list of all (unexpanded) sci objects that are grouped
        #now get a list of all the grouped objects.
        for g in groupList:
            gol=[]
            eol=[]
            groupedObjList.append(gol)
            expandedObjList.append(eol)
            for o in objList:
                if g.contains(o):
                    if o in groupedObjList:
                        raise Exception("Object grouped more than once")
                    gol.append(o)
                    allGroupedObjects.append(o)
                    eol.append([])
        print "All grouped object:",allGroupedObjects
        #now, for all the grouped objects, go through expanding them
        #groupedObjList is a list of lists.  Each of these lists is a list of grouped objects in a certain group.
        #expandedObjList contains a list for each entry of groupedObjList, ie the expanded objects...
        allObjList=[]
        idlistList=[]
        for i in range(len(groupList)):#for each group of objects...
            g=groupList[i]#this is the group.
            gol=groupedObjList[i]#this is a list of objects in the group.
            eol=expandedObjList[i]
            g.idstr=g.idstr.replace(","," ")
            tmp=g.idstr.split(" ")#get the group ID strings.
            idlist=[]
            for i in tmp:
                if len(i.strip())>0:
                    idlist.append(i.strip())
            #print "group ID list",idlist
            ningroup=len(idlist)
            idlistList.append(idlist)
            try:
                cpulist=eval(g.cpuList)
            except:
                print "ERROR Cannot evaluate group CPU list %s"%str(g.cpuList)
                cpulist=[]
            if type(cpulist)!=type([]):
                cpulist=[cpulist]
            if len(cpulist)==0:
                cpulist=[(1,1)]
            while len(cpulist)<ningroup:
                cpulist+=cpulist
            cpulist=cpulist[:ningroup]
            for j in range(len(cpulist)):
                if type(cpulist[j])==type(1):
                    cpulist[j]=(cpulist[j],1)
                if type(cpulist[j])!=type(()):
                    print cpulist,j
                    raise Exception("Unknown cpu type %s"%str(cpulist[j]))
            #the cpu list should now be correct for this group.
            #Now, for every object, expand, give each the correct idstr, work out parent tags and resource sharing.
            for j in range(len(gol)):
                o=gol[j]
                o.expandedTo=[]
                for k in range(ningroup):
                    cpu=cpulist[k]
                    if o.idstr==None or o.idstr=="":#if no idstr defined for this obj, use the group idstr.
                        idstr=idlist[k]
                    else:#if obj has an idstr defined, replace $ with the group idstr and use.
                        idstr=o.idstr
                        idstr=idstr.replace("$",idlist[k])
                    s=simObj(cpu=cpu,imp=o.imp,object=o.object,pos=(0,0),tag=self.newTag(),shortname=o.shortname,pixmap=None,feedback=o.feedback,pyname="",args=o.args,connectto=o.connectto[:],connectfrom=o.connectfrom[:],dictionary=None,mpiTag=None,shmname="",textcol="white",idstr=idstr,groupShare=o.groupShare)
                    eol[j].append(s)
                    allObjList.append(s)
                    s.precodetxt=o.precodetxt
                    s.postcodetxt=o.postcodetxt
                    s.parentNamesTxt=o.parentNamesTxt#.replace("$",idlist[k]+"=")
                    s.sharedToTxt=o.sharedToTxt
                    s.endlinestxt=o.endlinestxt
                    s.linestxt=o.linestxt
                    s.parentNames=[]
                    o.expandedTo.append(s)
                    o.idlist=idlist
                    for pn in o.parentNames:
                        s.parentNames.append(pn)#.replace("$",idlist[k]))
                    #s.parentNames=o.parentNames[:]
                    s.finalise()
        #now we've expanded all the objects, need to work through and update sharing and connections.
        for i in range(len(groupList)):#for each group of objects...
            eol=expandedObjList[i]
            gol=groupedObjList[i]
            idlist=idlistList[i]
            ningroup=len(idlist)
            for j in range(len(gol)):
                go=gol[j]
                ningroup=len(eol[j])
                for k in range(ningroup):
                    s=eol[j][k]
                    
                    #now we need to sort out objects that shareTo o.tag, and objects that connectto or connectfrom o.tag.
                    #We also need to sort out this objects shareTo and connectto/from.
                    sharedTo=[]
                    for stag in s.sharedTo:
                        shareObj=self.getObjWithTag(stag)
                        if shareObj in allGroupedObjects:#need to replace this with all the expanded ones...
                            for ii in range(len(groupList)):#first of all, find it...
                                eoltmp=expandedObjList[ii]
                                for gotmp in groupedObjList[ii]:
                                    #print shareObj,gotmp
                                    if shareObj is gotmp:
                                        tmplist=eoltmp[groupedObjList[ii].index(shareObj)]#this is the expanded list...
                                        for tmpobj in tmplist:
                                            if tmpobj.cpu==s.cpu:
                                                if (tmpobj.tag not in sharedTo) and (s.tag not in tmpobj.sharedTo):
                                                    sharedTo.append(tmpobj.tag)
                                
                        else:
                            if shareObj.cpu==s.cpu:#sharing allowed.
                                if stag not in sharedTo:
                                    sharedTo.append(stag)
                    s.sharedTo=sharedTo
                    if s.groupShare:#share with the other repeats of this object...
                        for tmpobj in eol[j]:
                            if s is not tmpobj and s.cpu==tmpobj.cpu:#if on same CPU
                                if (tmpobj.tag not in s.sharedTo) and (s.tag not in tmpobj.sharedTo):#if not already sharing
                                    print "Appending share tag %s to %s"%(str(tmpobj.tag),str(s.tag))
                                    s.sharedTo.append(tmpobj.tag)
        #now check whether other objs shareto any of the group objs.
        for i in range(len(self.objList)):
            remlist=[]
            o=self.objList[i]
            if o.type=="simObj":
                if o not in allGroupedObjects:#not a shared object
                    for stag in o.sharedTo:
                        so=self.getObjWithTag(stag)#get the objects it shares with
                        if so in allGroupedObjects:
                            #we are trying to share to a group object...
                            for ii in range(len(groupList)):#first of all find it
                                eoltmp=expandedObjList[ii]
                                for gotmp in groupedObjList[ii]:
                                    print gotmp
                                    print so
                                    if so in gotmp:
                                        tmplist=eoltmp[gotmp.index(so)]#this is the expanded list
                                        for tmpobj in tmplist:
                                            if tmpobj.cpu==o.cpu:
                                                if (tmpobj.tag not in o.sharedTo) and (o.tag not in tmpobj.sharedTo):
                                                    o.sharedTo.append(tmpobj.tag)
                            remlist.append(stag)
                    for t in remlist:
                        o.sharedTo.remove(t)

        #now update any connections...
        for i in range(len(groupList)):#for each group of objects...
            eol=expandedObjList[i]
            gol=groupedObjList[i]
            idlist=idlistList[i]
            ningroup=len(idlist)
            for j in range(len(gol)):#for each object with in this group...
                go=gol[j]
                ningroup=len(eol[j])
                for k in range(ningroup):
                    s=eol[j][k]
                    print "connect from",s.connectfrom,"connect to",s.connectto,"s.tag",s.tag,"go.tag",go.tag
                    #now look at the connections...
                    for l in range(len(s.connectfrom)):
                        ftag=s.connectfrom[l]
                        tmpobj=self.getObjWithTag(ftag,self.objList+allObjList)
                        if tmpobj==None:
                            print "Object not found with tag",ftag
                        #is tmpobj in the group?  If so, need to update to connnect from the correct object if its been created so far, and then change its connectto.
                        #If not, need to update tmpobj, so that its connectto includes this object.
                        if tmpobj in gol:#connecting from something in the group.
                            print "tmpobj in gol",tmpobj
                            expandedParents=eol[gol.index(tmpobj)]
                            if s.tag not in expandedParents[k].connectto:
                                expandedParents[k].connectto.append(s.tag)
                            if expandedParents[k].tag not in s.connectfrom:
                                s.connectfrom.append(expandedParents[k].tag)
                                s.parentNames.append(s.parentNames[s.connectfrom.index(ftag)])#go.parentNames[xxx])
                        elif tmpobj in allGroupedObjects:
                            #connects to a grouped object in another group.  Need to find everything this has been expanded too...
                            for m in range(len(tmpobj.expandedTo)):
                                tmpobj2=tmpobj.expandedTo[m]
                                id=tmpobj.idlist[m]
                                if s.tag not in tmpobj2.connectto:
                                    tmpobj2.connectto.append(s.tag)
                                if tmpobj2.tag not in s.connectfrom:
                                    s.connectfrom.append(tmpobj2.tag)
                                    #s.parentNames.append(s.parentNames[s.connectfrom.index(ftag)])
                                    s.parentNames.append(s.parentNames[s.connectfrom.index(ftag)].replace("$",id))
                        else:#update tmpobj connectto to connect to this obj.
                            if s.tag not in tmpobj.connectto:
                                tmpobj.connectto.append(s.tag)
                            if tmpobj.tag not in s.connectfrom:
                                s.connectfrom.append(tmpobj.tag)
                                s.parentNames.append(s.parentNames[s.connectfrom.index(ftag)])
                    for ttag in s.connectto:
                        tmpobj=self.getObjWithTag(ttag,self.objList+allObjList)
                        if tmpobj==None:
                            print "Object not found with tag",ttag
                        #is tmpobj in the group?  If so, and it has already been expanded, need to update to connect to the correct one, and update this one.
                        #If not, need to update tmpobj, so that its connectfrom includes this object.
                        if tmpobj in gol:#connecting to something in the group.
                            expandedParents=eol[gol.index(tmpobj)]
                            if s.tag not in expandedParents[k].connectfrom:
                                expandedParents[k].connectfrom.append(s.tag)
                                expandedParents[k].parentNames.append(expandedParents[k].parentNames[expandedParents[k].connectfrom.index(go.tag)])
                            if expandedParents[k].tag not in s.connectto:
                                s.connectto.append(expandedParents[k].tag)
                        elif tmpobj in allGroupedObjects:
                            #connects from another group.  Need to find everything this has been expanded too..
                            for tmpobj2 in tmpobj.expandedTo:
                                if s.tag not in tmpobj2.connectfrom:
                                    tmpobj2.connectfrom.append(s.tag)
                                    tmpobj2.parentNames.append(tmpobj2.parentNames[tmpobj2.connectfrom.index(go.tag)].replace("$",idlist[k]))#+"!"))
                                if tmpobj2.tag not in s.connectto:
                                    s.connectto.append(tmpobj2.tag)
                        else:#update tmpobj connectfrom to connect to this obj.
                            print "tmpobj not in gol",tmpobj,k,idlist,len(groupList),i
                            if s.tag not in tmpobj.connectfrom:
                                tmpobj.connectfrom.append(s.tag)
                                tmpobj.parentNames.append(tmpobj.parentNames[tmpobj.connectfrom.index(go.tag)].replace("$",idlist[k]))
                            if tmpobj.tag not in s.connectto:
                                s.connectto.append(tmpobj.tag)
                                
                        
        #finally, remove all the group objects from the list, and update everything connecting to/from them and sharing with them.
        for go in allGroupedObjects:
            for o in self.objList+allObjList:
                if o.type=="simObj" and o not in allGroupedObjects:
                    if go.tag in o.sharedTo:
                        o.sharedTo.remove(go.tag)
                    if go.tag in o.connectfrom:
                        indx=o.connectfrom.index(go.tag)
                        o.connectfrom.pop(indx)
                        o.parentNames.pop(indx)
                        #o.connectfrom.remove(go.tag)
                    if go.tag in o.connectto:
                        o.connectto.remove(go.tag)
        for go in allGroupedObjects:
            self.objList.remove(go)
        for g in groupList:
            self.objList.remove(g)
        for i in range(len(groupList)):
            eol=expandedObjList[i]
            for j in range(len(eol)):
                o=eol[j]
                self.objList+=o
        for obj in self.objList:
            obj.shareTags=None
        self.sortSharing()
        print "OBJLIST:",self.objList
        
    def makePython(self):
        """Make the python simulation file...
        objList.cpu is a tuple of node, processor.  Things on the same
        node but different processor communicate by SHM.
        Things on different nodes use MPI.  If on same node and
        processor will be part of same process.
        Note, this modifies the objects, so if needed for anything else, should be reloaded...
        When using group sharing, cpu for the group object is a list of cpus (either single numbers,
        the nodes, or tuples of node and cpu).
        
        """
        procDict={}#keys are cpus, values are list of objects on this cpu.
        orderedProcDict={}#keys are cpus, vaules are list of object on this cpu, in the order in which they will be created and exec'd.
        importList=["numpy","util.Ctrl","base.mpiGet","base.mpiSend","base.shmGet","base.shmSend","Scientific.MPI"]
        objectList=["newMPIGet","newMPISend","newSHMGet","newSHMSend"]
        MPIRank={}
        MPIRank2={}
        nprocesses=0
        self.freeTag=self.getMaxTag()+1

        #first sort out any grouped objects - ie convert the groups to full objects.
        self.expandGroups()
        
        for o in self.objList:
            if o.type=="simObj":
                if o.cpu not in procDict.keys():
                    procDict[o.cpu]=[o]
                    orderedProcDict[o.cpu]=[]
                    MPIRank[o.cpu]=nprocesses
                    MPIRank2[nprocesses]=o.cpu
                    nprocesses+=1
                else:
                    procDict[o.cpu].append(o)
                if o.imp not in importList:
                    importList.append(o.imp)
                if o.object not in objectList:
                    objectList.append(o.object)
        #First, insert the mpi/shm connections...
        self.insertRemoteConnections(procDict)
        #Now add things to the exec list order...
        self.createProcOrder(procDict,orderedProcDict)

        #Now, things in orderedProcDict are in the correct order for execution.
        ncpu=len(procDict.keys())
        hostlist=""
        for i in range(ncpu):
            cpu=MPIRank2[i]
            hostlist+="n%d-c437 "%cpu[0]
        txt="#mpirun -np %d -hostlist %s /usr/local/bin/mpipython $PWD/thisfile.py\n"%(ncpu,hostlist)
        txt+="#Python code created using the simulation setup GUI...\n"
        txt+="#Order of execution may not be quite optimal - you can always change by hand\n"
        txt+="#for large simulations - typically, the order of sends and gets may not be\n"
        txt+="#quite right.  Anyway, enjoy...\n"

        for i in importList:
            txt+="import %s\n"%i
        txt+="ctrl=util.Ctrl.Ctrl(globals=globals())\n"
        txt+='print "Rank %d imported modules"%ctrl.rank\n'
        txt+="#Set up the science modules...\n"
        for o in objectList:
            txt+="%sList=[]\n"%o
        if len(self.globalPreCode)>0:
            txt+=self.globalPreCode+"\n"
        rank=0
        txt+="#Add any personal code after this line and before the next, and it won't get overwritten\n"
        for key in orderedProcDict.keys():
            objlengths={}
            objtags={}
            rank=MPIRank[key]
            txt+="if ctrl.rank==%d:\n"%rank
            objlist=orderedProcDict[key]
            for obj in objlist:
                if obj.imp=="base.mpiGet":
                    parobjMpiSend=self.getObjWithTag(obj.connectfrom[0])
                    parobj=self.getObjWithTag(parobjMpiSend.connectfrom[0])
                    args=getArgsAsString(parobj.args)
                    sourceRank=MPIRank[parobj.cpu]
                    if parobj.idstr!=None and parobj.idstr!="":
                        if parobj.idstr[0]=='|' and parobj.idstr[-1]=='|':
                            idstr=parobj.idstr[1:-1]
                        else:
                            idstr='"%s"'%parobj.idstr
                    else:
                        idstr=None
                    txt+="    dims,dtype=%s.%s(None,ctrl.config,args=%s,forGUISetup=1,idstr=%s).outputData\n"%(parobj.imp,parobj.object,args,idstr)
                    txt+="    newMPIGetList.append(base.mpiGet.newMPIGet(dims,dtype,%d,%d,ctrl.mpiComm))\n"%(sourceRank,obj.mpiTag)
                    if objlengths.has_key(obj.object):
                        objlengths[obj.object]+=1
                    else:
                        objlengths[obj.object]=1
                    objtags[obj.tag]=objlengths[obj.object]-1#save the position in the list...
                    
                elif obj.imp=="base.mpiSend":
                    chobj=self.getObjWithTag(obj.connectto[0])
                    #parobj=self.getObjWithTag(obj.connectfrom[0])
                    if obj.newParentNeeded:
                        parent="None"
                    else:
                        parent=self.getParentString(obj,objtags)
                    destRank=MPIRank[chobj.cpu]
                    txt+="    newMPISendList.append(base.mpiSend.newMPISend(%s,%d,%d,ctrl.mpiComm))\n"%(parent,destRank,obj.mpiTag)
                    if objlengths.has_key(obj.object):
                        objlengths[obj.object]+=1
                    else:
                        objlengths[obj.object]=1
                    objtags[obj.tag]=objlengths[obj.object]-1#save the position in the list...
                elif obj.imp=="base.shmGet":
                    parobjShmSend=self.getObjWithTag(obj.connectfrom[0])
                    parobj=self.getObjWithTag(parobjShmSend.connectfrom[0])
                    args=getArgsAsString(parobj.args)
                    if parobj.idstr!=None and parobj.idstr!="":
                        if parobj.idstr[0]=='|' and parobj.idstr[-1]=='|':
                            idstr=parobj.idstr[1:-1]
                        else:
                            idstr='"%s"'%parobj.idstr
                    else:
                        idstr=None
                    txt+="    dims,dtype=%s.%s(None,ctrl.config,args=%s,forGUISetup=1,idstr=%s).outputData\n"%(parobj.imp,parobj.object,args,idstr)
                    txt+='    newSHMGetList.append(base.shmGet.newSHMGet(%s,dims,dtype))\n'%(obj.shmname)
                    if objlengths.has_key(obj.object):
                        objlengths[obj.object]+=1
                    else:
                        objlengths[obj.object]=1
                    objtags[obj.tag]=objlengths[obj.object]-1#save the position in the list...

                elif obj.imp=="base.shmSend":
                    parobj=self.getObjWithTag(obj.connectfrom[0])
                    args=getArgsAsString(parobj.args)
                    if parobj.idstr!=None and parobj.idstr!="":
                        if parobj.idstr[0]=='|' and parobj.idstr[-1]=='|':
                            idstr=parobj.idstr[1:-1]
                        else:
                            idstr='"%s"'%parobj.idstr
                    else:
                        idstr=None
                    if obj.newParentNeeded:
                        parent="None"
                    else:
                        parent=self.getParentString(obj,objtags)
                    txt+="    dims,dtype=%s.%s(None,ctrl.config,args=%s,forGUISetup=1,idstr=%s).outputData\n"%(parobj.imp,parobj.object,args,idstr)
                    txt+='    newSHMSendList.append(base.shmSend.newSHMSend(%s,%s,dims,dtype))\n'%(parent,obj.shmname)
                    if objlengths.has_key(obj.object):
                        objlengths[obj.object]+=1
                    else:
                        objlengths[obj.object]=1
                    objtags[obj.tag]=objlengths[obj.object]-1#save the position in the list...
                else:
                    #get the parent string...
                    if obj.newParentNeeded:
                        parent="None"
                    else:
                        parent=self.getParentString(obj,objtags)#string for objects parent...
                    args=getArgsAsString(obj.args)
                    if obj.idstr!=None and obj.idstr!="":
                        idstr='"%s"'%obj.idstr
                        if obj.idstr[0]=='|' and obj.idstr[-1]=='|':
                            idstr=obj.idstr[1:-1]
                        else:
                            idstr='"%s"'%obj.idstr
                    else:
                        idstr=None
                    txt+=obj.precode
                    if obj.pyname=="" or obj.pyname==None:
                        if objlengths.has_key(obj.object):
                            objlengths[obj.object]+=1
                        else:
                            objlengths[obj.object]=1
                        objtags[obj.tag]=objlengths[obj.object]-1#save the position in the list...
                        if obj.shareTags==None or len(obj.shareTags)==0 or type(obj.shareTags[-1])!=type(""):
                            #not resource sharing, or this is the first implementation...
                            txt+="    %sList.append(%s.%s(%s,ctrl.config,args=%s,idstr=%s))\n"%(obj.object,obj.imp,obj.object,parent,args,idstr)
                            if obj.shareTags!=None and len(obj.shareTags)>0:#save the first implementation position if resource sharing.
                                obj.shareTags.append("%sList[%d]"%(obj.object,objtags[obj.tag]))
                        else:#resource sharing...
                            txt+="    %sList.append(%s.addNewIdObject(%s,%s))\n"%(obj.object,obj.shareTags[-1],parent,idstr)
                    else:
                        if objlengths.has_key(obj.pyname):
                            print "warning - overwriting object %s"%obj.pyname
                        objlengths[obj.pyname]=1
                        objtags[obj.tag]=None
                        if obj.shareTags==None or len(obj.shareTags)==0 or type(obj.shareTags[-1])!=type(""):
                            #not resource sharing, or this is the first implementation
                            txt+="    %s=%s.%s(%s,ctrl.config,args=%s,idstr=%s)\n"%(obj.pyname,obj.imp,obj.object,parent,args,idstr)
                            if obj.shareTags!=None and len(obj.shareTags)>0:#save resource sharing position
                                obj.shareTags.append(obj.pyname)
                        else:
                            txt+="    %s=%s.addNewIdObject(%s,%s)\n"%(obj.pyname,obj.shareTags[-1],parent,idstr)
                    txt+=obj.postcode
            execOrder=""
            for obj in objlist:
                if objtags[obj.tag]==None:
                    execOrder+=obj.pyname+","
                else:
                    execOrder+="%sList[%d],"%(obj.object,objtags[obj.tag])
                if obj.newParentNeeded:
                    parent=self.getParentString(obj,objtags)
                    if objtags[obj.tag]==None:
                        thisobj=obj.pyname
                    else:
                        thisobj="%sList[%d]"%(obj.object,objtags[obj.tag])
                    if obj.idstr!=None and obj.idstr!="":
                        if obj.idstr[0]=='|' and obj.idstr[-1]=='|':
                            idstr=obj.idstr[1:-1]
                        else:
                            idstr='"%s"'%obj.idstr
                    else:
                        idstr="None"
                    txt+="    %s.newParent(%s,%s)\n"%(thisobj,parent,idstr)
            txt+="    execOrder=[%s]\n"%execOrder
            txt+="    ctrl.mainloop(execOrder)\n"
        txt+='print "Simulation finished..."\n'
        if len(self.globalPostCode)>0:
            txt+=self.globalPostCode+"\n"
        txt+="#Add any personal code after this, and it will not get overwritten\n"
        #txt+="Scientific.MPI.world.abort(0)\n"
        return txt


    def getParentString(self,obj,objtags):
        """Get the string that represents the parent.  This could be eg:
        {'recon':recon[0],'atmos':atmos[1]}, or atmos[0], or myatmos etc...
        """
        indx=0
        l1=len(obj.connectfrom)
        l2=len(obj.parentNames)
        useDict=1
        if l2==0:
            obj.parentNames=range(l1)
            useDict=0
        else:
            tmp=reduce(lambda x,y:str(x)+str(y),obj.parentNames)
            if type(tmp)==type(""):
                if l1==1 and len(tmp)==0:#no names given
                    useDict=0
            
        l2=l1
        if l2!=l1:
            print "ERROR: number of parent names does not equal number of connections (or zero)"
            raise Exception("ERROR: number of parent names does not equal number of connections (or zero)")
        if useDict==1 or l1>1:#create the dictionary of parents...
            parent="{"
            for fi in range(len(obj.connectfrom)):
                tag=obj.connectfrom[fi]#get tag connecting from
                #and find the name of this object...
                parobj=self.getObjWithTag(tag)
                if objtags[tag]==None:#the object has been given a name
                    pyobj=parobj.pyname
                else:
                    pyobj="%sList[%d]"%(parobj.object,objtags[tag])
                print obj
                if type(obj.parentNames[fi])==type(""):
                    if len(obj.parentNames[fi])>0:
                        strng='"%s"'%obj.parentNames[fi]
                    else:
                        indx+=1
                        strng='"%d"'%indx
                else:
                    strng="%s"%str(obj.parentNames[fi])
                parent+='%s:%s,'%(strng,pyobj)
            parent+="}"
        elif l1==1:#only 1 parent - don't bother with dict...
            tag=obj.connectfrom[0]
            parobj=self.getObjWithTag(tag)
            print tag,objtags.keys(),obj.object,obj.cpu,obj.newParentNeeded,parobj.newParentNeeded,parobj.object
            if obj.object in ["newSHMGet","newMPIGet"]:#these objects aren't asigned parents...
                parent="None"
            elif objtags[tag]==None:#the object has been given a name
                parent=parobj.pyname
            else:
                parent="%sList[%d]"%(parobj.object,objtags[tag])
        else:
            parent="None"#no parents
        return parent

    def createProcOrder(self,procDict,orderedProcDict):
        """Here, we move objects from procDict to orderedProcDict where they are
        ordered in the order that they should be executed.  This uses a recursive algorithm.
        Once an object is placed, all its children are then checked, and so on...
        """
        nadded1=1
        while nadded1>0:
            nadded2=1
            while nadded2>0:
                nadded3=1
                while nadded3>0:
                    nadded4=1
                    while nadded4>0:
                        nadded5=1
                        while nadded5>0:
                            nadded5=self.addObj(procDict,orderedProcDict,tag=None,parentless=1)
                            print "nadded5",nadded5
                        nadded4=self.addObj(procDict,orderedProcDict,tag=None,remote=0)
                        print "nadded4",nadded4
                    nadded3=self.addObj(procDict,orderedProcDict,tag=None,remote="withsend")
                    print "nadded3",nadded3
                nadded2=self.addObj(procDict,orderedProcDict,tag=None,remote="all")
                print "nadded2",nadded2
            nadded1=self.addObj(procDict,orderedProcDict,tag=None,remote="all",ignoreSharing=1)
            print "nadded1",nadded1
        nleft=len(reduce(lambda x,y:x+y,procDict.values()))
        if nleft!=0:
            print "Unallocated objects:"
            for key in procDict.keys():
                print key,procDict[key]
            raise Exception("Not all objects allocated")
        
    def addObj(self,procDict,orderedProcDict,tag=None,remote=0,ignoreSharing=0,parentless=0):
        """Here, add objects to the execution order list, depending on
        what allowed to do...  Remote specifies whether an object on a
        different processor is allowed or not.  ignoreSharing
        specifies whether a resource sharing object should be ignored
        (which might lead to wrong simulation results).

        If parentless is set, will only allocate objects with no parents or with parents who are feedback.
        
        This function adds only one object (per processor), and then
        tries to add all the children of this object recursively.
        
        """
        nadd=0
        for key in procDict.keys():#for each processor...
            procList=procDict[key]
            for obj in procList:
                if tag==None or obj.tag==tag:
                    #initially, we only allocate a resource sharing object if its children can be placed immediately.
                    ready=1
                    canPlaceIfIgnoreShared=0
                    allocated=orderedProcDict[key]
                    fobjs=[]
                    for tmp in procList:
                        if tmp.feedback:
                            fobjs.append(tmp)
                    searchlist=allocated+fobjs
                    #first look for objects with no parents, or feedback as parents...
                    if len(obj.connectfrom)>0:#if this ==0, then has no parents, so can probably be allocated
                        for parenttag in obj.connectfrom:
                            parentobj=self.getObjWithTag(parenttag)
                            if not parentobj.feedback:
                                ready=0#not a feedback object
                                break
                        if ready:
                            obj.newParentNeeded=1

                    if ready:#just check that if this is a resource sharer, children can be allocated.
                        if obj.shareTags!=None and len(obj.shareTags)>0:#a resource sharer...
                            canPlaceIfIgnoreShared=1
                            for childtag in obj.connectto:
                                if not self.canPlaceObj(childtag,searchlist+[obj]):
                                    ready=0
                                    obj.newParentNeeded=0
                                    break
                    #now look for objects with allocated parents...
                    if ready==0 and parentless==0:
                        ready=1
                        obj.newParentNeeded=0
                        for parenttag in obj.connectfrom:
                            parentobj=self.getObjWithTag(parenttag,searchlist)
                            if parentobj==None:#parent not yet allocated or a feedback object
                                ready=0
                                break
                            else:
                                if parentobj.feedback and parentobj.object not in ["newSHMSend","newMPISend"]:#shm/mpi get objects do not have a parent module, since data comes via mpi or shm.
                                    obj.newParentNeeded=1
                        if ready==0: # if not allowed, reset the newParentNeeded flag.
                            obj.newParentNeeded=0
                        else:#this object can be allocated - but if it is resource sharing, and its children can not be allocated immediately, then we won't add it.
                            if obj.shareTags!=None and len(obj.shareTags)>0:
                                canPlaceIfIgnoreShared=1
                                for childtag in obj.connectto:
                                    if not self.canPlaceObj(childtag,searchlist+[obj]):
                                        ready=0
                                        obj.newParentNeeded=0
                                        break

                    if ready==0 and remote!=0:
                        #now see if there is object with remote object as parent.
                        if obj.imp in ["base.mpiGet","base.shmGet"]:#these objects can't be shared, so ignore...
                            if remote=="all":
                                ready=1#the object can be used (we're desparate)
                            elif remote=="withsend":
                                #look to see if the parent has been allocated
                                allallocedlist=reduce(lambda x,y:x+y,orderedProcDict.values())
                                if len(obj.connectfrom)!=1:
                                    raise Exception("get send object with other than one parent...")
                                tmp=self.getObjWithTag(obj.connectfrom[0],allallocedlist)
                                if tmp!=None:
                                    ready=1
                    if ready==0 and ignoreSharing==1 and canPlaceIfIgnoreShared==1:#attempt to break the rules for placing object (we're stuck!).
                        ready=1

                    print ready,parentless,remote,tag,obj

                    if ready==1:#This object can be placed...
                        print "Appending object",obj
                        orderedProcDict[key].append(obj)
                        nadd+=1
                        #only do one object per cpu at a time...
                        procList.remove(obj)
                        #remlist.append(obj)
                        #now only add the harmless objects...
                        #if not parentless:
                        ignoreSharing=0
                        remote=0
                        parentless=0
                        #now add all the children of this object...
                        #First add just those that are mpiSend of shmSend objects, then add other children - we add the send objects first to improve performance (allows other things to begin processing)
                        for ctag in obj.connectto:
                            cobj=self.getObjWithTag(ctag,self.objList)
                            if cobj.imp in ["base.mpiSend","base.shmSend"]:
                                naddpart=self.addObj(procDict,orderedProcDict,tag=ctag,remote=0)#"withsend")
                                nadd+=naddpart
                                print "Added %d childs to %s object "%(naddpart,str(obj))
                        for ctag in obj.connectto:
                            remt=0
                            if obj.imp in ["base.mpiSend","base.shmSend"]:
                                remt=0#"withsend"removed withsend 090218
                            else:
                                remt=0
                            naddpart=self.addObj(procDict,orderedProcDict,tag=ctag,remote=remt)
                            nadd+=naddpart
                            if obj.imp in ["base.mpiSend","base.shmSend"]:
                                print "Added %d childs to %s object "%(naddpart,str(obj))
                        break
        return nadd
    def canPlaceObj(self,tag,objlist):
        """see if this object can be placed...
        objlist is the list of already allocated objects."""
        obj=self.getObjWithTag(tag)
        if obj in objlist:#already allocated...
            return 1
        place=1
        for ptag in obj.connectfrom:
            pobj=self.getObjWithTag(ptag)
            if pobj not in objlist and pobj.feedback==0:
                place=0
        if place:
            if obj.shareTags!=None and len(obj.shareTags)>0:#this is a resource sharing object - treat as allocatable if all its children can also be allocated...
                for ctag in obj.connectto:
                    place=self.canPlaceObj(ctag,objlist+[obj])
                    if place==0:
                        break
        return place

    def createProcOrderOld(self,procDict,orderedProcDict):
        #First those with no parents, or existing parents.
        #Then those with parents that give feedback.
        #Then those with remote parents that have already been allocated.
        #Then others...
        nadded1=1
        while nadded1>0:
            nadded2=1
            while nadded2>0:
                nadded3=1
                while nadded3>0:
                    nadded4=1
                    while nadded4>0:
                        nadded4=self.addObjOld(procDict,orderedProcDict,remote=0,feedback=0)
                    nadded3=self.addObjOld(procDict,orderedProcDict,remote=0,feedback=1)
                nadded2=self.addObjOld(procDict,orderedProcDict,remote="withsend",feedback=1)
            nadded1=self.addObjOld(procDict,orderedProcDict,remote="all",feedback=1)        

    
    def addObjOld(self,procDict,orderedProcDict,remote=0,feedback=0):
        """Add objects... moving from procDict to orderedProcDict - this will then be the order in which they are executed."""
        nadd=0
        for key in procDict.keys():#for each processor...
            procList=procDict[key]
            remlist=[]
            #first look to see if the previously allocated object has children that can be allocated straight away (this is necessary for resource sharing, to ensure that the output is still present!).
            
            if len(orderedProcDict[key])>0:
                cadd=1
                remlist=[]
                while cadd==1:
                    cadd=0
                    parobj=orderedProcDict[key][-1]
                    for cobjtag in parobj.connectto:
                        cobj=self.getObjWithTag(cobjtag,procList)
                        if cobj!=None and cobj not in orderedProcDict[key]:#hasn't currently been allocated...
                            #check that all its parents have been allocated...
                            nalloc=len(cobj.connectfrom)
                            fobjs=[]
                            for tmp in procList:
                                if tmp.feedback:
                                    fobjs.append(tmp)
                            feedbackNeeded=0
                            for partag in cobj.connectfrom:
                                parentObj=self.getObjWithTag(partag,fobjs+orderedProcDict[key])
                                if parentObj!=None:
                                    if parentObj in fobjs:
                                        feedbackNeeded=1
                                    nalloc-=1
                            if nalloc==0:#all parents allocated... so allocate cobj...
                                orderedProcDict[key].append(cobj)
                                nadd+=1
                                cadd+=1
                                remlist.append(cobj)
                                cobj.newParentNeeded=feedbackNeeded
                                break
                for obj in remlist:
                    procList.remove(obj)
            remlist=[]
            for obj in procList:
                ready=1
                allocated=orderedProcDict[key]#reduce(lambda x,y:x+y,orderedProcDict.values())
                
                # now look to see if there is an object with allocated parents
                for parenttag in obj.connectfrom:
                    parentobj=self.getObjWithTag(parenttag,allocated)
                    if parentobj==None:
                        ready=0
                        break
                if ready==0 and feedback==1:
                    #now see if there is object with feedback as parent.
                    fobjs=[]
                    for tmp in procList:
                        if tmp.feedback:
                            fobjs.append(tmp)
                    searchlist=allocated+fobjs
                    ready=1
                    for parenttag in obj.connectfrom:
                        parentobj=self.getObjWithTag(parenttag,searchlist)
                        if parentobj==None:
                            ready=0
                            break
                    if ready==1:#first check we're not a send object (no sense sending before its been made... even for a feedback object...
                        if obj.object in ["newMPISend","newSHMSend"]:
                            ready=0
                        else: # make a note that it was from an unallocated feedback object (a call to newParent() will be needed during sim setup).
                            obj.newParentNeeded=1
                if ready==0 and remote!=0:
                    #now see if there is object with remote object as parent.
                    if obj.imp in ["base.mpiGet","base.shmGet"]:
                        if remote=="all":
                            ready=1#the object can be used (we're desparate)
                        elif remote=="withsend":
                            #look to see if the parent has been allocated
                            allallocedlist=reduce(lambda x,y:x+y,orderedProcDict.values())
                            if len(obj.connectfrom)!=1:
                                raise Exception("get send object with other than one parent...")
                            tmp=self.getObjWithTag(obj.connectfrom[0],allallocedlist)
                            if tmp!=None:
                                ready=1
                if ready==1:#all parent for this object are allocated...
                    #so this object can be added to the list...
                    orderedProcDict[key].append(obj)
                    nadd+=1
                    #only do one object per cpu at a time...
                    remlist.append(obj)
                    if ready==1:# now only add the harmless objects...
                        feedback=0
                        remote=0
##                     #Should also check whether the object just added has any
##                     #remote send children...
##                     for cobjtag in obj.connectto:
##                         cobj=self.getObjWithTag(cobjtag)
##                         if cobj.object in ["newMPISend","newSHMSend"]:
##                             orderedProcDict[key].append(cobj)
##                             nadd+=1
##                             remlist.append(cobj)
##                             if cobj.cpu!=obj.cpu:
##                                 raise Exception("Different CPU for sending object... - algorithm error...")
                    #Also: if the object is an send object, should add the corresponding get object too, so that they come in order... (this may not be the most efficient way, but at least stops simulation freezes - I think).
                    if obj.object in ["newMPISend","newSHMSend"]:
                        for cobjtag in obj.connectto:
                            cobj=self.getObjWithTag(cobjtag)
                            if cobj.object not in ["newMPIGet","newSHMGet"]:
                                raise Exception("Send object not connected to get object")
                            #add the object to the exec order...
                            orderedProcDict[cobj.cpu].append(cobj)
                            cprocList=procDict[cobj.cpu]
                            cprocList.remove(cobj)
                            nadd+=1
                    break
                
            for obj in remlist:
                print obj.cpu,obj.object
                procList.remove(obj)
        print nadd
        return nadd


        
    def insertRemoteConnections(self,procDict):
        procDictAdded={}
        for key in procDict.keys():
            procDictAdded[key]=[]
        for key in procDict.keys():
            print key
            procList=procDict[key]#list of processes on this node (unordered).
            for obj in procList:
                print obj
                childs={}#list of children separated into CPUs
                for childTag in obj.connectto:
                    childObj=self.getObjWithTag(childTag)
                    if childObj==None:
                        print "Couldn't get object with tag",childTag
                    if childs.has_key(childObj.cpu):
                        childs[childObj.cpu].append(childObj)
                    else:
                        childs[childObj.cpu]=[childObj]
                for key2 in childs.keys():#do children on each CPU in turn
                    print "key2:",key2
                    roList=childs[key2]
                    nChilds=len(roList)
                    roCPU=roList[0].cpu
                    if roCPU[0]==obj.cpu[0]:#same node
                        if roCPU[1]!=obj.cpu[1]:#different processor
                            #need SHM objects...
                            ptag=self.newTag()
                            ctag=self.newTag()
                            childTags=[]
                            for c in roList:
                                childTags.append(c.tag)
                                obj.connectto[obj.connectto.index(c.tag)]=ptag
                                c.connectfrom[c.connectfrom.index(obj.tag)]=ctag
                            print "Adding SHM connection",ptag,ctag
                            shmname=self.newSHMName()
                            shmParent=simObj(cpu=obj.cpu,imp="base.shmSend",object="newSHMSend",feedback=obj.feedback,connectto=[ctag],connectfrom=[obj.tag],tag=ptag,shmname=shmname)
                            shmChild=simObj(cpu=roCPU,imp="base.shmGet",object="newSHMGet",feedback=obj.feedback,connectto=childTags,connectfrom=[ptag],tag=ctag,shmname=shmname)
                            procDictAdded[obj.cpu].append(shmParent)
                            procDictAdded[roCPU].append(shmChild)
                    else: #different node, needs MPI
                        ptag=self.newTag()
                        ctag=self.newTag()
                        childTags=[]
                        for c in roList:
                            childTags.append(c.tag)
                            obj.connectto[obj.connectto.index(c.tag)]=ptag
                            c.connectfrom[c.connectfrom.index(obj.tag)]=ctag
                        print "Adding mpi connection",obj.tag,childTags
                        mpiTag=self.getNextMPITag()
                        mpiParent=simObj(cpu=obj.cpu,imp="base.mpiSend",object="newMPISend",feedback=0,connectto=[ctag],connectfrom=[obj.tag],tag=ptag,mpiTag=mpiTag)
                        mpiChild=simObj(cpu=roCPU,imp="base.mpiGet",object="newMPIGet",feedback=0,connectto=childTags,connectfrom=[ptag],tag=ctag,mpiTag=mpiTag)
                        #mpiParent=simObj(cpu=obj.cpu,imp="base.mpiSend",object="newMPISend",feedback=obj.feedback,connectto=[ctag],connectfrom=[obj.tag],tag=ptag,mpiTag=mpiTag)
                        #mpiChild=simObj(cpu=roCPU,imp="base.mpiGet",object="newMPIGet",feedback=obj.feedback,connectto=childTags,connectfrom=[ptag],tag=ctag,mpiTag=mpiTag)
                        procDictAdded[obj.cpu].append(mpiParent)
                        procDictAdded[roCPU].append(mpiChild)
        for key in procDict.keys():
            procDict[key]+=procDictAdded[key]
            self.objList+=procDictAdded[key]
            for obj in procDictAdded[key]:
                print "Added %s object (parent=%s, child=%s)"%(obj.object,str(obj.connectfrom),str(obj.connectto))
    def newTag(self):
        tag=self.freeTag
        self.freeTag+=1
        return tag
    def getNextMPITag(self):
        self.mpiTag+=1
        return self.mpiTag
    def newSHMName(self):
        name='"/aosimSHMMem%%d_%d"%%ctrl.shmtag'%self.shmCnt
        self.shmCnt+=1
        return name
                              
    def getObjWithTag(self,tag,objList=None):
        if objList==None:
            objList=self.objList
        if type(objList)==type({}):
            objList=objList.values()
        obj=None
        for o in objList:
            if o.type=="simObj":
                if o.tag==tag:
                    obj=o
                    break
        return obj
    def getMaxTag(self):
        maxtag=None
        for o in self.objList:
            if o.type=="simObj":
                if maxtag==None or maxtag<o.tag:
                    maxtag=o.tag
        return maxtag

class groupShare:
    def __init__(self,dictionary=None):
        self.type="groupShare"
        self.idstr=""
        self.coordList=[(10,10),(20,10),(20,20),(10,20)]
        self.cpuList="[]"
        if dictionary!=None:
            if dictionary.has_key("cpu"):
                self.cpuList=dictionary["cpu"]
            if dictionary.has_key("coordlist"):
                self.coordList=eval(dictionary["coordlist"])
            if dictionary.has_key("idstr"):
                self.idstr=dictionary["idstr"]
    def finalise(self):
        pass
    def contains(self,o):
        minx=-1
        maxx=0
        miny=-1
        maxy=0
        for x,y in self.coordList:
            if maxx<x:
                maxx=x
            if maxy<y:
                maxy=y
            if minx==-1 or minx>x:
                minx=x
            if miny==-1 or miny>y:
                miny=y
        if o.pos[0]>minx and o.pos[0]<maxx and o.pos[1]>miny and o.pos[1]<maxy:
            return True
        return False
        
    
class simObj:
    def __init__(self,cpu=(1,1),imp=None,object=None,pos=(0,0),tag=None,shortname=None,pixmap=None,feedback=0,pyname="",args="",connectto=[],connectfrom=[],dictionary=None,mpiTag=None,shmname="",textcol="white",idstr=None,groupShare=0):
        self.type="simObj"
        self.precodetxt=""
        self.postcodetxt=""
        self.precode=""
        self.postcode=""
        self.linestxt=""
        self.endlinestxt=""
        self.parentNamesTxt=""
        self.sharedToTxt="[]"
        self.lines=[]
        self.endlines=[]
        self.parentNames=[]
        self.sharedTo=[]
        self.shareTags=None
        self.cpu=cpu
        self.imp=imp
        self.object=object
        self.pos=pos
        self.tag=tag
        self.shortname=shortname
        self.pixmap=pixmap
        self.feedback=feedback
        self.pyname=pyname
        self.groupShare=groupShare
        self.idstr=idstr
        self.args=args
        self.connectto=connectto
        self.connectfrom=connectfrom
        self.textcol=textcol
        self.newParentNeeded=0#call to newParent() for this object - eg for recon.
        self.mpiTag=mpiTag#the MPI tag used: only used by mpiGet and mpiSend...
        self.shmname=shmname#shm file name used - only used by shmget and send.
        self.dictionary=dictionary
        if dictionary!=None:
            if dictionary.has_key("cpu"):
                self.cpu=eval(dictionary["cpu"])
            if dictionary.has_key("import"):
                self.imp=dictionary["import"]
            if dictionary.has_key("object"):
                self.object=dictionary["object"]
            if dictionary.has_key("pos"):
                self.pos=eval(dictionary["pos"])
            if dictionary.has_key("tag"):
                self.tag=int(dictionary["tag"])
            if dictionary.has_key("shortname"):
                self.shortname=dictionary["shortname"]
            if dictionary.has_key("pixmap"):
                self.pixmap=dictionary["pixmap"]
            if dictionary.has_key("feedback"):
                self.feedback=int(dictionary["feedback"])
            if dictionary.has_key("pyname"):
                self.pyname=dictionary["pyname"]
            if dictionary.has_key("groupshare"):
                self.groupShare=int(dictionary["groupshare"])
            if dictionary.has_key("args"):
                self.args=dictionary["args"]
            if dictionary.has_key("connectto"):
                self.connectto=eval(dictionary["connectto"])
            if dictionary.has_key("connectfrom"):
                self.connectfrom=eval(dictionary["connectfrom"])
            if dictionary.has_key("textcol"):
                self.textcol=dictionary["textcol"]
            if dictionary.has_key("idstr"):
                self.idstr=dictionary["idstr"]
        if self.idstr==None or self.idstr=="":
            try:
                d=eval(getArgsAsString(self.args))
                if d.has_key("idstr"):
                    self.idstr=d["idstr"]
            except:
                pass
    def __repr__(self):
        txt="%s args=%s, cpu=%s connectto=%s connectfrom=%s shareto=%s parentNames=%s tag=%s idstr=%s\n"%(self.object,str(self.args),str(self.cpu),str(self.connectto),str(self.connectfrom),str(self.sharedTo),str(self.parentNames),str(self.tag),str(self.idstr))
        return txt
    
    def getConnectionsFromTxt(self):
        """use linestxt etc to fill lines..."""
        print "TODO"
    def finalise(self):
        self.parentNames=eval(self.parentNamesTxt)
        if type(self.parentNames)!=type([]):
            print "WARNING - parentNames not correct: %s"%self.parentNamesTxt
            self.parentNames=[]
        self.endlines=eval(self.endlinestxt)
        if type(self.endlines)!=type([]):
            self.endlines=[]
        if type(self.lines)!=type([]):
            self.lines=[]
        self.lines=eval(self.linestxt)
        self.postcode=self.postcodetxt.replace("\t","    ")
        self.postcode="    "+self.postcode.replace("\n","\n    ").strip()#add extra indent
        if len(self.postcode.strip())>0:
            self.postcode+="\n"
        else:
            self.postcode=""
        self.precode=self.precodetxt.replace("\t","    ")
        self.precode="    "+self.precode.replace("\n","\n    ").strip()#add extra indent
        if len(self.precode.strip())>0:
            self.precode+="\n"
        else:
            self.precode=""
        self.sharedTo=eval(self.sharedToTxt)
def getArgsAsString(args):
    if args=="" or args==None:
        args="{}"
    else:
        args=str(args)
    try:
        tmp=eval(args)
    except:
        tmp=args
    if type(tmp)!=type({}):
        #args='{"idstr":"%s"}'%args
        args="%s"%args
    return args
