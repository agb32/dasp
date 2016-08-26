#!/usr/bin/env python
#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pygtk
pygtk.require("2.0")
import gtk, gobject
import gtk.glade as glade
import gui.dialog.dialog
import base.readConfig,os,types,sys,string,socket
import gui.selectionbox.selectionbox as selectionbox
from gui.myFileSelection.myFileSelection import myFileSelection
import base.dataType
#import Pyro.core,Pyro.naming

class paramGUI:#(Pyro.core.ObjBase):
    def __init__(self,root=""):
##         Pyro.core.ObjBase.__init__(self)
        try:
            gladefile=__file__.split("/")
            if gladefile[-2]=="bin":
                gladefile[-2]="paramgui"
            gladefile[-1]="paramgui.glade"
            gladefile=string.join(gladefile,"/")
        except:
            gladefile="paramgui.glade"
        self.gladetree=glade.XML(gladefile)
        self.sigdict={"on_new1_activate":self.on_new1_activate,
                      "on_open1_activate":self.on_open1_activate,
                      "on_save1_activate":self.on_save1_activate,
                      "on_save_as1_activate":self.on_save_as1_activate,
                      "on_quit1_activate":self.on_quit,
                      "on_duplicate1_activate":self.duplicate,
                      "on_copy1_activate":self.copy,
                      "on_paste1_activate":self.paste,
                      "on_add_from_file1_activate":self.on_add_from_file1_activate,
                      "on_skeleton_from_schema_file_activate":self.on_skel_from_schema_activate,
                      "on_delete2_activate":self.on_delete2_activate,
                      "on_about1_activate":self.on_about1_activate,
                      "on_buttonReCalculate_clicked":self.on_buttonReCalculate_clicked,
                      "on_buttonInsertModule_clicked":self.on_buttonInsertModule_clicked,
                      "on_buttonInsertVar_clicked":self.on_buttonInsertVar_clicked,
                      "on_buttonGetGlobals_clicked":self.getGlobals,
                      "on_window1_delete_event":self.on_quit,
                      "on_window1_destroy_event":self.on_quit,
                      "on_buttonQuit_clicked":self.on_quit,
                      "on_buttonDelete_clicked":self.on_buttondelete_activate,
                      "on_buttonMoveUp_clicked":self.on_buttonMoveUp_clicked,
                      "on_buttonMoveDown_clicked":self.on_buttonMoveDown_clicked,
                      "on_buttonShrink_clicked":self.shrinkAll,
                      "on_import_from_param_file1_activate":self.importModule,
                      "on_buttonSearch_clicked":self.on_buttonsearch_clicked,
                      "on_buttonDuplicate_clicked":self.duplicate,
                      "on_buttonCopy_clicked":self.copy,
                      "on_buttonPaste_clicked":self.paste
                      
                      }
        self.gladetree.signal_autoconnect(self.sigdict)
##         self.pyroSockets={}
##         self.pyroD=None
##         self.gotPyro=0
        sys.stdout=myStdOut(self.gladetree.get_widget("textviewStdout"))
        sys.stderr=sys.stdout
        self.aoxml=base.readConfig.AOXml(batchno=None)
        self.filename=None
        self.searchpos=0
        self.pasteData=None
        self.pymodfilename=None
        self.existfilename=None
        self.schemafilename=None
        self.treeView=self.gladetree.get_widget("treeview1")
        self.treeView.connect("key-press-event",self.keyPress,"treeView")
        self.treeStore=gtk.TreeStore(str,str,str,str,str,str,str)
        self.updateTree()
        self.host=socket.gethostname()
        self.port=8990
        self.clientList=[]
        self.openListenSocket()

    def openListenSocket(self):
        """so that simsetup can connect"""
        self.lsock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.lsock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
        bound=0
        while bound==0:
            try:
                self.lsock.bind((self.host,self.port))
            except:
                print "Couldn't bind to port %d, trying %d.  "%(self.port,self.port+1)
                self.port+=1#inc by number of processes in this mpi run
            else:
                bound=1
        self.lsock.listen(1)
        tag=gobject.io_add_watch(self.lsock,gobject.IO_IN,self.handleSocket)
        self.gladetree.get_widget("labelPortNumber").set_text("Port %d"%self.port)
    def handleSocket(self,sock=None,condition=None):
        print "handle socket"
        rtval=True
        if sock==self.lsock:
            print "accepting new client"
            conn,raddr=sock.accept()
            self.clientList.append(conn)
            tag=gobject.io_add_watch(conn,gobject.IO_IN,self.handleSocket)
        else:
            print "reading socket"
            if self.readsock(sock)==-1:
                sock.close()
                self.clientList.remove(sock)
                rtval=False
        return rtval
    
    def readsock(self,sock):
        rtval=0
        try:
            cmd=sock.recv(1024)
        except:
            rtval=-1
            print "error on socket"
        if len(cmd)==0:
            print "zero length cmd - closing socket"
            rtval=-1
        if rtval!=-1:
            d={"self":self}
            try:
                exec cmd in d,d
            except:
                print "remote command failed",cmd
        return rtval
    
##     def handlePyro(self,source=None,condition=None):
##         print "handlePyro"
##         if source!=None:
##             self.pyroD.handleRequests()
##         newPyroSockets=self.pyroD.getServerSockets()
##         for sock in newPyroSockets:
##             if sock not in self.pyroSockets.values():
##                 tag=gobject.io_add_watch(sock,gobject.IO_IN,self.handlePyro)
##                 self.pyroSockets[tag]=sock
##         for key in self.pyroSockets.keys():#now remove any old ones...
##             if self.pyroSockets[key] not in newPyroSockets:
##                 gobject.source_remove(key)
##                 del(self.pyroSockets[key])
##         return True

    def shrinkAll(self,w,a=None):
        self.updateTree()
    def updateTree(self):
        titleList=["Name","Value","Batch","Type","Comment","Calculated value","Other args"]
        ncols=len(titleList)
        #self.treeStore=gtk.TreeStore(str,str,str,str,str,str)
        self.treeStore.clear()
        #self.treeStoreData=[]#this is used to change when cells are edited.
        for module in self.aoxml.fileData.modules:
            pbatchno="All"
            comment=None
            otherargs={}
            for key in module.args.keys():
                if key=="batchno":
                    pbatchno=module.args[key]
                elif key=="comment":
                    comment=module.args[key]
                elif key in ["name"]:
                    pass
                else:
                    otherargs[key]=module.args[key]
##             if module.args.has_key("batchno"):
##                 pbatchno=module.args["batchno"]
##             else:
##                 pbatchno="All"
##             if module.args.has_key("comment"):
##                 comment=module.args["comment"]
##             else:
##                 comment=None
            parent=self.treeStore.append(None,(module.name,None,pbatchno,None,comment,None,str(otherargs)))
            self.showVariables(parent,module,pbatchno)
        self.treeView.set_model(self.treeStore)
        self.renderer={}
        for col in self.treeView.get_columns():
            self.treeView.remove_column(col)
        self.treeColumn={}
        for i in range(ncols):
            self.renderer[i]=gtk.CellRendererText()
            if i!=5:#calculated value isn't editable.
                self.renderer[i].set_property("editable",True)
            self.renderer[i].connect("edited",self.colEdited,i)
            self.treeColumn[i]=gtk.TreeViewColumn(titleList[i],self.renderer[i],text=i)
            self.treeColumn[i].set_resizable(True)
            self.treeView.append_column(self.treeColumn[i])

    def showVariables(self,parent,module,pbatchno=None):
        if pbatchno==None:
            if module.args.has_key("batchno"):
                pbatchno=module.args["batchno"]
            else:
                pbatchno="All"
            
        for vars in module.variableList:
            self.showVar(parent,module,vars,pbatchno)

    def showVar(self,parent,module,vars,pbatchno=None,pos=None):
        if pbatchno==None:
            if module.args.has_key("batchno"):
                pbatchno=module.args["batchno"]
            else:
                pbatchno="All"
        var=vars[0]
        batchno=pbatchno
        comment=None
        otherargs={}
        for key in var.keys():
            if key=="batchno":
               batchno=var[key]
            elif key=="comment":
                comment=var[key]
            elif key=="type":
                typ=var[key]
            elif key in ["name","value","extendedTxt"]:
                pass
            else:
                otherargs[key]=var[key]
##         if var.has_key("batchno"):
##             batchno=var["batchno"]
##         else:
##             batchno=pbatchno
##         if var.has_key("comment"):
##             comment=var["comment"]
##         else:
##             comment=None
##         if var.has_key("type"):
##             typ=var["type"]
##         else:
##             typ=None
                
        cvalstr=str(vars[1])
        if "\n" in cvalstr:
            cval="MULTILINE"
            xcval=cvalstr
        else:
            cval=cvalstr
            xcval=None
                    
        if "\n" in var["value"]:
            val="MULTILINE"
            xval=var["value"]
        else:
            val=var["value"]
            xval=None
        if pos==None:
            parent2=self.treeStore.append(parent,(var["name"],val,pbatchno,typ,comment,cval,str(otherargs)))
        else:
            parent2=self.treeStore.insert(parent,pos,(var["name"],val,pbatchno,typ,comment,cval,str(otherargs)))
            
        #self.treeStoreData.append(vars)
        if cval=="MULTILINE" or val=="MULTILINE":
            self.treeStore.append(parent2,(None,xval,None,None,None,xcval,None))
            #self.treeStoreData.append("MULTILINE")




    def colEdited(self,widget,rowno,new_text,column):
        #print "Edited row %s, col %s"%(rowno,column),new_text,type(widget)
        #note, rowno is the visible row number starting at 0, ignoring any compressed columns.  Need to get hold of which module/variable this is.
        #if self.treeStoreData[
        #cols=self.treeView.get_columns()
        print "old val",self.treeStore[rowno][column]
        rowinfo=map(int,rowno.split(":"))
        if rowinfo[0]==len(self.aoxml.fileData.modules):
            print "Empty module"
            return
        module=self.aoxml.fileData.modules[rowinfo[0]]
        if len(rowinfo)==1:#a module
            if column==0:#change module name...
                module.name=new_text
                module.args["name"]=new_text
                self.treeStore[rowno][0]=new_text
            elif column==2:#change batch...
                if new_text=="All" or new_text=="all" or new_text=="":
                    if module.args.has_key("batchno"):
                        del(module.args["batchno"])
                    self.treeStore[rowno][2]="All"
                else:#test its an allowed batch.
                    ok=1
                    try:
                        batch=eval(new_text)
                    except:
                        ok=0
                        print "Error - could not interpret batch"
                    else:
                        if type(batch)!=types.ListType and type(batch)!=types.TupleType:
                            batch=[batch]
                        for i in batch:
                            if type(i)!=types.IntType:
                                print "Error - batch type must be integer:",i
                                ok=0
                                break
                    if ok==1:
                        module.args["batchno"]=new_text
                        self.treeStore[rowno][2]=new_text
            elif column==4:#change a comment...
                if new_text=="":
                    if module.args.has_key("comment"):
                        del(module.args["comment"])
                    self.treeStore[rowno][4]=""
                else:
                    self.treeStore[rowno][4]=new_text
                    module.args["comment"]=new_text
            elif column==6:#change other args
                newargs=eval(new_text)
                for key in module.args.keys():
                    if key not in ["name","batchno","comment"]:
                        del(module.args[key])
                for key in newargs.keys():
                    module.args[key]=newargs[key]
                self.treeStore[rowno][6]=new_text
        elif len(rowinfo)==2:#a var
            var=module.variableList[rowinfo[1]][0]
            if column==0:#change name
                var["name"]=new_text
                self.treeStore[rowno][0]=new_text
            elif column==1:#change value
                if new_text!="MULTILINE":
                    var["value"]=new_text
                    self.treeStore[rowno][1]=new_text
            elif column==2:#change batch
                if new_text=="All" or new_text=="all" or new_text=="":
                    if var.has_key("batchno"):
                        del(var["batchno"])
                    self.treeStore[rowno][2]="All"
                else:#test its an allowed batch.
                    ok=1
                    try:
                        batch=eval(new_text)
                    except:
                        ok=0
                        print "Error - could not interpret batch"
                    else:
                        if type(batch)!=types.ListType and type(batch)!=types.TupleType:
                            batch=[batch]
                        for i in batch:
                            if type(i)!=types.IntType:
                                print "Error - batch type must be integer:",i
                                ok=0
                                break
                    if ok==1:
                        var["batchno"]=new_text
                        self.treeStore[rowno][2]=new_text
            elif column==3:#change type
                if new_text in ["i","f","eval","code","copy","list","string",""]:
                    var["type"]=new_text
                    self.treeStore[rowno][3]=new_text
            elif column==4:#change comment
                if new_text=="":
                    if var.has_key("comment"):
                        del(var["comment"])
                    self.treeStore[rowno][4]=""
                else:
                    self.treeStore[rowno][4]=new_text
                    var["comment"]=new_text
            elif column==6:#other args
                newargs=eval(new_text)
                for key in var.keys():
                    if key not in ["name","type","value","batchno","comment"]:
                        del(var[key])
                for key in newargs.keys():
                    var[key]=newargs[key]
                self.treeStore[rowno][6]=new_text
                
        elif len(rowinfo)==3:#a multi-line var
            var=module.variableList[rowinfo[1]][0]
            if column==1:#change the value...
                var["value"]=new_text
                self.treeStore[rowno][1]=new_text
        else:
            print "Warning - don't know what row this is!!! (shouldn't ever see this message)"

    def on_new1_activate(self,mitem):
        self.aoxml=base.readConfig.AOXml(batchno=None)
        self.filename=None
        self.updateTree()
        


    def on_open1_activate(self,mitem):
        fn=myFileSelection("Open an xml parameter file",self.filename).fname
        print "Got filename",fn
        if fn!=None:
            self.filename=fn
            print "Opening file %s"%fn
            self.aoxml.reset()
            self.aoxml.open(fn,ignoreError=1)
            if self.aoxml.error:
                if gui.dialog.dialog.myDialog(msg="Contains python errors - continue?").resp=="ok":
                    pass
                else:
                    self.aoxml.reset()
            self.updateTree()
    def on_save1_activate(self,mitem):
        print "save1_activated"
        if self.filename==None:
            self.on_save_as1_activate(mitem)
        else:
            if gui.dialog.dialog.myDialog(msg="Overwrite existing file?\n%s"%self.filename).resp=="ok":
                self.saveXML(self.filename)
    def on_save_as1_activate(self,mitem):
        fn=myFileSelection("Save XML parameter file as",self.filename).fname
        if fn!=None:
            self.saveXML(fn)
            self.filename=fn
            
    def on_quit(self,widget,arg2=None):
        gtk.main_quit()
    def on_cut1_activate(self,mitem):
        print "cut1_activated"
    def on_paste1_activate(self,mitem):
        print "paste1_activated"
    def on_delete1_activate(self,mitem):
        print "delete1_activated"
    def on_add_from_file1_activate(self,mitem,fn=None,comment=None,idstr=""):
        """Get the params needed by a given science module"""
        if fn==None:
            fn=myFileSelection("Open a simulation python module",self.pymodfilename).fname
            if fn!=None:
                self.pymodfilename=fn
        if fn!=None:
            if fn[-3:]=="pyc":
                fn=fn[:-1]
            varList=self.grepValues(fn)
            modname=fn.split("/")[-1][:-3]
            if len(idstr)>0:
                modname+="_"+idstr
            d={"name":modname}
            if comment!=None:
                d["comment"]=comment
            newmod=base.readConfig.ModuleClass(modname,d)
            self.aoxml.fileData.modules.append(newmod)
            parent=self.treeStore.append(None,(modname,None,"All",None,comment,None,"{}"))
            for var in varList:
                self.treeStore.append(parent,(var.description,var.val,"All",var.type,var.comment,"NOT EVALUATED","{}"))
                newmod.variableList.append(({"name":var.description,"value":var.val,"type":var.type,"comment":var.comment},None))
                
    def addBlankTemplate(self,modname,comment):
        d={"name":modname}
        if comment!=None:
            d["comment"]=comment
        newmod=base.readConfig.ModuleClass(modname,d)
        self.aoxml.fileData.modules.append(newmod)
        parent=self.treeStore.append(None,(modname,None,"All",None,comment,None,"{}"))
        

    def grepValues(self,filename):
        varList=[]
        module=None
        for p in sys.path:
            if len(p)>0 and filename.find(p)>-1:
                print p
                module=filename[len(p)+1:-3]
                m=module.split("/")
                module=string.join(m,".")
                print "Got module",module
                break
        lines=open(filename).readlines()
        for line in lines:
            if module!=None:
                classpos=line.find("class")
                bpos=line.find("(base.aobase.aobase)")
                if classpos>-1 and bpos>-1: #got a science object
                    objname=line[classpos+6:bpos].strip()
                    try:
                        d={"config":config()}
                        execstr="import %s;obj=%s.%s(None,config,forGUISetup=1)"%(module,module,objname)
                        print execstr
                        exec execstr in d,d
                        obj=d["obj"]
                        paramList=obj.getParams()
                    except:
                        print "Unable to access getParams for module %s"%module
                        paramList=[]
                    varList+=paramList
            
            pos=line.find("config.getVal(")
            hpos=line.find("#")
            ipos=line.find("backwards compatibility")
            if ipos<0:
                ipos=line.find("depreciated")
            if ipos<0:
                ipos=line.find("backwards compatible")
            if pos>=0 and (hpos<0 or hpos>pos):
                txt=line[pos+14:]
                #print txt
                pos=txt.find('"')
                if pos>=0:
                    txt=txt[pos+1:]
                pos=txt.find('"')
                if pos>=0:
                    val=txt[:pos]
                    found=0
                    for var in varList:
                        if var.description==val:
                            found=1
                            break
                    if found==0 and ipos<0:#not already included in list... and not a variable kept for backwards compatibility.
                        varList.append(base.dataType.dataType(description=val,typ="eval",val="None",comment=line.strip()))
                    #if val not in varList:
                    #    varList.append((txt[:pos],line.strip()))
                        #print "Got var:",varList[-1]
        #print varList
        return varList            

    def on_buttondelete_activate(self,w,arg2=None):
        """delete module or variable"""
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)>0:
            path=pathlist[0]
            if len(path)==1:#delete module
                del(self.aoxml.fileData.modules[path[0]])
                titer=self.treeStore.get_iter((path[0],))
                self.treeStore.remove(titer)
                self.treeView.get_selection().select_path((path[0],))
            elif len(path)>1:
                module=self.aoxml.fileData.modules[path[0]]
                del(module.variableList[path[1]])
                titer=self.treeStore.get_iter(path[0:2])
                self.treeStore.remove(titer)
                self.treeView.get_selection().select_path((path[0],path[1]))
        
    def on_delete2_activate(self,mitem):
        """Delete currently selected module"""
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)>0:
            path=pathlist[0]
            del(self.aoxml.fileData.modules[path[0]])
            titer=self.treeStore.get_iter((path[0],))
            self.treeStore.remove(titer)
            self.treeView.get_selection().select_path((path[0],))
    def on_about1_activate(self,mitem):
        print "AOSim Parameter GUI"
    def on_buttonReCalculate_clicked(self,button):
        tmpaoxml=base.readConfig.AOXml(batchno=None)
        self.aoxml.writeOut(".tmpaoxml.xml")
        tmpaoxml.reset()
        try:
            tmpaoxml.open(".tmpaoxml.xml")
        except:
            print "Error evaluating code"
            os.remove(".tmpaoxml.xml")
            raise
        else:
            self.aoxml=tmpaoxml
            os.remove(".tmpaoxml.xml")
        self.updateTree()
        
        
    def on_buttonInsertModule_clicked(self,button):
        """Insert a new module after current selection"""
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)==0:
            pathlist=[(0,)]
        path=pathlist[0]
        self.aoxml.fileData.modules.insert(path[0]+1,base.readConfig.ModuleClass("New",{"name":"New"}))
        self.treeStore.insert(None,path[0]+1,("New",None,"All",None,None,None,"{}"))
        self.treeView.get_selection().select_path((path[0]+1,))
    def on_buttonInsertVar_clicked(self,button):
        """Insert new variable after current selection"""
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)==0:
            pathlist=[(0,)]
        path=pathlist[0]
        module=self.aoxml.fileData.modules[path[0]]
        if len(path)>1:
            pos=path[1]+1
        else:
            pos=0#insert at beginning.
        module.variableList.insert(pos,({"name":"New","value":"None","type":"eval"},"NOT EVALUATED"))
        parent=self.treeStore.get_iter((path[0],))
        self.treeStore.insert(parent,pos,("New","None","All","eval","","NOT EVALUATED","{}"))
        self.treeView.get_selection().select_path((path[0],pos))
    def saveXML(self,filename):
        print "Saving",filename
        self.aoxml.writeOut(filename)

    def keyPress(self,widget,event,arg3=None):
        #print "keypress",type(event),arg3,dir(event)
        if widget==self.treeView:
            print "treeview",event.keyval
            if event.keyval in [65535,65288,65439]:#delete event
                self.on_buttondelete_activate(None)
            elif event.keyval==65379:#insert
                pathlist=[]
                self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
                if len(pathlist)==0:
                    pathlist=[(0,)]
                path=pathlist[0]
                if len(path)==1:
                    self.on_buttonInsertModule_clicked(None)
                else:
                    self.on_buttonInsertVar_clicked(None)
    def on_buttonMoveUp_clicked(self,button):
        self.moveObj("up")
    def on_buttonMoveDown_clicked(self,button):
        self.moveObj("down")
    def moveObj(self,dirn):
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)>0:
            path=pathlist[0]
            module=self.aoxml.fileData.modules[path[0]]
            pbatchno="All"
            comment=None
            otherargs={}
            for key in module.args.keys():
                if key=="batchno":
                    pbatchno=module.args[key]
                elif key=="comment":
                    comment=module.args[key]
                elif key in ["name"]:
                    pass
                else:
                    otherargs[key]=module.args[key]
##             if module.args.has_key("batchno"):
##                 pbatchno=module.args["batchno"]
##             else:
##                 pbatchno="All"
##             if module.args.has_key("comment"):
##                 comment=module.args["comment"]
##             else:
##                 comment=None
            if len(path)==1:#move module
                del(self.aoxml.fileData.modules[path[0]])
                if dirn=="up" and path[0]>0:
                    newpos=path[0]-1
                elif dirn=="down":
                    newpos=path[0]+1
                else:
                    newpos=path[0]
                self.aoxml.fileData.modules.insert(newpos,module)
                titer=self.treeStore.get_iter((path[0],))
                self.treeStore.remove(titer)
                parent=self.treeStore.insert(None,newpos,(module.name,None,pbatchno,None,comment,None,str(otherargs)))
                self.showVariables(parent,module)
                self.treeView.get_selection().select_path((newpos,))
                    
            else:#move variable
                newmodule=None
                if path[1]==0 and dirn=="up":#move to prev module
                    if path[0]>0:
                        newmodule=self.aoxml.fileData.modules[path[0]-1]
                        modpath=path[0]-1
                        newpos=len(newmodule.variableList)
                elif path[1]==len(module.variableList)-1 and dirn=="down":#move to next module
                    if path[0]<len(self.aoxml.fileData.modules)-1:
                        newmodule=self.aoxml.fileData.modules[path[0]+1]
                        modpath=path[0]+1
                        newpos=0
                else:
                    newmodule=module
                    if dirn=="up":
                        newpos=path[1]-1
                    elif dirn=="down":
                        newpos=path[1]+1
                    modpath=path[0]
                if newmodule!=None:
                    var=module.variableList[path[1]]
                    del(module.variableList[path[1]])
                    newmodule.variableList.insert(newpos,var)
                    titer=self.treeStore.get_iter((path[0],path[1]))
                    self.treeStore.remove(titer)
                    titer=self.treeStore.get_iter((modpath,))
                    value="None"
                    pbatchno="All"
                    comment=None
                    typ="eval"
                    otherargs={}
                    for key in var[0].keys():
                        if key=="batchno":
                            pbatchno=var[0][key]
                        elif key=="comment":
                            comment=var[0][key]
                        elif key in ["name"]:
                            pass
                        elif key=="value":
                            value=var[0][key]
                        elif key=="type":
                            typ=var[0][key]
                        else:
                            otherargs[key]=var[0][key]
                    if pbatchno=="All" and newmodule.args.has_key("batchno"):
                        pbatchno=newmodule.args["batchno"]
                    parent=self.treeStore.insert(titer,newpos,(var[0]["name"],value,pbatchno,typ,comment,"NOT EVALUATED",str(otherargs)))
                    self.treeView.get_selection().select_path((modpath,newpos))


    def importModule(self,w,arg2=None):
        """Import a single module from an existing parameter file"""
        fn=myFileSelection("Open an xml parameter file",self.existfilename).fname
        if fn!=None:
            self.existfilename=fn
            tmpaoxml=base.readConfig.AOXml(batchno=None)
            tmpaoxml.reset()
            try:
                tmpaoxml.open(fn)
            except:
                print "Error evaluating parameter file %s"%fn
            else:
                mlist=[]
                for module in tmpaoxml.fileData.modules:
                    mlist.append(module.name)
                s=selectionbox.selectionbox(mlist,txt="Choose a module",title="Choose a module",multiple=1)
                if s.choice!=None:
                    for selection in s.choice:#selection is integer...
                        module=tmpaoxml.fileData.modules[selection]
                        pbatchno="All"
                        comment=None
                        otherargs={}
                        for key in module.args.keys():
                            if key=="batchno":
                                pbatchno=module.args[key]
                            elif key=="comment":
                                comment=module.args[key]
                            elif key in ["name"]:
                                pass
                            else:
                                otherargs[key]=module.args[key]
    ##                     if module.args.has_key("batchno"):
    ##                         pbatchno=module.args["batchno"]
    ##                     else:
    ##                         pbatchno="All"
    ##                     if module.args.has_key("comment"):
    ##                         comment=module.args["comment"]
    ##                     else:
    ##                         comment=None
                        self.aoxml.fileData.modules.append(module)
                        parent=self.treeStore.append(None,(module.name,None,pbatchno,None,comment,None,str(otherargs)))
                        self.showVariables(parent,module)
                    self.updateTree()
    def copy(self,w,arg2=None):
        """copy module or variable into clipboard"""
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)>0:
            path=pathlist[0]
            module=self.aoxml.fileData.modules[path[0]].copy()
            if len(path)==1:#copy module
                self.pasteData=("module",module)
            elif len(path)>1:#copy variable
                self.pasteData=("var",module.variableList[path[1]])
    def paste(self,w,arg2=None):
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)>0 and self.pasteData!=None:
            path=pathlist[0]
            if self.pasteData[0]=="module":
                module=self.pasteData[1]
                pbatchno="All"
                comment=None
                otherargs={}
                for key in module.args.keys():
                    if key=="batchno":
                        pbatchno=module.args[key]
                    elif key=="comment":
                        comment=module.args[key]
                    elif key in ["name"]:
                        pass
                    else:
                        otherargs[key]=module.args[key]
                self.aoxml.fileData.modules.insert(path[0]+1,module)
                parent=self.treeStore.insert(None,path[0]+1,(module.name,None,pbatchno,None,comment,None,str(otherargs)))
                self.showVariables(parent,module)
                self.treeView.get_selection().select_path((path[0]+1,))
            else:
                var=self.pasteData[1]
                module=self.aoxml.fileData.modules[path[0]]
                newpos=0
                if len(path)>1:
                    newpos=path[1]+1
                module.variableList.insert(newpos,var)
                parent=self.treeStore.get_iter((path[0],))
                self.showVar(parent,module,var,pbatchno=None,pos=newpos)
                self.treeView.get_selection().select_path((path[0],newpos))
    def duplicate(self,w,arg2=None):
        """duplicate a variable or module"""
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        if len(pathlist)>0:
            path=pathlist[0]
            module=self.aoxml.fileData.modules[path[0]].copy()
            pbatchno="All"
            comment=None
            otherargs={}
            for key in module.args.keys():
                if key=="batchno":
                    pbatchno=module.args[key]
                elif key=="comment":
                    comment=module.args[key]
                elif key in ["name"]:
                    pass
                else:
                    otherargs[key]=module.args[key]
##             if module.args.has_key("batchno"):
##                 pbatchno=module.args["batchno"]
##             else:
##                 pbatchno="All"
##             if module.args.has_key("comment"):
##                 comment=module.args["comment"]
##             else:
##                 comment=None
            if len(path)==1:#copy module
                self.aoxml.fileData.modules.insert(path[0]+1,module)
                parent=self.treeStore.insert(None,path[0]+1,(module.name,None,pbatchno,None,comment,None,str(otherargs)))
                self.showVariables(parent,module)
                self.treeView.get_selection().select_path((path[0]+1,))
            elif len(path)>1:#copy variable
                vars=module.variableList[path[1]][0]
                #vars=(vars[0].copy(),vars[1])
                module.variableList.insert(path[1]+1,vars)
                parent=self.treeStore.get_iter((path[0],))
                self.showVar(parent,module,vars,pbatchno=None,pos=path[1]+1)
                self.treeView.get_selection().select_path((path[0],path[1]+1))

    def getGlobals(self,w,arg2=None):
        """Method to parse all modules in the tree, and place any common variables into the globals object..."""
        freqDict={}#this will be a dict of names, each of which will be a dict of types, each of which will be a dict of values.  If this final dict then contains more than one var, they can be put into globals...
        globmod=None
        for module in self.aoxml.fileData.modules:
            if module.name=="globals":
                globmod=module
            else:
                for var,calcvar in module.variableList:
                    name=None
                    typ="eval"
                    value=""
                    if var.has_key("name"):
                        name=var["name"]
                    if var.has_key("type"):
                        typ=var["type"]
                    if var.has_key("value"):
                        value=var["value"]
                    if name!=None and freqDict.has_key(name):
                        if freqDict[name].has_key(typ):
                            if freqDict[name][typ].has_key(value):
                                freqDict[name][typ][value].append((var,calcvar,module))
                            else:
                                freqDict[name][typ][value]=[(var,calcvar,module)]
                        else:
                            freqDict[name][typ]={value:[(var,calcvar,module)]}
                    elif name!=None:
                        freqDict[name]={typ:{value:[(var,calcvar,module)]}}
        if globmod==None:
            d={"name":"globals"}
            globmod=base.readConfig.ModuleClass("globals",d)
            self.aoxml.fileData.modules.insert(0,globmod)
        #Now if there is more than one identical entry in usagelist append this entry to globals, and remove it from the module
        for name in freqDict.keys():
            for typ in freqDict[name].keys():
                for value in freqDict[name][typ].keys():
                    vars=freqDict[name][typ][value]
                    globvar=None
                    for gv,gcv in globmod.variableList:
                        if gv.has_key("name") and name==gv["name"]:#variable already exists in globals...
                            globvar=gv
                            break
                    if globvar!=None:#variable already exists in globals...
                        if globvar.has_key("type") and typ==globvar["type"] and globvar.has_key("value") and value==globvar["value"]:
                            #same variable already in globals - remove from the other modules.
                            for v,c,m in vars:
                                m.variableList.remove((v,c))
                        else:
                            #has a different value in globals, so leave in other modules... ie do nothing
                            pass
                    else:#variable doesn't yet exist in globals...
                        if len(vars)>1:#can be put in globals if not there already.
                            newvar=vars[0][0]
                            calcvar=vars[0][1]
                            comment=""
                            for v,c,m in vars:
                                if v.has_key("comment"):
                                    comment+="%s from %s, "%(v["comment"],m.name)
                                else:
                                    comment+="from %s, "%m.name
                                m.variableList.remove((v,c))
                            newvar["comment"]=comment
                            globmod.variableList.append((newvar,calcvar))
        #may need to worry about the order in which the variables are put into globals, ie if they depend on each other etc... maybe this should be left to the user?
        #and now update the tree...
        self.updateTree()

    def findOrAdd(self,module,imp="science",idstr=None):
        """Called from simsetup gui - will search for a module, if find it
        will display there, if not, will add it.
        Module is the name, eg el_dm.  imp is where to import from.
        """
        print "findoradd",module,imp,idstr
        foundgen=self.search(module)#found the general version
        found=foundgen
        if foundgen==None:
            #add module, and associated variables
            modimport="%s.%s"%(imp,module)
            try:
                d={}
                exec "import %s;file=%s.__file__"%(imp,imp) in d,d
                modfile=d["file"]
                print "import worked for %s"%modfile
            except:#couldn't get filename... get from pythonpath
                for p in sys.path:
                    apos=p.find("/aosim")
                    if apos>-1:
                        modfile=p[:apos+6]+"/science/%s.py"%module
                        break
                print "import didn't work for",modfile
            print "modulename",module,modfile


            self.on_add_from_file1_activate(None,fn=modfile,comment="added from setup GUI")

        idstrlist=[]
        if idstr!=None and len(idstr)>0:
            if idstr[0]=="[":
                idstrlist=eval(idstr)
            else:
                idstrlist=[idstr]
        for id in idstrlist:
            if len(str(id))>0:
                modname=module+"_"+str(id)
            else:
                modname=module
            found=self.search(modname)#search for the specific version.
            if found==None:
                print "Adding module",module
                if len(str(id))>0:
                    #add empty module - also see whether we need to add a general module first...
                    print "adding empty module",modname
                    found=self.search(module)

                    self.addBlankTemplate(modname=modname,comment="added from setup GUI")
                #now get the position...
                found=self.search(modname)
        self.displayPosition(found)
        return True
    
    def on_buttonsearch_clicked(self,w,arg2=None):
        sterm=self.gladetree.get_widget("entrySearch").get_text()
        found=self.search(sterm,searchVars=1)
        self.displayPosition(found)
    def search(self,sterm,searchVars=0):
        poslist=[]
        pos=0
        modno=0
        varno=0
        found=None
        for module in self.aoxml.fileData.modules:
            if module.name==sterm:
                poslist.append((pos,modno,-1))
            pos+=1
            varno=0
            if searchVars:
                for var,calcvar in module.variableList:
                    if var.has_key("name") and var["name"]==sterm:
                        poslist.append((pos,modno,varno))
                    pos+=1
                    varno+=1
            modno+=1
        for p,m,v in poslist:
            if p>self.searchpos:
                found=(p,m,v)
                self.searchpos=p
                break
        if found==None:
            if len(poslist)>0:
                self.searchpos=poslist[0][0]
                found=poslist[0]
        return found

    def displayPosition(self,found):
        #now display this position (found).
        #print "Found: %s"%str(found)
        if found!=None:
            self.treeView.expand_row((found[1],),True)
            if found[2]>=0:
                self.treeView.expand_row((found[1],found[2]),True)
                self.treeView.get_selection().select_path((found[1],found[2]))
                self.treeView.scroll_to_cell((found[1],found[2]))
            else:
                self.treeView.get_selection().select_path((found[1],))
                self.treeView.scroll_to_cell((found[1],))
    def on_skel_from_schema_activate(self,w,arg2=None):
        fn=myFileSelection("Open an python schema file",self.schemafilename,complete="*.py").fname
        if fn!=None and fn[-3:]==".py":
            modulesUsed=[]
            self.schemafilename=fn
            lines=open(fn).readlines()
            for line in lines:
                pos=line.find("#")
                if pos>=0:
                    line=line[:pos]
                pos=line.find("science.")
                if line.find("import")==-1:#not an import line
                    instantiateTxt=None
                    modobj=None
                    args=None
                    if pos>=0:
                        instantiateTxt=line[pos+8:]
                        ppos=instantiateTxt.find("(")
                        if ppos>=0:
                            modobj=instantiateTxt[:ppos]
                            args=instantiateTxt[ppos:]
                    if instantiateTxt!=None and modobj!=None and args!=None:
                        #print modobj,args
                        targs=args[:]
                        objname=line[:pos-1]
                        modulename=modobj.split(".")[0]
                        modimport="science."+modulename
                        try:
                            d={}
                            exec "import %s;file=%s.__file__"%(modimport,modimport) in d,d
                            modfile=d["file"]
                            print "import worked for",modfile
                        except:#couldn't get filename... get from pythonpath
                            for p in sys.path:
                                apos=p.find("/aosim")
                                if apos>-1:
                                    modfile=p[:apos+6]+"/science/%s.py"%modulename
                                    break
                            print "import didn't work for",modfile
                        print "modulename",modulename,modfile
                        idstr=None
                        # args=line[ppos+1:]
                        #get the id string...
                        idstrlist=None
                        #look for an idstr in the args dictionary...
                        argpos=args.find("args={")
                        if argpos>-1:
                            args=args[argpos+5:]
                            endargpos=args.find("}")
                            if endargpos>-1:
                                args=args[:endargpos+1]
                                idpos=args.find('"idstr":')
                                if idpos>-1:
                                    id=args[idpos+8:]
                                    if id[0]=="[":#a list...
                                        end=id.find("]")
                                        if end>0:
                                            id=id[:end+1]
                                            idstrlist=eval(id)
                                    else:#a string or int etc...
                                        eid=id.find(",")
                                        eeid=id.find("}")
                                        if eid>-1:
                                            if eeid>-1:
                                                if eid<eeid:
                                                    end=eid
                                                else:
                                                    end=eeid
                                            else:
                                                end=eid
                                        else:
                                            end=eeid
                                        idstr=id[:end].strip().strip('"')
                                    print idstr,idstrlist
    ##                             try:
    ##                                 args=eval(args)
    ##                                 if args.has_key("idstr"):
    ##                                     idstr=args["idstr"]
    ##                             except:#unable to get idstr - probably an evaluated string.
    ##                                 idstr="CHANGE"
                        idpos=args.find("idstr=")
                        if idpos>-1:
                            idtxt=args[idpos+6:].strip()
                            if idtxt[0]=='[':#a list of idstrs
                                endpos=idtxt.find(']')
                                if endpos>0:
                                    idtxt=idtxt[:end+1]
                                idstrlist=eval(idtxt)
                            else:#a string or eg int/float...
                                endposa=idtxt.find(",")
                                endposb=idtxt.find(")")
                                if endposa>0:
                                    if endposb>0:
                                        if endposa<endposb:
                                            end=endposa
                                        else:
                                            end=endposb
                                    else:
                                        end=endposa
                                else:
                                    end=endposb
                                idtxt=idtxt[:end].strip().strip('"').strip("'")
                                idstr=idtxt
                            print idstr,idstrlist
                        #now look for a parameter of the form idstr="xxx".
                        idpos=targs.find("idstr=")
                        if idpos>=-1:
                            idtxt=targs[idpos+6:].strip()
                            if idtxt[0]=="[":
                                idpos=idtxt.find("]")
                                if idpos>=-1:
                                    idtxt=idtxt[:idpos+1]
                                    idstrlist=eval(idtxt)
                            elif idtxt[0]=="'":
                                idpos=idtxt[1:].find("'")
                                if idpos>=-1:
                                    idtxt=idtxt[1:idpos+1]
                                    idstr=idtxt
                            elif idtxt[0]=='"':
                                idpos=idtxt[1:].find('"')
                                if idpos>=-1:
                                    idtxt=idtxt[1:idpos+1]
                                    idstr=idtxt
                                print "GOT idstr",idstr,idstrlist
                                    
                        if modulename not in modulesUsed:
                            modulesUsed.append(modulename)
                            comment=None
                            if idstr==None and idstrlist==None:
                                comment=objname.strip()
                            self.on_add_from_file1_activate(None,fn=modfile,comment=comment)
                        if idstr!=None:
                            modulename2=modulename+"_"+idstr
                            if modulename2 not in modulesUsed:
                                modulesUsed.append(modulename2)
                                self.addBlankTemplate(modname=modulename2,comment="%s (add variables specific to this object)"%objname.strip())
                        if idstrlist!=None:
                            for idstr in idstrlist:
                                modulename2=modulename+"_"+str(idstr)
                                if modulename2 not in modulesUsed:
                                    modulesUsed.append(modulename2)
                                    self.addBlankTemplate(modname=modulename2,comment="%s (add variables specific to this object)"%objname.strip())
                                
                                #self.on_add_from_file1_activate(None,fn=modfile,comment="%s (delete variables shared between all %s objects)"%(objname.strip(),modulename),idstr=idstr)

##         if fn!=None:
##             modulesUsed=[]
##             self.schemafilename=fn
##             lines=open(fn).readlines()
##             for line in lines:
##                 pos=line.find("#")
##                 if pos>=0:
##                     line=line[:pos]
##                 pos=line.find("science.")
##                 eqpos=line.find("=")
##                 ppos=line.find("(")
##                 if pos>=0 and line.find("import")==-1 and eqpos>=0 and ppos>=0:
##                     objname=line[:eqpos]
##                     modulename=line[pos+8:ppos].split(".")[0]
##                     modimport="science."+modulename
##                     try:
##                         d={}
##                         exec "import %s;file=%s.__file__"%(modimport,modimport) in d,d
##                         modfile=d["file"]
##                     except:#couldn't get filename... get from pythonpath
##                         for p in sys.path:
##                             apos=p.find("/aosim")
##                             if apos>-1:
##                                 modfile=p[:apos+6]+"/science/%s.py"%modulename
##                                 break
##                     idstr=None
##                     args=line[ppos+1:]
##                     argpos=args.find("args={")
##                     if argpos>-1:
##                         args=args[argpos+5:]
##                         endargpos=args.find("}")
##                         if endargpos>-1:
##                             args=args[:endargpos+1]
##                             try:
##                                 args=eval(args)
##                                 if args.has_key("idstr"):
##                                     idstr=args["idstr"]
##                             except:#unable to get idstr - probably an evaluated string.
##                                 idstr="CHANGE"
##                     if modulename not in modulesUsed:
##                         modulesUsed.append(modulename)
##                         comment=None
##                         if idstr==None:
##                             comment=objname.strip()
##                         self.on_add_from_file1_activate(None,fn=modfile,comment=comment)
##                     if idstr!=None:
##                         modulename2=modulename+"_"+idstr
##                         if modulename2 not in modulesUsed:
##                             modulesUsed.append(modulename2)
##                             self.on_add_from_file1_activate(None,fn=modfile,comment="%s (delete variables shared between all %s objects)"%(objname.strip(),modulename),idstr=idstr)
        elif fn!=None and fn[-4:]==".xml":#create from a simsetup xml file.
            print "TODO: implement setup from simsetup xml file (currently only works for a .py file)"

class myStdOut:
    def __init__(self,textView,logfilename="paramgui.log"):
        self.textView=textView
        if self.textView!=None:
            self.buffer=self.textView.get_buffer()
        else:
            self.buffer=None
        if logfilename!=None:
            self.logfile=open(logfilename,"w")
        else:
            self.logfile=None

    def write(self,txt):
        if self.logfile!=None:
            self.logfile.write(txt)
            self.logfile.flush()
        if self.textView!=None:
            self.buffer.insert(self.buffer.get_end_iter(),txt)
            self.textView.scroll_to_mark(self.buffer.get_insert(), 0.4)
            #self.textView.scroll_to_iter(self.buffer.get_end_iter(),0.5)
        
class config:
    def getVal(self,val,default=None):
        return default
    def setSearchOrder(self,so):
        pass

## def connectPyro(paramgui):
##     ns=None
##     id=""
##     print "Connecting to pyro nameserver..."
##     try:
##         ns=Pyro.naming.NameServerLocator().getNS()
##         gotPyro=1
##     except:
##         print "Couldn't connect to pyro nameserver (pyrons)"
##         gotPyro=0
##     if gotPyro:
##         daemon=Pyro.core.Daemon()
##         daemon.useNameServer(ns)
##         id="paramgui_%s_%d"%(os.environ["USER"],os.getpid())
##         uri=daemon.connect(paramgui,id)
##         paramgui.gotPyro=1
##         paramgui.pyroD=daemon
##         paramgui.handlePyro()
##     return ns,id

def run():
    #Pyro.core.initServer(0)
    g=paramGUI()
    #ns,id=connectPyro(g)
    gtk.main()
    #if ns!=None:
    #    ns.unregister(id)
        
if __name__=="__main__":
    run()
