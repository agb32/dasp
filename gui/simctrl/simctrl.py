#!/usr/bin/env python

import pygtk
pygtk.require("2.0")
import gtk, gobject
import gtk.glade as glade
import string,socket,types
import util.serialise as serialise
import base.readConfig
import util.analyse
import sys,os
import gui.dialog.dialog
import popen2,signal
import gui.simctrl.simdata as simdata
from gui.myFileSelection.myFileSelection import myFileSelection
import gui.selectionbox.selectionbox as selectionbox
class simctrlGUI:
    """Class for simulation control GUI"""
    def __init__(self,root=""):
        try:
            gladefile=__file__.split("/")
            if gladefile[-2]=="bin":
                gladefile[-2]="simctrl"
            gladefile[-1]="simctrl.glade"
            gladefile=string.join(gladefile,"/")
        except:#older versions of python fail here...
            gladefile="simctrl.glade"
        self.gladetree=glade.XML(gladefile)
        self.sigdict={"on_buttonQuit_clicked":self.on_quit,
                      "on_buttonLoad_clicked":self.loadGUIFile,
                      "on_buttonFromSim_clicked":self.loadFromSim,
                      "on_buttonGetPlots_clicked":self.loadFromSimRestricted,
                      "on_buttonClearPlots_clicked":self.clearPlotControl,
                      "on_buttonStop_clicked":self.stop,
                      "on_buttonKill_clicked":self.kill,
                      "on_buttonPause_clicked":self.pause,
                      "on_buttonRun_clicked":self.run,
                      "on_buttonConnect_clicked":self.connect,
                      "on_buttonDisconnect_clicked":self.disconnect,
                      "on_buttonExecute_clicked":self.execute,
                      "on_cmd1_activate":self.set_cmd,
                      "on_rpt1_activate":self.set_rpt,
                      "on_now1_activate":self.set_now,
                      "on_del1_activate":self.set_del,
                      "on_buttonModules_clicked":self.listModules,
                      "on_buttonList_clicked":self.listSim,
                      "on_buttonClear_clicked":self.clearText,
                      "on_buttonSpawn_clicked":self.spawn,
                      "on_buttonHelp_clicked":self.help,
                      "on_togglebuttonPlot_toggled":self.plotDataToggled,
                      "on_window1_delete_event":self.on_quit,
                      "on_entryCommand_button_press_event":self.doCommandWindow,
                      "on_entryCommand_key_press_event":self.doCommandWindowKey,
                      "on_entryCommand2_key_press_event":self.doCommandWindowKey,
                      "on_entryReturn_key_press_event":self.doCommandWindowKey,
                      "on_entryReturn2_key_press_event":self.doCommandWindowKey,
                      "on_entryCommand_activate":self.execute,
                      "on_entryCommand2_activate":self.preexecute,
                      "on_entryReturn_activate":self.execute,
                      "on_entryReturn2_activate":self.preexecute,
                      "on_window2_delete_event":self.hideShowWidget,
                      "on_buttonDisplayAll_clicked":self.toggleDisplayAll,
                      "on_treeviewCommands_cursor_changed":self.treeSelect,
                      }
        self.sockidDict={}
        self.displayAll=0
        self.sim=util.analyse.analyse()
        self.simdata=None
        self.actionType="cmd"
        self.commandHistoryList=[]
        self.commandHistoryListPos=0
        self.gladetree.signal_autoconnect(self.sigdict)
        self.gladetree.get_widget("entryHostPort").set_text("9000")

        #set up the socket connection list
        self.sockView=self.gladetree.get_widget("treeviewConnections")
        self.sockView.get_selection().set_mode(gtk.SELECTION_MULTIPLE)
        self.sockListStore=gtk.ListStore(str,int)
        self.sockView.set_model(self.sockListStore)
        self.oldfilename=None
        title=["Host","Port"]
        for i in range(2):
            r=gtk.CellRendererText()
            t=gtk.TreeViewColumn(title[i],r,text=i)
            self.sockView.append_column(t)
        self.updateSockList()
        #set up the data plotting list
        dataTitles=["Title","Command","Return","Type","Preprocess","Post","Dimensions","Xaxis","When","Options"]
        self.dataView=self.gladetree.get_widget("treeviewCommands")
        self.dataTreeStore=gtk.TreeStore(gobject.TYPE_BOOLEAN,str,str,str,str,str,str,str,str,str,str)
        self.lasttreeiter=None
        self.dataView.set_model(self.dataTreeStore)
        r=gtk.CellRendererToggle()
        r.set_property("activatable",True)
        r.connect("toggled",self.plotToggled)
        c=gtk.TreeViewColumn("Button",r)
        c.add_attribute(r,"active",0)
        self.dataView.append_column(c)
        self.plotUserData=0
        i=1
        for i in range(len(dataTitles)):
            title=dataTitles[i]
            r=gtk.CellRendererText()
            r.set_property("editable",True)
            r.connect("edited",self.colEdited,i+1)
            c=gtk.TreeViewColumn(title,r,text=i+1)
            c.set_resizable(True)
            self.dataView.append_column(c)
        sys.stdout=myStdOut(self.gladetree.get_widget("textviewStdout"))
        sys.stderr=sys.stdout

    def plotDataToggled(self,w,a=None):
        """Should we plot data returned from the users command?"""
        self.plotUserData=w.get_active()

            
    def help(self,w,arg2=None):
        """Pseudo help function"""
        print "Help?  Try the pdf file in the docs directory"

    def spawn(self,w,arg2=None):
        print "Attempting to run simctrl.py (check that it can be found in your path)"
        os.system("simctrl.py&")
    def doCommandWindow(self,w,e,d=None):
        if e.button==3:
            self.hideShowWidget(self.gladetree.get_widget("window2"))
            return True
    def doCommandWindowKey(self,w,e,d=None):
        if e.keyval==65362:#up pressed
            if len(self.commandHistoryList)>0:
                self.commandHistoryListPos-=1
                if self.commandHistoryListPos<0:
                    self.commandHistoryListPos=0
                self.redisplayCommand(self.commandHistoryListPos)
            return True
        elif e.keyval==65364:#down pressed
            l=len(self.commandHistoryList)
            if l>0:
                self.commandHistoryListPos+=1
                if self.commandHistoryListPos>=l:
                    self.commandHistoryListPos=l-1
                self.redisplayCommand(self.commandHistoryListPos)
            return True
    def redisplayCommand(self,pos):
        cmd,ret,freq,self.actionType,self.plotUserData=self.commandHistoryList[pos]
        self.gladetree.get_widget("entryCommand").set_text(cmd)
        self.gladetree.get_widget("entryCommand2").set_text(cmd)
        self.gladetree.get_widget("entryReturn").set_text(ret)
        self.gladetree.get_widget("entryReturn2").set_text(ret)
        self.gladetree.get_widget("spinbuttonFreq").set_value(int(freq))
        self.gladetree.get_widget("togglebuttonPlot").set_active(self.plotUserData)
        m=["cmd","rpt","now","del"].index(self.actionType[:3])
        self.gladetree.get_widget("optionmenu1").set_history(m)
        tv=self.gladetree.get_widget("textviewCommand")
        tv.scroll_to_iter(tv.get_buffer().get_iter_at_line(pos),0.)
    def hideShowWidget(self,w,d=None):
        fw=w#.get_parent().get_children()[0]
        if fw.flags()&gtk.VISIBLE:
            fw.hide()
        else:
            fw.show_all()
        return True
    def preexecute(self,w,e=None):
        self.gladetree.get_widget("entryCommand").set_text(self.gladetree.get_widget("entryCommand2").get_text())
        self.gladetree.get_widget("entryReturn").set_text(self.gladetree.get_widget("entryReturn2").get_text())
        self.execute(w)
        
    def on_quit(self,widget,arg2=None):
        """Quit the gui"""
        if gui.dialog.dialog.myDialog(msg="Quit?\nAre you sure?\n",title="Really quit?!?",parent=self.gladetree.get_widget("window1")).resp=="ok":
            self.sim.closeConnection()
            gtk.main_quit()
        return True
    def stop(self,widget,arg2=None):
        """Stop simulation pressed"""
        #print "Stopping simulation"
        if gui.dialog.dialog.myDialog(msg="Stop the simulation?\nAre you sure?\n",title="Stop the simulation?!?",parent=self.gladetree.get_widget("window1")).resp=="ok":
            selectedConns=self.getConns()
            self.sim.stop(selectedConns)
    def kill(self,widget,arg2=None):
        """Kill simulation pressed"""
        #print "Stopping simulation"
        resp=gui.dialog.dialog.myDialog(msg="Kill all mpipython processes?\nAre you sure?\n",title="Kill the simulation?!?",parent=self.gladetree.get_widget("window1"),buttons=(gtk.STOCK_OK,gtk.RESPONSE_ACCEPT,gtk.STOCK_CANCEL,gtk.RESPONSE_REJECT,"All MPI jobs",-4)).resp
        if resp=="ok":
            selectedConns=self.getConns()
            self.sim.kill(selectedConns)
        elif resp=="All MPI jobs":
            self.forceMPIKill()
    def forceMPIKill(self):
        o,i=popen2.popen2("ps axu | grep %s"%os.environ["USER"])
        lines=o.readlines()
        print "Killing all MPI processes for user %s"%os.environ["USER"]
        for line in lines:
            if "mpirun" in line or "mpich" in line or "MPIRUN" in line or "mpipython" in line:
                pid=int(line.split()[1])
                print "Killing PID %d"%pid
                os.kill(pid,signal.SIGKILL)
       
        
    def pause(self,w,arg2=None):
        """Pause simulation pressed"""
        #print "Pausing simulation"
        selectedConns=self.getConns()
        self.sim.pause(selectedConns)
    def run(self,w,arg2=None):
        """run simulation pressed"""
        niter=self.gladetree.get_widget("spinbuttonRun").get_value_as_int()
        #print niter,type(niter)
        if niter==0:
            niter=None
        selectedConns=self.getConns()
        self.sim.run(niter,selectedConns)

    def connect(self,w,arg2=None):
        """connect pressed"""
        txt=self.gladetree.get_widget("entryHostPort").get_text()
        notowned=0
        if len(txt)>0:
            #already selected port and hostname.
            txt=txt.split()
            port=9000
            host="localhost"
            if len(txt)>1:
                host=txt[0]
                port=int(txt[1])
            elif len(txt)==1:
                port=int(txt[0])
            owned=self.makeConnection(host,port)
            if owned==-1:#not owned by you...
                notowned=1
            if owned==0:#not connected
                txt=""
        if len(txt)==0:
            #print "Trying default simulation port (9000)"
            owned=0#self.makeConnection("localhost",9000)
            if owned==-1:#not owned by you...
                notowned=1
            elif owned==0:#not connected
                try:

                    data=self.sim.queryPortDict()
                except:
                    data=[]
                if data==[]:
                    data=[["localhost",9000,"Guess",os.environ.get("USER","?")]]
                    # try:
                    #     data=self.sim.queryPortDict("192.168.3.30")
                    # except:
                    #     data=[]
                if len(data)==0:
                    gui.dialog.dialog.myDialog(msg="No simulation found to connect to.\nPlease check it is running.\nIf util/portdict.py is not running, try running it\nand restarting your simulation so that it is registered.\n",title="No simulations found",buttons=(gtk.STOCK_OK,gtk.RESPONSE_ACCEPT),parent=self.gladetree.get_widget("window1"))
                else:
                    mysims=[]
                    simulations=[]
                    mydata=[]
                    odata=[]
                    for d in data:
                        if d[3]==os.environ.get("USER","?"):
                            mysims.append(str(d))
                            mydata.append(d)
                        else:
                            simulations.append(str(d))
                            odata.append(d)
                    simulations=mysims+simulations
                    data=mydata+odata
                    s=selectionbox.selectionbox(simulations,txt="Choose the connections",title="Choose the connections",multiple=1,parent=self.gladetree.get_widget("window1"))
                    if s.choice!=None:
                        if type(s.choice)==types.ListType and len(s.choice)==0:
                            for d in data:#connect to all...
                                if self.makeConnection(d[0],d[1],d[3])==-1:
                                    notowned=1
                        else:
                            for selection in s.choice:
                                if self.makeConnection(data[selection][0],data[selection][1],data[selection][3])==-1:
                                    notowned=1
                        
        if notowned:
            gui.dialog.dialog.myDialog(msg="WARNING: Connected to simulation(s) not owned by you.\nPlease check that this is what you expected\n",title="Owner warning",buttons=(gtk.STOCK_OK,gtk.RESPONSE_ACCEPT),parent=self.gladetree.get_widget("window1"))
                
        self.updateSockList()
    def makeConnection(self,host,port,user=None):
        owned=1
        print "Connecting to",host,port
        if user!=None and user!=os.environ["USER"]:
            print "WARNING - not your simulation"
            owned=-1
        s=self.sim.openConnection(host,port)
        if s==None:
            print "Couldn't connect"
            owned=0
        else:
            #self.sockidDict[s]=gtk.input_add(s,gtk.gdk.INPUT_READ,self.handleData)
            self.sockidDict[s]=gobject.io_add_watch(s,gtk.gdk.INPUT_READ,self.handleData)
            tag=self.sim.execute("import os;d=(os.environ['USER'],ctrl.rank,'%s',ctrl.simID,__file__)"%os.environ["USER"],rt='d',action="now",connList=[s])
            self.sim.addCallback(s,tag,self.printWelcome)
        return owned

        
    def disconnect(self,w,arg2=None):
        """disconnect pressed"""
        selectedConns=self.getConns()
        for conn in selectedConns:
            gtk.input_remove(self.sockidDict[conn])
            del(self.sockidDict[conn])
        self.sim.closeConnection(selectedConns)
        self.updateSockList()
    def printWelcome(self,data):
        """Handle welcome message when connecting"""
        if data[1]=="warning" or data[1]=="error":
            print "Error retrieving connection data",data[1:],data[0]
        else:
            simid=""
            if len(data[3]["d"][3])>0:
                simid=" ID %s "%data[3]["d"][3]
            simfile=data[3]["d"][4]
            print "Connected to user %s's simulation %s, rank %d %s%s"%(data[3]["d"][0],simfile,data[3]["d"][1],simid,data[0].getpeername())
    def handleData(self,source,condition):
        """Handle data from simulation"""
        #source is a fd
        #print "Handling %s"%str(source)
        l=len(self.sim.connList)
        #print "len %d"%l
        #print self.sim.connList
        handledList=self.sim.process()
        if l!=len(self.sim.connList):#connections have been closed...
            #print "updateSockList"
            #print self.sim.connList
            for key in self.sockidDict.keys():
                if key not in self.sim.connList:#connection has been closed
                    print "Connection closed remotely"
                    gtk.input_remove(self.sockidDict[key])
                    del(self.sockidDict[key])
            self.updateSockList()
        #self.readGUIFeedback(handledList)
        for data in self.sim.recDataList:#unhandled things...
            for d in data:
                if type(d)==types.DictType:
                    for key in d.keys():
                        print key,":\n",d[key]
                elif type(d)==socket.SocketType:
                    print "Socket:",d.getpeername()
                else:
                    print d
            #print data
        self.sim.recDataList=[]
            
##         data=serialise.ReadMessage(self.sock.fileno())
##         print "Got data",type(data)
##         if type(data)==types.NoneType:
##             print "Disconnecting"
##             self.gladetree.get_widget("togglebuttonConnect").set_active(False)
        return True
    def clearPlotControl(self,w,arg2=None):
        """Clear the plot control pane"""
        if self.simdata!=None:
            #delete any previous ones
            self.sim.execute("",action="del",tag=None)
            #for plot in self.simdata.dataList:
            #    plot.cancel()
            self.cancelPlots()
            self.simdata=None
            self.updateSimDataList()

    def cancelPlots(self,plot=None):
        #for plot in self.simdata.dataList:
        #    plot.cancel()
        if plot==None:
            plot=self.simdata.dataList
        if plot!=None:
            plot.cancel()
            for p in plot.childList:
                self.cancelPlots(p)

    def loadFromSim(self,w,arg2=None):
        """Load GUI control functions from simulation"""
        connList=[self.getConns()[0]]
        tag=self.sim.execute("d=ctrl.simctrlXML",rt='d',action="cmd",connList=connList)
        self.sim.addCallback(connList[0],tag,self.loadFromSimXML)
    def loadFromSimRestricted(self,w,arg2=None):
        """Load GUI control functions from simulation"""
        connList=self.getConns()
        for conn in connList:
            tag=self.sim.execute("d=ctrl.simctrlXMLRestricted",rt='d',action="cmd",connList=[conn])
            self.sim.addCallback(conn,tag,self.loadFromSimXMLNoOverwrite)
        
    def loadFromSimXMLNoOverwrite(self,data):
        """handle loading of gui info from simulation without overwriting existing."""
        self.loadFromSimXML(data,cancel=0)
        
    def loadFromSimXML(self,data,cancel=1):
        """Handle the loading of GUI info from simulation..."""
        if data[1]=="warning" or data[1]=="error":
            print "Error retrieving simulation XML data",data[1:],data[0]
        else:
            simid=""
            xmltxt=data[3]["d"]
            if len(xmltxt)>0:
                open(".simctrltmp.xml","w").write(xmltxt)
                if self.simdata!=None and cancel==1:
                    #delete any previous ones.
                    self.sim.execute("",action="del",tag=None)
                    #for plot in self.simdata.dataList:
                    #    plot.cancel()
                    self.cancelPlots()
                if self.simdata!=None and cancel==0:
                    print "Receiving plot info"
                    #add new to old...
                    tmp=simdata.simdata(None,xmltxt,cancelData=self.sim.execute)
                    self.simdata.dataList.childList+=tmp.dataList.childList
                else:#overwrite old.
                    self.simdata=simdata.simdata(None,xmltxt,cancelData=self.sim.execute)
                self.updateSimDataList()
    def loadGUIFile(self,w,arg2=None):
        """Load a file for GUI control functions"""
        fn=myFileSelection("Open an xml parameter file",self.oldfilename,parent=self.gladetree.get_widget("window1")).fname
        
        if fn!=None:
            self.oldfilename=fn
            if self.simdata!=None:
                self.sim.execute("",action="del",tag=None)
                #for plot in self.simdata.dataList:
                #    plot.cancel()
                self.cancelPlots()
##                     if plot.repeating:#stop it repeating...
##                         self.sim.execute(plot.cmd,action="del")
            self.simdata=simdata.simdata(fn,cancelData=self.sim.execute)
            self.updateSimDataList()
    def updateSimDataList(self):
        """update the GUI functionality list"""
        self.dataTreeStore.clear()
        self.lasttreeiter=None
        if self.simdata!=None:
            for node in self.simdata.dataList.childList:
                self.addNode(node,None)
            #for plot in self.simdata.dataList:
            #    self.dataTreeStore.append((None,plot.title,plot.cmd,plot.ret,plot.gettype(),plot.preprocess,plot.post,str(plot.dim),plot.xaxis,plot.when,str(plot.info)))
    def addNode(self,node,parent):
        #print "addnode",node.title
        tu=[None]+node.aslist(self.displayAll)
        newparent=self.dataTreeStore.append(parent,tu)
        for child in node.childList:
            self.addNode(child,newparent)
        
    def treeSelect(self,w,arg2=None):
        """Called when a new row of the data list is selected.
        """
        treeiter=w.get_selection().get_selected()[1]
        if treeiter==None:
            return
        path=self.dataTreeStore.get_path(treeiter)
        plot=self.simdata.dataList
        for i in path:
            plot=plot.childList[i]
        tu=plot.aslist(1)
        for i in range(2,self.dataTreeStore.get_n_columns()):#len(self.dataTitles)):
            self.dataTreeStore.set(treeiter,i,tu[i-1])
        if self.lasttreeiter!=None:
            plot=self.simdata.dataList
            path=self.dataTreeStore.get_path(self.lasttreeiter)
            for i in path:
                plot=plot.childList[i]
            tu=plot.aslist(self.displayAll)
            for i in range(2,self.dataTreeStore.get_n_columns()):#len(self.dataTitles)):
                self.dataTreeStore.set(self.lasttreeiter,i,tu[i-1])
        self.lasttreeiter=treeiter
    def toggleDisplayAll(self,w,arg2=None):
        self.displayAll=1-self.displayAll
        self.updateSimDataDisplay()
    def updateSimDataDisplay(self):
        pathList=self.makePathList()
        for path in pathList:
            plot=self.simdata.dataList
            for i in path:
                plot=plot.childList[i]
            tu=plot.aslist(self.displayAll)
            for j in range(2,self.dataTreeStore.get_n_columns()):#len(self.dataTitles)):
                iter=self.dataTreeStore.get_iter(path)
                self.dataTreeStore.set(iter,j,tu[j-1])


    def makePathList(self):
        iterlist=[]
        def func(mode,path,iter,data):
            data.append(path)
        self.dataTreeStore.foreach(func,iterlist)
        return iterlist
    
    def colEdited(self,widget,path,new_text,column):
        """column has been edited in functionality list"""
        #print "Edited row %s, col %s"%(path,column),new_text,type(widget)
        if not self.dataTreeStore[path][0]:#not currently running - so can edit
            plot=self.simdata.dataList#[int(path)]
            epath=eval("["+path.replace(":",",")+"]")
            for i in epath:
                plot=plot.childList[i]
            self.dataTreeStore[path][column]=new_text
            if column==1:
                plot.title=new_text
            elif column==2:
                plot.command=new_text
            elif column==3:
                plot.ret=new_text
            elif column==4:
                ss=new_text.split()
                plot.filetype=plot.texttype=plot.pylabtype=plot.gisttype=plot.feedbacktype=0
                for s in ss:
                    if s=="file":
                        plot.filetype=1
                    if s=="text":
                        plot.texttype=1
                    if s=="pylab":
                        plot.pylabtype=1
                    if s=="gist":
                        plot.gisttype=1
                    if s=="feedback":
                        plot.feedbacktype=1
                self.dataTreeStore[path][column]=plot.gettype()
                #plot.type=new_text
            elif column==5:
                plot.preprocess=new_text
            elif column==6:
                plot.post=new_text
            elif column==7:
                try:
                    plot.dim=int(new_text)
                except:
                    self.dataTreeStore[path][column]=str(plot.dim)
            elif column==8:
                plot.xaxis=new_text
            elif column==9:
                plot.when=new_text
            elif column==10:
                try:
                    plot.info=eval(new_text)
                    if type(plot.info)!=types.DictType:
                        raise
                except:
                    self.dataTreeStore[path][column]=str(plot.info)
                
        else:
            print "Cannot edit while currently active"
##     def readGUIFeedback(self,handledList):
##         """Method to determine whether anything in the GUI xml section needs
##         changing - ie whether a button needs setting or not.
##         HandledList is a list of lists, which have:
##         [socket,"data"|"warning"|"error",tag,datadict].
##         This is not currently used. And hasn't been tested.
##         """
##         for plot in self.simdata.dataList:
##             if plot.feedbackFunc!=None:#expecting some sort of feedback
##                 for handledData in handledList:
##                     if handledData[1]=="data" and plot.tag==handledData[2] and handledData[3].has_key("feedback") and handledData[0] in plot.connList:
##                         #print "Setting feedback"
##                         path=str(self.simdata.dataList.index(plot))
##                         self.dataTreeStore[path][0]=int(handledData[3]["feedback"])
                        
            
    def plotToggled(self,cell,path):
        """button toggled by user"""
        #print "Toggled",cell,path,type(path)
        epath=eval("["+path.replace(":",",")+"]")
        plot=self.simdata.dataList#[int(path)]
        for i in epath:
            plot=plot.childList[i]
        plot.button=self.dataTreeStore[path]
        connList=self.getConns()
        plot.connList=connList
        if plot.cmd!=None and len(plot.cmd)>0:
            if plot.when in ["cmd","now"]:#send command but leave in same state
                tag=self.sim.execute(plot.cmd,rt=plot.ret,action=plot.when,connList=plot.connList)
                for conn in plot.connList:
                    self.sim.addCallback(conn,tag,plot.handle)
            elif plot.when[:3]=="rpt":
                self.dataTreeStore[path][0]=not self.dataTreeStore[path][0]
                curval=self.dataTreeStore[path][0]
                if curval:
                    plot.tag=self.sim.execute(plot.cmd,rt=plot.ret,action=plot.when,connList=plot.connList)
                    plot.repeating=1
                    for conn in plot.connList:
                        self.sim.addCallback(conn,plot.tag,plot.handle)
                else:
                    self.sim.execute(plot.cmd,tag=plot.tag,action="del",connList=plot.connList)
                    plot.repeating=0
                    plot.cancel()
                
    def execute(self,w,arg2=None):
        """Tell simulation to execute a command"""
        cmd=self.gladetree.get_widget("entryCommand").get_text()
        ret=self.gladetree.get_widget("entryReturn").get_text()
        freq=self.gladetree.get_widget("spinbuttonFreq").get_value_as_int()
        if ret=="":
            ret=" "
        if self.actionType[:3]=="rpt":
            self.actionType="rpt%d"%freq
        #print "Executing...[%s, %s, %s]"%(self.actionType,cmd,ret)
        self.commandHistoryList.append((cmd,ret,freq,self.actionType,self.plotUserData))
        tv=self.gladetree.get_widget("textviewCommand")
        tv.set_editable(1)
        buf=tv.get_buffer()
        buf.insert(buf.get_end_iter(),"%s -> %s\n"%(cmd,ret))
        tv.set_editable(0)
        tv.scroll_to_iter(buf.get_end_iter(),0.)
        self.commandHistoryListPos=len(self.commandHistoryList)-1
        connList=self.getConns()
        tag=self.sim.execute(cmd,rt=ret,action=self.actionType,connList=connList)
        if self.plotUserData:
            for conn in connList:
                plot=simdata.Plot()
                plot.ret=ret
                plot.cmd=cmd
                plot.title=cmd
                plot.pylabtype=1
                self.sim.addCallback(conn,tag,plot.handle)
        
    def listModules(self,w,arg2=None):
        """list simulation modules"""
        connList=self.getConns()
        for conn in connList:
            self.sim.dataProcessDict[(conn,"modulelist")]=self.printData
        self.sim.execute("data=ctrl.compListNames",tag="modulelist",rt="data",action="now",connList=connList)
    def listSim(self,w,arg2=None):
        """list button pressed by user"""
        txt=self.gladetree.get_widget("entryList").get_text()
        if txt=="":
            txt=None
        connList=self.getConns()
        for conn in connList:
            self.sim.dataProcessDict[(conn,"getObj")]=self.printList
        self.sim.getObjects(regexp=txt,connList=connList,readsock=0,indent="\t")
    def printData(self,data):
        """print out only the data"""
        if len(data)>3:
            if type(data[3])==types.DictType:
                for key in data[3].keys():
                    print data[3][key]
            else:
                print data[3]
    def printList(self,data):
        """print out the data"""
        #print "Got printList",data
        if len(data)>3 and type(data[3])==types.DictType and data[3].has_key("txt"):
            print data[3]["txt"]
        else:
            print "Unable to get list"
        
    def saveText(self,w,arg2=None):
        """not yet implemented"""
        pass
    def clearText(self,w,arg2=None):
        """Clear button has been pressed"""
        buf=self.gladetree.get_widget("textviewStdout").get_buffer()
        buf.delete(buf.get_start_iter(),buf.get_end_iter())
        
    def set_cmd(self,widget,arg2=None):
        self.actionType="cmd"
    def set_rpt(self,widget,arg2=None):
        self.actionType="rpt"
    def set_now(self,widget,arg2=None):
        self.actionType="now"
    def set_del(self,widget,arg2=None):
        self.actionType="del"

                   
    def send(self,lst):
        """Send stuff over socket"""
        if self.sock==None:
            print "Must connect first"
        else:
            serialise.Send(lst,self.sock)
    def updateSockList(self):
        """update the socket connected list"""
        self.sockListStore.clear()
        for conn in self.sim.connList:
            self.sockListStore.append(conn.getpeername())
    
    def getConns(self):
        """Get current sockets selected in selection box"""
        #read the selection from treeViewConnections.
        #Split this into hostname,port pairs, make into list.
        pathlist=[]
        self.sockView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)

        if len(pathlist)==0:
            connList=self.sim.connList
        else:
            selectedConns=[]
            connList=[]
            for path in pathlist:
                #print path
                connList.append(self.sim.connList[path[0]])
                #selectedConns.append(self.sim.connList[path[0]].getpeername())
##             for conn in self.sim.connList:
##                 if conn.getpeername() in selectedConns:
##                     connList.append(conn)
        return connList

class myStdOut:
    """Class for redirecting std output and error"""
    def __init__(self,textView,logfilename="simgui.log"):
        """Initialise"""
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
        """Update the text"""
        if self.logfile!=None:
            self.logfile.write(txt)
            self.logfile.flush()
        if self.textView!=None:
            self.buffer.insert(self.buffer.get_end_iter(),txt)
            self.textView.scroll_to_mark(self.buffer.get_insert(), 0.4)
            #self.textView.scroll_to_iter(self.buffer.get_end_iter(),0.5)

## class myFileSelection:
##     def __init__(self,title="Open an XML parameter file"):
##         self.fname=None
##         self.fs=gtk.FileSelection(title)
##         self.fs.connect("destroy", self.fsdestroy)
##         self.fs.ok_button.connect("clicked", self.file_ok_sel)
##         self.fs.cancel_button.connect("clicked",self.cancel)
##         self.fs.set_filename("*.xml")
##         self.fs.show()
##         gtk.main()
        
##     def file_ok_sel(self, w):
##         self.fname=self.fs.get_filename()
##         #print self.fname
##         self.fs.destroy()
##         gtk.main_quit()
##     def cancel(self,w):
##         self.fs.destroy()
##         gtk.main_quit()
##     def fsdestroy(self,w):
##         gtk.main_quit()

def run():
    """Run the gui"""
    g=simctrlGUI()
    gtk.main()

if __name__=="__main__":
    #import profile
    #profile.run("run()","profile.out")
    run()
