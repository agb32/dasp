#!/usr/bin/env python

import pygtk
pygtk.require("2.0")
import gtk, gobject
import gtk.glade as glade
import pango
import string,os,socket,os.path
import xml.parsers.expat
import util.parseSimXml as parseSimXml
import gui.dialog.dialog
from gui.myFileSelection.myFileSelection import myFileSelection
from sys import argv
class simsetup:
    def __init__(self):
        self.drawMode="full"#whether to do full or partial redraws.  full will be slower, but ensures all is correct...
        self.cwd=os.getcwd()
        try:
            gladefile=__file__.split("/")
            if gladefile[-2]=="bin":
                gladefile[-2]="simsetup"
            self.filepath=string.join(gladefile[:-1],"/")+"/"
            os.chdir(self.filepath)
        except:#older versions of python fail here...
            self.filepath="./"
            print "Guessing current directory for filepath (did you run this from the simsetup directory?)"
        self.drawlist=[]
        gladefile="simsetup.glade"
        self.gladetree=glade.XML(gladefile)
        self.sigdict={"on_togglebuttonNewModule_toggled":self.buttonNewModule,
                      "on_togglebuttonDelete_toggled":self.buttonDeleteModule,
                      "on_togglebuttonSelect_toggled":self.buttonSelect,
                      "on_togglebuttonConnect_toggled":self.buttonConnect,
                      "on_togglebuttonRePlace_toggled":self.buttonRePlace,
                      "on_togglebuttonShare_toggled":self.buttonShare,
                      "on_togglebuttonGroup_toggled":self.buttonGroup,
                      "on_buttonGenerate_clicked":self.buttonGenerate,
                      "on_window1_delete_event":self.deleteevent,
                      "on_window2_delete_event":self.deleteevent,
                      "on_window3_delete_event":self.closeParamGui,
                      "on_window4_delete_event":self.doWidgetGui,
                      "on_buttonWidOK_clicked":self.doWidgetGui,
                      "on_buttonWidCancel_clicked":self.doWidgetGui,
                      "on_add_widget1_activate":self.runWidgetGui,
                      "on_buttonParamGuiOk_clicked":self.closeParamGui,
                      "on_buttonParamGuiCancel_clicked":self.closeParamGui,
                      "on_window1_key_press_event":self.keypress,
                      "on_window1_key_release_event":self.keyrelease,
                      "on_window2_key_press_event":self.keypress,
                      "on_window2_key_release_event":self.keyrelease,
#                      "on_drawingarea1_configure_event":self.daconfig,
                      "on_drawingarea1_expose_event":self.daexpose,
                      "on_drawingarea1_realize":self.darealise,
                      "on_drawingarea1_button_press_event":self.dabuttonpressed,
                      "on_drawingarea1_key_press_event":self.dakey,
                      "on_drawingarea1_motion_notify_event":self.damotion,
                      "on_drawingarea1_enter_notify_event":self.daenter,
                      "on_drawingarea1_leave_notify_event":self.daleave,
                      "on_drawingarea1_button_release_event":self.dabuttonreleased,
                      "on_drawingarea1_drag_begin":self.dadragbegin,
                      "on_buttonFit_clicked":self.fitarea,
                      "on_togglebuttonSwallow_toggled":self.buttonSwallow,
                      "on_button2_clicked":self.button2,
                      "on_reload_widgets1_activate":self.loadSciWidgets,
                      "on_connect_to_paramgui_activate":self.connectParamGui,
                      "on_entryName_changed":self.nameChanged,
                      "on_checkbuttonGroupShare_toggled":self.groupShareToggled,
                      "on_entryCPU_changed":self.cpuChanged,
                      "on_entryArgs_changed":self.argsChanged,
                      "on_entryIDStr_changed":self.idstrChanged,
                      "on_checkbuttonFeedback_toggled":self.feedbackToggled,
                      "on_eventboxObjectPrefs_button_press_event":self.hideShowWidgets,
                      "on_eventboxInitTxt_button_press_event":self.hideShowWidgets,
                      "on_eventboxFinalTxt_button_press_event":self.hideShowWidgets,
                      "on_eventboxFinalisation_button_press_event":self.eventReleaseTxt,
                      "on_eventboxInitialisation_button_press_event":self.eventReleaseTxt,
                      "on_eventboxReleaseObjectPrefs_button_press_event":self.eventReleaseTxt,
                      "on_handleboxWidgets_child_detatched":self.childDetached,
                      "on_save1_activate":self.saveSetup,
                      "on_save_as1_activate":self.saveAsSetup,
                      "on_new1_activate":self.newSetup,
                      "on_open1_activate":self.openSetup,
                      "on_quit1_activate":self.deleteevent,
                      "on_print_drawlist1_activate":self.printDrawlist
                      
                     }
        self.gladetree.signal_autoconnect(self.sigdict)
        self.da=self.gladetree.get_widget("drawingarea1")
        #self.da.set_double_buffered(True)
        self.da.set_events(gtk.gdk.EXPOSURE_MASK |
                           gtk.gdk.LEAVE_NOTIFY_MASK |
                           gtk.gdk.BUTTON_PRESS_MASK |
                           gtk.gdk.POINTER_MOTION_MASK |
                           #gtk.gdk.POINTER_MOTION_HINT_MASK |
                           gtk.gdk.KEY_PRESS_MASK |
                           gtk.gdk.BUTTON_MOTION_MASK |
                           gtk.gdk.BUTTON_RELEASE_MASK |
                           gtk.gdk.KEY_RELEASE_MASK |
                           gtk.gdk.ENTER_NOTIFY_MASK |
                           gtk.gdk.LEAVE_NOTIFY_MASK |
                           gtk.gdk.VISIBILITY_NOTIFY_MASK)#these are set in glade - Events selection.
        self.gladetree.get_widget("textviewPrecode").get_buffer().connect("changed",self.textPrecodeChanged)
        self.gladetree.get_widget("textviewPostcode").get_buffer().connect("changed",self.textPostcodeChanged)
        self.gladetree.get_widget("textviewInitTxt").get_buffer().connect("changed",self.textInitChanged)
        self.gladetree.get_widget("textviewFinalTxt").get_buffer().connect("changed",self.textFinalChanged)

        #self.mincanvsizex=300
        #self.mincanvsizey=300
        self.canvsizex=300
        self.canvsizey=300
        self.currTagID=0
        self.cursors={}
        self.cursors["object"]=gtk.gdk.Cursor(gtk.gdk.GUMBY)
        self.cursors["select"]=None
        self.cursors["delete"]=gtk.gdk.Cursor(gtk.gdk.PIRATE)
        self.cursors["connect"]=gtk.gdk.Cursor(gtk.gdk.TARGET)
        self.cursors["replace"]=gtk.gdk.Cursor(gtk.gdk.TREK)
        self.cursors["share"]=gtk.gdk.Cursor(gtk.gdk.SB_H_DOUBLE_ARROW)
        self.cursors["group"]=gtk.gdk.Cursor(gtk.gdk.BOGOSITY)
        #self.da.window.set_cursor(self.cursors[self.doing]) to change cursor
        self.prefsVisible=1
        self.da.connect("configure-event",self.adjustCanvasSize)
        self.da.connect("expose-event", self.daexpose)
        self.gladetree.get_widget("window1").show()
        self.drawable=self.da.window
        self.style = self.da.get_style()
        self.gc = self.style.fg_gc[gtk.STATE_NORMAL]
        self.gc=self.style.black_gc
        #self.drawable.draw_line(self.gc,10,10,20,20)
        self.oldpyfilename="tmp.py"
        self.oldxmlfilename=""
        self.doing="object"
        self.paramguiSock=None
        self.grouplist=[]#list of group objects.
        self.ctrlPressed=0
        self.shiftPressed=0
        self.buttonpressed=0
        self.buttonpresspos=(0,0)
        self.selobj=None
        self.connectObject=None
        self.blitimg=None
        self.initTxt=""
        self.finalTxt=""
        self.vp=self.gladetree.get_widget("viewport1")
        self.vpx=self.vp.get_hadjustment()
        self.vpy=self.vp.get_vadjustment()
        self.tooltip=gtk.Tooltips()
        self.loadSciWidgets(fn="simwidgets.xml")
        self.modified(0)
        self.newObj=self.createNewSimObj()
        self.hideShowWidgets(self.gladetree.get_widget("eventboxInitTxt"))
        self.hideShowWidgets(self.gladetree.get_widget("eventboxFinalTxt"))
        self.gladetree.get_widget("checkbuttonGroupShare").hide()
        os.chdir(self.cwd)
        self.autoOpenFilename=None
        if len(argv)>1:
            self.autoOpenFilename=argv[1]
            if self.autoOpenFilename[-3:]==".py":
                self.autoOpenFilename=self.autoOpenFilename[:-3]+".xml"
            gobject.idle_add(self.openSetup,None)

    def prepareDrawable(self):
        self.da=self.gladetree.get_widget("drawingarea1")
        #self.da.set_double_buffered(True)
##         self.da.set_events(gtk.gdk.EXPOSURE_MASK |
##                            gtk.gdk.LEAVE_NOTIFY_MASK |
##                            gtk.gdk.BUTTON_PRESS_MASK |
##                            gtk.gdk.POINTER_MOTION_MASK |
##                            #gtk.gdk.POINTER_MOTION_HINT_MASK |
##                            gtk.gdk.KEY_PRESS_MASK |
##                            gtk.gdk.BUTTON_MOTION_MASK |
##                            gtk.gdk.BUTTON_RELEASE_MASK |
##                            gtk.gdk.KEY_RELEASE_MASK |
##                            gtk.gdk.ENTER_NOTIFY_MASK |
##                            gtk.gdk.LEAVE_NOTIFY_MASK |
##                            gtk.gdk.VISIBILITY_NOTIFY_MASK)#these are set in glade - Events selection.
        self.da.connect("configure-event",self.adjustCanvasSize)
        self.da.connect("expose-event", self.daexpose)
        self.da.connect("realize",self.darealise)
        self.drawable=self.da.window
        self.style = self.da.get_style()
        self.gc = self.style.fg_gc[gtk.STATE_NORMAL]
        self.gc=self.style.black_gc
        for obj in self.drawlist:
            obj.gc=self.gc
            obj.drawable=self.drawable
        self.newObj=self.createNewSimObj()
        

    def loadSciWidgets(self,w=None,fn=None):
        if fn==None:
            fn=myFileSelection("Choose the node/widget XML file","simwidgets.xml",complete="*.xml",parent=self.gladetree.get_widget("window2")).fname
        if fn==None:
            fn="simwidgets.xml"
        vbox=self.gladetree.get_widget("vboxSelections")
        vbox2=self.gladetree.get_widget("vboxSimSelect")
        try:
            txt=open(fn).read()
        except:
            txt=""
            print "file %s not found"%fn
        p = xml.parsers.expat.ParserCreate()
        self.simButtonList=[]
        self.cpuList=[]
        self.xmlInList=[]
        self.defaultObjAttribs={"name":"Unknown","shortname":"Unknown","comment":"Does nothing","pixmap":None,"textcol":"black","import":None,"object":None,"givesFeedback":0}
        self.defaultCpuAttribs={"pe":"1","processor":"1","colour":"white"}
        p.StartElementHandler = self.start_element
        p.EndElementHandler = self.end_element
        #p.CharacterDataHandler = self.char_data
        p.returns_unicode=0
        p.Parse(txt)
        if len(self.cpuList)==0:
            self.cpuList.append(self.defaultCpuAttribs.copy())
        if len(self.simButtonList)==0:
            self.simButtonList.append(self.defaultObjAttribs.copy())
        vbox.forall(vbox.remove)
        vbox2.forall(vbox2.remove)
        radButtonGroup=None
        indx=0
        for cpu in self.cpuList:
            #colours can be in form #rrggbb or name found in /usr/X11R6/lib/X11/rgb.txt

            #print "adding cpu",cpu["colour"]
            radButtonGroup=gtk.RadioButton(radButtonGroup,label="Node "+cpu["pe"]+", CPU "+cpu["processor"])
            radButtonGroup.connect("toggled", self.newCPUSelected,indx)
            eb = gtk.EventBox()
            eb.add(radButtonGroup)
            eb.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse(cpu["colour"]))
            eb.modify_bg(gtk.STATE_SELECTED, gtk.gdk.color_parse(cpu["colour"]))
            eb.modify_bg(gtk.STATE_ACTIVE, gtk.gdk.color_parse(cpu["colour"]))
            eb.modify_bg(gtk.STATE_PRELIGHT, gtk.gdk.color_parse(cpu["colour"]))
            vbox.pack_start(eb)
            indx+=1
        indx=0
        self.sciWidGroup=None
        for sim in self.simButtonList:
            self.addSciWidget(sim,indx)
##             hbox=gtk.HBox(False,0)
##             if sim["pixmap"]!=None:
##                 img=gtk.Image()
##                 img.set_from_file(sim["pixmap"])
##                 hbox.pack_start(img)
##             label=gtk.Label(sim["name"])
##             hbox.pack_start(label)
                
##             radButtonGroup=gtk.RadioButton(radButtonGroup)
##             if sim.has_key("comment"):
##                 self.tooltip.set_tip(radButtonGroup,sim["comment"])
##             radButtonGroup.add(hbox)
##             radButtonGroup.connect("toggled", self.newSimObjSelected,indx)
##             vbox2.pack_start(radButtonGroup)
            indx+=1
        self.selectedSimObj=0
        self.selectedCPU=0
        self.gladetree.get_widget("window2").show_all()

    def addSciWidget(self,sim,indx):
        """Add a science widget to the selection panel..."""
        vbox2=self.gladetree.get_widget("vboxSimSelect")
        hbox=gtk.HBox(False,0)
        if sim["pixmap"]!=None:
            img=gtk.Image()
            failed=0
##             print "trying",self.filepath,sim["pixmap"]
##             try:
##                 os.stat(sim["pixmap"])
##             except:
##                 failed=1
##             if failed==1:
##                 failed=0
##                 try:
##                     os.stat(self.filepath+sim["pixmap"])
##                 except:
##                     failed=1
##             if failed:
##                 print "File did not exist",self.filepath,sim["pixmap"]
            try:
                pb=gtk.gdk.pixbuf_new_from_file(sim["pixmap"])
                img.set_from_pixbuf(pb)
                hbox.pack_start(img)
            except:
                failed=1
            if failed:
                failed=0
                try:
                    pb=gtk.gdk.pixbuf_new_from_file(self.filepath+sim["pixmap"])
                    img.set_from_pixbuf(pb)
                    hbox.pack_start(img)
                except:
                    failed=1
                    
        label=gtk.Label(sim["name"])
        hbox.pack_start(label)

        self.sciWidGroup=gtk.RadioButton(self.sciWidGroup)
        if sim.has_key("comment"):
            self.tooltip.set_tip(self.sciWidGroup,sim["comment"])
        self.sciWidGroup.add(hbox)
        self.sciWidGroup.connect("toggled", self.newSimObjSelected,indx)
        vbox2.pack_start(self.sciWidGroup)
        vbox2.show_all()
        
    def runWidgetGui(self,w,d2=None):
        self.gladetree.get_widget("window4").show_all()
        gtk.main()
        self.gladetree.get_widget("window4").hide_all()

    def doWidgetGui(self,w,d2=None):
        if w==self.gladetree.get_widget("buttonWidOK"):
            wdict={"name":"entryWidName",
                   "comment":"entryWidComment",
                   "pixmap":"entryWidPixmap",
                   "shortname":"entryWidShort",
                   "textcol":"entryWidColour",
                   "import":"entryWidImport",
                   "object":"entryWidObject"
                   }
            d={}
            for key in wdict.keys():
                d[key]=self.gladetree.get_widget(wdict[key]).get_text()
            d["givesFeedback"]=str(int(self.gladetree.get_widget("checkbuttonWidFeedback").get_active()))
            self.simButtonList.append(d)
            self.addSciWidget(d,len(self.simButtonList)-1)
            print "ok",d
        gtk.main_quit()
    def printDrawlist(self,w,d2=None):
        print "drawlist",self.drawlist
                      
    def connectParamGui(self,w,d2=None):
        if w.get_active():
            self.gladetree.get_widget("window3").show_all()
            gtk.main()
            self.gladetree.get_widget("window3").hide_all()
            host=socket.gethostname()
            port=8990
            if self.paramGuiConnect:
                details=self.gladetree.get_widget("entryParamGuiPort").get_text().split(",")
                if len(details)==1:
                    port=int(details[0])
                elif len(details)==2:
                    port=int(details[1])
                    host=details[0]
                print "Connecting",host,port
                self.paramguiSock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                err=0
                try:
                    self.paramguiSock.connect((host,port))
                except:
                    print "Couldn't connect"
                    err=1
                if err:
                    w.set_active(0)
            else:
                w.set_active(0)
        else:
            if self.paramguiSock!=None:
                print "Disconnecting"
                self.paramguiSock.close()
                self.paramguiSock=None
                
    def closeParamGui(self,w=None,d2=None):
        self.paramGuiConnect=0
        if w==self.gladetree.get_widget("buttonParamGuiOk"):
            self.paramGuiConnect=1
        gtk.main_quit()
        return True
        
    def textPrecodeChanged(self,w,data1=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            if o.name=="sci":
                txt=w.get_text(w.get_start_iter(),w.get_end_iter()).strip()
                if txt!="#Python code to be placed\n#before the object is created\n#(use tab to indent)":
                    o.precode=txt
                    self.modified(1)
    def textPostcodeChanged(self,w,data1=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            if o.name=="sci":
                txt=w.get_text(w.get_start_iter(),w.get_end_iter()).strip()
                if txt!="#Python code to be placed\n#after the object is created\n#(use tab to indent)":
                    o.postcode=txt
                    self.modified(1)
    def textFinalChanged(self,w,data1=None):
        self.finalTxt=w.get_text(w.get_start_iter(),w.get_end_iter()).strip()
        self.modified(1)
    def textInitChanged(self,w,data1=None):
        self.initTxt=w.get_text(w.get_start_iter(),w.get_end_iter()).strip()
        self.modified(1)
    def modified(self,val):
        self.simmodified=val
        txt="Simulation setup: "
        if self.oldxmlfilename=="":
            txt+="Untitled"
        else:
            txt+=self.oldxmlfilename
        if val:
            txt+=" *"
        self.gladetree.get_widget("window1").set_title(txt)
    def saveSetup(self,w=None,d=None):
        if self.oldxmlfilename=="":
            self.saveAsSetup()
        if self.oldxmlfilename!="":#now have valid name...
            fn=self.oldxmlfilename
            try:
                oldtxt=open(fn).read()
            except:
                oldtxt=""
            if len(oldtxt)>0:#create backup...
                open(fn+"~","w").write(oldtxt)
                print "todo: warning message - overwrite?"
            spos=oldtxt.find("<aosim")
            epos=oldtxt.find("</aosim")
            sspos=oldtxt.find("<simSetup")
            sepos=oldtxt.find("</simSetup")
            if spos==-1 or epos==-1 or sspos==-1 or sepos==-1:
                if len(oldtxt)>0:
                    print "Warning - overwriting current file %s"%fn
                pretxt="<aosim>\n<simSetup>\n"
                posttxt="</simSetup>\n</aosim>\n"
            else:
                l=len(oldtxt)
                while sspos<l and oldtxt[sspos]!=">":
                    sspos+=1
                pretxt=oldtxt[:sspos+1]+"\n"#ends with "<simSetup>"
                posttxt=oldtxt[sepos:]
            midtxt=self.createSetupXml(all=0)
            newtxt=pretxt+midtxt+posttxt
            open(fn,"w").write(newtxt)
            self.modified(0)
    def saveAsSetup(self,w=None,d=None):
        fn=myFileSelection("Choose an output XML file",self.oldxmlfilename,complete="*.xml",parent=self.gladetree.get_widget("window2")).fname
        if fn!=None:
            self.oldxmlfilename=fn
            if w!=None:
                self.saveSetup()
                
    def openSetup(self,w=None,d=None):
        if self.autoOpenFilename==None:
            fn=myFileSelection("Choose a simulation setup XML file",self.oldxmlfilename,complete="*.xml",parent=self.gladetree.get_widget("window2")).fname
        else:
            fn=self.autoOpenFilename
        self.autoOpenFilename=None
        if fn!=None:
            self.drawlist=[]
            self.selobj=None
            self.oldxmlfilename=fn
            txt=open(fn).read()
            p=parseSimXml.parseSimXml(txt)
            self.initTxt=p.globalPreCode
            self.finalTxt=p.globalPostCode
            self.gladetree.get_widget("textviewInitTxt").get_buffer().set_text(self.initTxt)
            self.gladetree.get_widget("textviewFinalTxt").get_buffer().set_text(self.finalTxt)

            maxtag=0
            for obj in p.objList:
                if obj.type=="simObj":
                    if obj.tag>maxtag:
                        maxtag=obj.tag
                    colour="white"
                    for cpu in self.cpuList:
                        #print cpu["pe"],cpu["processor"],obj.cpu,type(obj.cpu),type(cpu["pe"])
                        if int(cpu["pe"])==obj.cpu[0] and int(cpu["processor"])==obj.cpu[1]:
                            colour=cpu["colour"]
                            break
                    #dcpu=self.cpuList[obj.cpu]
                    os.chdir(self.filepath)
                    nobj=simObj(self.drawable,self.gc,"sci",obj.tag,x=obj.pos[0],y=obj.pos[1],msg=obj.shortname,textcol=obj.textcol,pango_context=self.da.create_pango_context(),pixmap=obj.pixmap,colour=colour,rank=obj.cpu,imp=obj.imp,object=obj.object,feedback=obj.feedback)
                    os.chdir(self.cwd)
                    nobj.lines=obj.lines
                    nobj.endlines=obj.endlines
                    nobj.connectto=obj.connectto
                    nobj.connectfrom=obj.connectfrom
                    nobj.sharedTo=obj.sharedTo
                    if obj.dictionary.has_key("idstr"):
                        nobj.idstr=obj.dictionary["idstr"]
                    nobj.precode=obj.precode.replace("    ","\t")
                    nobj.precode=nobj.precode.replace("\n\t","\n").strip()#remove indent
                    nobj.postcode=obj.postcode.replace("    ","\t")
                    nobj.postcode=nobj.postcode.replace("\n\t","\n").strip()
                    nobj.pyname=obj.pyname
                    nobj.groupShare=obj.groupShare
                    nobj.args=obj.args
                    nobj.parentNames=obj.parentNames
                    nobj.setText()
                    self.drawlist.append(nobj)
                elif obj.type=="groupShare":
                    nobj=boxObj(self.drawable,self.drawable.new_gc(self.drawable.get_colormap().alloc_color("blue")),obj.coordList)
                    nobj.cpuList=obj.cpuList
                    nobj.idstr=obj.idstr
                    nobj.calcArea()
                    self.drawlist.append(nobj)
                    self.grouplist.append(nobj)
            #now go through and add the line objects, and update the lines, endlines, tags etc...
            linelist=[]
            for obj in self.drawlist:
                if obj.name=="sci":
                    for i in range(len(obj.connectto)):
                        tag=obj.connectto[i]
                        cobj=self.getObjWithTag(tag,self.drawlist)
                        obj.lines[i].append(cobj.getHandlePos("top"))
                        obj.lines[i].insert(0,obj.getHandlePos("bottom"))
                        #cobj.endlines[cobj.connectfrom.index(obj.tagID)]=obj.lines
                        l=lineObj(self.drawable,self.gc,obj.lines[i])#(int(e.x),int(e.y))])
                        obj.lines[i]=l
                        cobj.endlines[cobj.connectfrom.index(obj.tagID)]=l
                        l.idstr=cobj.parentNames[cobj.connectfrom.index(obj.tagID)]
                        l.calcArea()
                        linelist.append(l)
                    sharedTo=obj.sharedTo
                    obj.sharedTo=[]
                    for i in range(len(sharedTo)):
                        tag=sharedTo[i]
                        cobj=self.getObjWithTag(tag,self.drawlist)
                        if obj.rank==cobj.rank:
                            l=lineObj(self.drawable,self.drawable.new_gc(self.drawable.get_colormap().alloc_color("green")),[obj.getHandlePos("main"),cobj.getHandlePos("main")])#(int(e.x),int(e.y))])
                            l.calcArea()
                            linelist.append(l)
                            obj.sharedTo.append((tag,l))
                            cobj.sharedFrom.append((obj.tagID,l))
                    
            self.drawlist+=linelist
            self.currTagID=maxtag+1
            self.newObj=self.createNewSimObj()
            self.fitarea()
            self.modified(0)
            self.daexpose()
    def newSetup(self,w=None,d=None):
        if self.simmodified:
            msg=" (current simulation is modified)"
        else:
            msg=""
        if gui.dialog.dialog.myDialog(msg="New project?\nAre you sure%s?\n"%msg,title="New projet?",parent=self.gladetree.get_widget("window1")).resp=="ok":
            self.drawlist=[]
            self.selobj=None
            self.currTagID=0
            self.newObj=self.createNewSimObj()
            self.oldxmlfilename=""
            self.modified(0)
            self.daexpose()
    def getObjWithTag(self,tag,objlist):
        rt=None
        for obj in objlist:
            try:
                t=obj.tagID
            except:
                t=None
                #print "Error in getObjWith Tag, obj has no tagID attribute (probably a line!)"
                #print obj
                #raise
            
            if t==tag:
                rt=obj
                break
        if rt==None:
            print "Warning - no object with tag %s found"%str(tag)
        return rt
    def keypress(self,w,e=None):
        #print "todo",e,e.type,e.keyval,e.string,hex(e.state)

        if e.state & gtk.gdk.CONTROL_MASK:
            #pressed with control
            if e.keyval>0 and e.keyval<256:
                ch=chr(e.keyval)
            else:
                ch=None
            f={"s":self.saveSetup,
               "n":self.newSetup,
               "o":self.openSetup,
               }
            d={"e":"togglebuttonSelect",
               "c":"togglebuttonConnect",
               "m":"togglebuttonNewModule",
               "r":"togglebuttonRePlace",
               "d":"togglebuttonDelete",
               "h":"togglebuttonShare",
               "g":"togglebuttonGroup"
               }

            if f.has_key(ch):
                f[ch]()
            elif ch=="S":
                self.saveAsSetup(1)
            elif d.has_key(ch):
                self.gladetree.get_widget(d[ch]).set_active(1)
                self.da.window.set_cursor(self.cursors[self.doing])#change cursor
            elif ch=="p":
                self.gladetree.get_widget("buttonGenerate").clicked()


                
##             if ch=='s':#115:#s
##                 self.saveSetup()
##             elif ch=='n':#110:#n
##                 self.newSetup()
##             elif ch=='o':#111:#o
##                 self.openSetup()
##             elif ch=='d':#100:#d
##                 pass
##             elif ch=='c':#99:#c
##                 pass
##             elif ch=='g':#103:#g
##                 pass
##             elif ch=='p':#112:#p
##                 pass
##             elif ch=='r':#114:#r
##                 pass
##             elif ch=='e':#101:#e
##                 pass
##             elif ch=='m':#109:#m
##                 pass
##             elif ch=='S':
##                 print "saveas"
            
##             else:
##                 print "keyval=",e.keyval
##             if e.state&gtk.gdk.SHIFT_MASK:
##                 if e.keyval==83:
##                     self.saveAsSetup(1)
##         else:
##             if not (e.state&gtk.gdk.SHIFT_MASK):
##                 d={"s":"togglebuttonSelect",
##                    "c":"togglebuttonConnect",
##                    "n":"togglebuttonNewModule",
##                    "r":"togglebuttonRePlace",
##                    "d":"togglebuttonDelete",
##                    "m":"togglebuttonShare",
##                    "g":"togglebuttonGroup"
##                    }
                
##                 if d.has_key(e.string):
##                     self.gladetree.get_widget(d[e.string]).set_active(1)
##                     self.da.window.set_cursor(self.cursors[self.doing])#change cursor
                    
##                 elif e.string=="p":
##                     self.gladetree.get_widget("buttonGenerate").clicked()
                    
        #if e.keyval==65507 or e.keyval==65508:
        #    self.ctrlPressed=1
        #if e.keyval==65505 or e.keyval==65506:
        #    self.shiftPressed=1
    def keyrelease(self,w,e=None):
        #if e.keyval==65507 or e.keyval==65508:
        #    self.ctrlPressed=0
        #if e.keyval==65505 or e.keyval==65506:
        #    self.shiftPressed=0
        pass

    def nameChanged(self,w,data=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            o.pyname=w.get_text().strip()
            self.modified(1)
    def groupShareToggled(self,w,data=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            o.groupShare=int(w.get_active())
            self.modified(1)
        
    def cpuChanged(self,w,data=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            if o.name=="sci":
                self.modified(1)
                txt=w.get_text().split(",")
                if len(txt)>0:
                    try:
                        node=int(txt[0])
                    except:
                        node=None
                cpu=1
                if len(txt)>1:
                    try:
                        cpu=int(txt[1])
                    except:
                        pass
                if node!=None:
                    for currcpu in self.cpuList:
                        if int(currcpu["pe"])==node and int(currcpu["processor"])==cpu:
                            #change the colour...
                            #Note, even if not in current list, should still
                            #allow the user to change it...
                            o.setColour(currcpu["colour"])
                            break
                    o.rank=(node,cpu)
                    o.setText()
                    o.draw()
                    o.drawSelected()
            elif o.name=="box":
                self.modified(1)
                txt=w.get_text()
                o.cpuList=txt
    def idstrChanged(self,w,data=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            self.modified(1)
            o.idstr=w.get_text().strip()
            o.setText()
            o.draw()
            o.drawSelected()
    def argsChanged(self,w,data=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            if o.name=="sci":
                self.modified(1)
                o.args=w.get_text().strip()
                o.setText()
                o.draw()
                o.drawSelected()
    def feedbackToggled(self,w,data=None):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            if o.name=="sci":
                self.modified(1)
                o.feedback=w.get_active()
    def newSimObjSelected(self,w,data1=None,data2=None):
        if w.get_active():
            self.selectedSimObj=data1
            self.newObj=self.createNewSimObj()
            self.gladetree.get_widget("togglebuttonNewModule").set_active(1)
            
    def newCPUSelected(self,w,data1=None):
        if w.get_active():
            self.selectedCPU=data1
            self.newObj=self.createNewSimObj()

    def hideShowWidgets(self,w,d=None):
        fw=w.get_parent().get_children()[0]
        if fw.flags()&gtk.VISIBLE:
            fw.hide()
        else:
            fw.show_all()
        return True
    def start_element(self,name,attrs):
        self.xmlInList.append(name)
        if "aosim" in self.xmlInList:
            if "simsetup" in self.xmlInList:
                if name=="simobj":
                    for key in self.defaultObjAttribs.keys():
                        if not attrs.has_key(key):
                            attrs[key]=self.defaultObjAttribs[key]
                    self.simButtonList.append(attrs)
                elif name=="cpu":
                    for key in self.defaultCpuAttribs.keys():
                        if not attrs.has_key(key):
                            attrs[key]=self.defaultCpuAttribs[key]
                    self.cpuList.append(attrs)
    def end_element(self,name):
        if self.xmlInList[-1]!=name:
            print "XML error..."
            raise Exception("XML parsing error")
        else:
            self.xmlInList.pop()
                
    def button2(self,w=None,a=None):
        print "**********************************button2"
        print self.drawable.get_size()
        print self.da.size_request()
        print self.vp.get_pointer()
        print self.vp.size_request()
        print self.vp.get_size_request()
        print self.vpx.value,self.vpx.lower,self.vpx.upper,self.vpx.page_size
        print self.vpx.get_value()
        print self.vpy.get_value()
        self.da.set_size_request(600,600)
        print self.da.window
        self.testbutton=gtk.Button(label="test")
        self.testbutton.set_parent_window(self.da.window)
        self.testbutton.show()
        #self.da.queue_resize()
    def adjustCanvasSize(self,w=None,a=None):
        #print "adjustCanvasSize",w,a,self.canvsizex,self.canvsizey
        if self.canvsizex<300:
            self.canvsizex=300
        if self.canvsizey<300:
            self.canvsizey=300
        if self.canvsizex>3000:
            self.canvsizex=3000
        if self.canvsizey>3000:
            self.canvsizey=3000
        if self.da.get_size_request()!=(int(self.canvsizex),int(self.canvsizey)):
            self.da.set_size_request(int(self.canvsizex),int(self.canvsizey))
            # need to see what current size is and adjust canvsizex/y to this.
            # ie since when initiated will be >300,300.
            self.gladetree.get_widget("labelInfo").set_text("%d x %d"%(self.canvsizex,self.canvsizey))
##     def drawobject(self,x,y):
##         x=int(x)
##         y=int(y)
##         self.drawable.draw_rectangle(self.gc,True,x-10,y-10,20,20)
##         self.drawlist.append(simObj(self.drawable,self.gc,"rect",x-10,y-10,20,20))
##         return (20,20)#return object width, height.
    
    def getobject(self,x,y,remove=1,name=None):
        """Gets the object currently under position x,y."""
        obj=None
        for i in range(len(self.drawlist)-1,-1,-1):
            d=self.drawlist[i]
            if d!=None:
                if name==None or d.name==name:
                    if d.atPosition(x,y):
                        #if d.x<x and d.x+d.w>x and d.y<y and d.y+d.h>y:
                        if remove:
                            obj=self.drawlist.pop(i)
                        else:
                            obj=self.drawlist[i]
                        break
        return obj
    def shiftAll(self,x,y):
        """Shift everythings position by x,y"""
        #print "shifting by",x,y
        for d in self.drawlist:
            if d!=None:
                d.shift(x,y)
                #d.x+=x
                #d.y+=y
    def fitarea(self,w=None,a=None):
        """Fit the canvas size to widgets on it."""
        minx=None
        miny=None
        maxx=0
        maxy=0
        for d in self.drawlist:
            if d!=None:
                #print d.x,d.y,d.x+d.w,d.y+d.h
                if minx==None or d.x<minx:
                    minx=d.x
                if miny==None or d.y<miny:
                    miny=d.y
                if d.x+d.w>maxx:
                    maxx=d.x+d.w
                if d.y+d.h>maxy:
                    maxy=d.y+d.h
        
        minx-=30
        miny-=30
        maxx+=30
        maxy+=30
        #print "max, min",minx,maxx,miny,maxy
        #if minx<0: minx=0
        #if miny<0: miny=0
        #if minx>0 or miny>0:
        #    self.shiftAll(-minx,-miny)
        self.shiftAll(-minx,-miny)
        maxx-=minx
        maxy-=miny
        self.canvsizex=maxx
        self.canvsizey=maxy
        self.adjustCanvasSize()
        self.modified(1)
        self.daexpose()
        
    def buttonNewModule(self,w,a=None):
        w=self.gladetree.get_widget("togglebuttonNewModule")
        if w.get_active():
            self.doing="object"
            self.newObj=self.createNewSimObj()
            #self.newObj=simObj(self.drawable,self.gc,"sci",colour="red",pango_context=self.da.create_pango_context())
            self.gladetree.get_widget("togglebuttonDelete").set_active(0)
            self.gladetree.get_widget("togglebuttonConnect").set_active(0)
            self.gladetree.get_widget("togglebuttonSelect").set_active(0)
            #self.gladetree.get_widget("togglebuttonNewModule").set_active(0)
            self.gladetree.get_widget("togglebuttonRePlace").set_active(0)
            self.gladetree.get_widget("togglebuttonShare").set_active(0)
            self.gladetree.get_widget("togglebuttonGroup").set_active(0)
    def buttonDeleteModule(self,w,a=None):
        w=self.gladetree.get_widget("togglebuttonDelete")
        if w.get_active():
            self.doing="delete"
            #self.gladetree.get_widget("togglebuttonDelete").set_active(0)
            self.gladetree.get_widget("togglebuttonConnect").set_active(0)
            self.gladetree.get_widget("togglebuttonSelect").set_active(0)
            self.gladetree.get_widget("togglebuttonNewModule").set_active(0)
            self.gladetree.get_widget("togglebuttonRePlace").set_active(0)
            self.gladetree.get_widget("togglebuttonShare").set_active(0)
            self.gladetree.get_widget("togglebuttonGroup").set_active(0)

    def buttonSelect(self,w,a=None):
        w=self.gladetree.get_widget("togglebuttonSelect")
        if w.get_active():
            self.doing="select"
            self.gladetree.get_widget("togglebuttonDelete").set_active(0)
            self.gladetree.get_widget("togglebuttonConnect").set_active(0)
            #self.gladetree.get_widget("togglebuttonSelect").set_active(0)
            self.gladetree.get_widget("togglebuttonNewModule").set_active(0)
            self.gladetree.get_widget("togglebuttonRePlace").set_active(0)
            self.gladetree.get_widget("togglebuttonShare").set_active(0)
            self.gladetree.get_widget("togglebuttonGroup").set_active(0)
            
    def buttonConnect(self,w,a=None):
        w=self.gladetree.get_widget("togglebuttonConnect")
        if w.get_active():
            self.doing="connect"
            self.connectObject=None
            self.gladetree.get_widget("togglebuttonDelete").set_active(0)
            #self.gladetree.get_widget("togglebuttonConnect").set_active(0)
            self.gladetree.get_widget("togglebuttonSelect").set_active(0)
            self.gladetree.get_widget("togglebuttonNewModule").set_active(0)
            self.gladetree.get_widget("togglebuttonRePlace").set_active(0)
            self.gladetree.get_widget("togglebuttonShare").set_active(0)
            self.gladetree.get_widget("togglebuttonGroup").set_active(0)

    def buttonRePlace(self,w,a=None):
        w=self.gladetree.get_widget("togglebuttonRePlace")
        if w.get_active():
            self.doing="replace"
            self.gladetree.get_widget("togglebuttonDelete").set_active(0)
            self.gladetree.get_widget("togglebuttonConnect").set_active(0)
            self.gladetree.get_widget("togglebuttonSelect").set_active(0)
            self.gladetree.get_widget("togglebuttonNewModule").set_active(0)
            #self.gladetree.get_widget("togglebuttonRePlace").set_active(0)
            self.gladetree.get_widget("togglebuttonShare").set_active(0)
            self.gladetree.get_widget("togglebuttonGroup").set_active(0)

    def buttonShare(self,w,a=None):
        w=self.gladetree.get_widget("togglebuttonShare")
        if w.get_active():
            self.doing="share"
            self.connectObject=None
            self.gladetree.get_widget("togglebuttonDelete").set_active(0)
            self.gladetree.get_widget("togglebuttonConnect").set_active(0)
            self.gladetree.get_widget("togglebuttonSelect").set_active(0)
            self.gladetree.get_widget("togglebuttonNewModule").set_active(0)
            self.gladetree.get_widget("togglebuttonRePlace").set_active(0)
            #self.gladetree.get_widget("togglebuttonShare").set_active(0)
            self.gladetree.get_widget("togglebuttonGroup").set_active(0)

    def buttonGroup(self,w,a=None):
        w=self.gladetree.get_widget("togglebuttonGroup")
        if w.get_active():
            self.doing="group"
            self.connectObject=None
            self.gladetree.get_widget("togglebuttonDelete").set_active(0)
            self.gladetree.get_widget("togglebuttonConnect").set_active(0)
            self.gladetree.get_widget("togglebuttonSelect").set_active(0)
            self.gladetree.get_widget("togglebuttonNewModule").set_active(0)
            self.gladetree.get_widget("togglebuttonRePlace").set_active(0)
            self.gladetree.get_widget("togglebuttonShare").set_active(0)
            #self.gladetree.get_widget("togglebuttonGroup").set_active(0)


    def buttonGenerate(self,w,a=None):
        fn=myFileSelection("Choose an output python file",self.oldxmlfilename[:-3]+"py",complete="",parent=self.gladetree.get_widget("window2")).fname
        if fn!=None:
            self.oldfilename=fn
            txt=self.createSetupXml()
            print txt
            p=parseSimXml.parseSimXml(txt)
            pycode=p.makePython()
            if os.path.exists(fn):
                #save any user added parts...
                lines=open(fn).readlines()
                stxt=""
                etxt=""
                adds=0
                adde=0
                for line in lines:
                    if adds==1 and line[:14]=="if ctrl.rank==":
                        adds=0
                    if adds==1:
                        stxt+=line
                    if line=="#Add any personal code after this line and before the next, and it won't get overwritten\n":
                        adds=1
                    if adde==1:
                        etxt+=line
                    if line=="#Add any personal code after this, and it will not get overwritten\n":
                        adde=1
                #now add stxt and etxt in the right places...
                etxt=etxt.strip()
                if len(etxt)==0:
                    etxt="Scientific.MPI.world.abort(0)\n"
                if len(stxt)>0 or len(etxt)>0:
                    lines=pycode.split("\n")
                    f=open(fn,"w")
                    for line in lines:
                        f.write(line+"\n")
                        if line=="#Add any personal code after this line and before the next, and it won't get overwritten":
                            f.write(stxt)
                        if line=="#Add any personal code after this, and it will not get overwritten":
                            f.write(etxt)
                            if len(etxt)>0 and etxt[-1]!="\n":
                                f.write("\n")
                    f.close()
                else:
                    open(fn,"w").write(pycode)
                    
            else:
                open(fn,"w").write(pycode)
    def buttonSwallow(self,w,a=None):
        if w.get_active():
            print "Swallowing"
            vbox=self.gladetree.get_widget("vbox1")
            win=self.gladetree.get_widget("window1")
            hpane=self.gladetree.get_widget("hpanedToolbar")
            #self.da.unrealize()
            win.remove(vbox)
            hpane.add2(vbox)
            win.hide()
            self.prepareDrawable()
            hpane.show_all()

        else:
            print "Retching"
            vbox=self.gladetree.get_widget("vbox1")
            hpane=self.gladetree.get_widget("hpanedToolbar")
            win=self.gladetree.get_widget("window1")
            hpane.remove(vbox)
            win.add(vbox)
            win.show_all()
            self.prepareDrawable()
    def createSetupXml(self,all=1,prepost=1):
        remlist=[]
        txt=""
        #first check that all sharings are legal...
        for d in self.drawlist:
            if d.name=="sci":
                remlist=[]
                for t,l in d.sharedTo:
                    o=self.getObjWithTag(t,self.drawlist)
                    if o==None or o.rank!=d.rank:
                        remlist.append((t,l))
                if len(remlist)>0:
                    print "Removing resource sharing for object not on same node (or not existing): tags: %s"%str(remlist)
                for r in remlist:
                    d.sharedTo.remove(r)
        if all:
            txt+="<aosim>\n<simSetup>\n"
        if prepost:
            if self.initTxt!="":
                txt+="<precode>\n%s\n</precode>\n"%self.initTxt
            if self.finalTxt!="":
                txt+="<postcode>\n%s\n</postcode>\n"%self.finalTxt
        remlist=[]
        for d in self.drawlist:
            if d.removed:
                remlist.append(d)
            else:
                if d.name=="sci":
                    txt+=d.__repr__()
                elif d.name=="box":
                    txt+=d.__repr__()
        for d in remlist:
            self.drawlist.remove(d)
        if all:
            txt+="</simSetup>\n</aosim>\n"
        return txt
    def deleteevent(self,w,a=None):
        print "deleteevent"
        if self.simmodified:
            msg=" (current simulation is modified)"
        else:
            msg=""
        if gui.dialog.dialog.myDialog(msg="Quit?\nAre you sure%s?\n"%msg,title="Really quit?!?",parent=self.gladetree.get_widget("window1")).resp=="ok":
            gtk.main_quit()
            return False
        return True
    def daconfig(self,w,a=None):
        print "daconfig"
    def daexpose(self,wid=None,a=None,rect=None):
        """Rect can be a tuple of the rectangle to clear, (x0,y0,x1,y1)
        """
        x=y=0
        w,h=self.drawable.get_size()
        w=int(w)
        h=int(h)
        if self.drawMode!="full" and rect!=None:
            try:
                x,y,w,h=rect
                x1=x+w
                y1=y+h
            except:
                print "error with expose rectangle",rect
                x=y=0
                pass
        self.drawable.draw_rectangle(self.style.white_gc,True,x,y,w,h)#clear it
        #self.drawable.draw_line(self.gc,10,10,20,20)
        remlist=[]
        for d in self.drawlist:
            if d.removed:
                remlist.append(d)
            else:
                d.draw(rect=(x,y,w,h))#redraw the object if any part lies in the area.
            #d.draw()
        for d in remlist:
            self.drawlist.remove(d)
        self.updateSelectedObject()#draw mark at the object
    def fillPrefs(self):
        if self.selobj!=None and self.selobj.has_key("obj"):
            o=self.selobj["obj"]
            self.gladetree.get_widget("entryIDStr").set_text(o.idstr)
            if o.name=="sci":
                self.gladetree.get_widget("entryName").set_sensitive(1)
                self.gladetree.get_widget("checkbuttonGroupShare").set_sensitive(1)
                if self.isInGroup(o)!=None:
                    self.gladetree.get_widget("labelPrefName").set_text("Share")
                    self.gladetree.get_widget("entryName").hide()
                    self.gladetree.get_widget("checkbuttonGroupShare").set_active(o.groupShare)
                    self.gladetree.get_widget("checkbuttonGroupShare").show()

                else:
                    self.gladetree.get_widget("labelPrefName").set_text("Name")
                    self.gladetree.get_widget("entryName").set_text(o.pyname)
                    self.gladetree.get_widget("checkbuttonGroupShare").hide()
                    self.gladetree.get_widget("entryName").show()
                w=self.gladetree.get_widget("entryCPU")
                w.set_sensitive(1)
                w.set_text("%s, %s"%(o.rank[0],o.rank[1]))
                w=self.gladetree.get_widget("entryArgs")
                w.set_sensitive(1)
                w.set_text(o.args)
                w=self.gladetree.get_widget("checkbuttonFeedback")
                w.set_sensitive(1)
                w.set_active(int(o.feedback))
                w=self.gladetree.get_widget("textviewPrecode")
                w.set_sensitive(1)
                w.get_buffer().set_text(o.precode)
                w=self.gladetree.get_widget("textviewPostcode")
                w.set_sensitive(1)
                w.get_buffer().set_text(o.postcode)
            elif o.name=="box":
                #self.gladetree.get_widget("labelPrefName").set_text("Share")
                self.gladetree.get_widget("entryName").set_sensitive(0)
                w=self.gladetree.get_widget("entryCPU")
                w.set_sensitive(1)
                w.set_text("%s"%(o.cpuList))
                w=self.gladetree.get_widget("entryArgs")
                w.set_sensitive(0)
                w=self.gladetree.get_widget("checkbuttonFeedback")
                w.set_sensitive(0)
                w=self.gladetree.get_widget("textviewPrecode")
                w.set_sensitive(0)
                w=self.gladetree.get_widget("textviewPostcode")
                w.set_sensitive(0)
                self.gladetree.get_widget("checkbuttonGroupShare").set_sensitive(0)
                
            else:
                self.gladetree.get_widget("labelPrefName").set_text("Name")
                self.gladetree.get_widget("entryName").set_sensitive(1)
                self.gladetree.get_widget("entryName").set_text(o.pyname)
                w=self.gladetree.get_widget("entryCPU")
                w.set_sensitive(0)
                w=self.gladetree.get_widget("entryArgs")
                w.set_sensitive(0)
                w=self.gladetree.get_widget("checkbuttonFeedback")
                w.set_sensitive(0)
                w=self.gladetree.get_widget("textviewPrecode")
                w.set_sensitive(0)
                w=self.gladetree.get_widget("textviewPostcode")
                w.set_sensitive(0)
                self.gladetree.get_widget("checkbuttonGroupShare").set_sensitive(0)


    def updateSelectedObject(self):
        if self.selobj!=None and self.selobj.has_key("obj"):
            self.selobj["obj"].drawSelected()
            self.fillPrefs()

    def darealise(self,w,a=None):
        #print "darealise"
        pass

    def doParamGui(self,obj):
        if self.paramguiSock!=None and obj.name=="sci":
            try:
                self.paramguiSock.send("self.findOrAdd('%s','%s','%s')"%(obj.imp.split(".")[-1],obj.imp,obj.getIDStr()))
            except:
                self.gladetree.get_widget("connect_to_paramgui").set_active(0)
        
    def dabuttonpressed(self,w,e=None):
        """If select, getobject(), and prepare to place it at end of list
        once button is released, possibly in new position or something"""
        #print "dabuttonpressed",type(w),type(e)
        self.buttonpressed=1
        
        self.buttonpresspos=(e.x,e.y)
        if self.doing=="object":#about to place object
            #print "objpressed"
            pass
        elif self.doing=="replace":#re-rank an object...
            if self.drawMode=="full":
                tmp=self.getobject(e.x,e.y,remove=0)
            else:
                tmp=self.getobject(e.x,e.y)#removes the object from drawlist.
            if tmp!=None:
                if tmp.name=="sci":
                    self.modified(1)
                    dcpu=self.cpuList[self.selectedCPU]
                    rank=(int(dcpu["pe"]),int(dcpu["processor"]))
                    col=dcpu["colour"]
                    tmp.rank=rank
                    tmp.setText()
                    tmp.setColour(col)
                self.doParamGui(tmp)
                self.selobj={"obj":tmp}
                self.daexpose()
        elif self.doing=="select":#select an object...
            #find out which object (if any) has been selected.
            if self.drawMode=="full":
                tmp=self.getobject(e.x,e.y,remove=0)
            else:
                tmp=self.getobject(e.x,e.y)#removes the object from the drawlist.
            if tmp!=None:
                self.selobj={"obj":tmp}
                rect=self.selobj["obj"].getRect()
                #rect=None
                self.daexpose(None,rect=rect)#redraw area without the object.
                #rect=self.selobj["obj"].getRect()
                self.selobj["blit"]=self.drawable.get_image(*rect)#copy the image.
                if self.drawMode!="full":# now redraw the object.
                    self.selobj["obj"].draw()
                    self.updateSelectedObject()#draw mark at the object
                # self.drawable.draw_rectangle(self.gc,True,self.selobj["obj"].x,self.selobj["obj"].y,self.selobj["obj"].w,self.selobj["obj"].h)
                self.selobj["coords"]=(e.x,e.y)
                self.selobj["offset"]=(e.x-self.selobj["obj"].x,e.y-self.selobj["obj"].y)
                self.doParamGui(tmp)
            else:
                self.selobj=None
        elif self.doing=="delete":
            pass
        elif self.doing=="connect":
            if self.connectObject==None:#start new connection...
                tmp=self.getobject(e.x,e.y,remove=0,name="sci")
                if tmp!=None:
                    if tmp.name=="sci":
                        if tmp.selected=="main":
                            tmp.selected="bottom"
                            tmp.centreCoords=tmp.getHandlePos("bottom")
                        #have selected top or bottom of object.
                        self.connectObject=tmp
                        if self.connectObject.selected=="top":
                            #draw from top to mouse...
                            l=lineObj(self.drawable,self.gc,[self.connectObject.centreCoords,self.connectObject.getHandlePos(self.connectObject.selected)])#(int(e.x),int(e.y))])
                            l.drawTmp()
                            self.connectObject.currentLine=l#a line ends at the top, ie dataflow ends at the input of an object.
                        elif self.connectObject.selected=="bottom":
                            l=lineObj(self.drawable,self.gc,[self.connectObject.centreCoords,self.connectObject.getHandlePos(self.connectObject.selected)])#(int(e.x),int(e.y))])
                            l.drawTmp()
                            self.connectObject.currentLine=l#a line begins at the bottom - ie data flows out of the object.
                        else:
                            self.connectObject=None
            else:#continue with a connection...
                #print "todo continue connection",self.connectObject.currentLine.coordList
                #check to see if connected to another object - if so, end connection.
                tmp=self.getobject(e.x,e.y,remove=0,name="sci")
                if tmp!=None:
                    if tmp.name=="sci" and tmp.tagID!=self.connectObject.tagID:
                        #have selected top or bottom of object - end connection
                        l=self.connectObject.currentLine
                        l.eraseTmp()
                        l.pop()
                        if self.connectObject.selected=="top":
                            tmp.lines.append(l)
                            tmp.connectto.append(self.connectObject.tagID)
                            self.connectObject.endlines.append(l)
                            self.connectObject.connectfrom.append(tmp.tagID)
                            l.append(tmp.getHandlePos("bottom"))
                            l.coordList.reverse()
                        else:
                            tmp.endlines.append(l)
                            tmp.connectfrom.append(self.connectObject.tagID)
                            self.connectObject.lines.append(l)
                            self.connectObject.connectto.append(tmp.tagID)
                            l.append(tmp.getHandlePos("top"))
                        self.drawlist.append(l)
                        l.draw()
                        self.selobj={"obj":l}
                        self.connectObject=None
                        self.modified(1)
                        self.daexpose()
                    else:
                        print "Invalid connection"
                        l=self.connectObject.currentLine
                        l.eraseTmp()
                        self.connectObject=None
                else:
                    l=self.connectObject.currentLine
                    l.eraseTmp()
                    l.append(l.coordList[-1])
                    l.drawTmp()
        elif self.doing=="share":
            print "sharing..."
            if self.connectObject==None:#starting new share...
                tmp=self.getobject(e.x,e.y,remove=0,name="sci")
                if tmp!=None:
                    if tmp.name=="sci":
                        tmp.selected="main"
                        tmp.centreCoords=tmp.getHandlePos("main")
                        self.connectObject=tmp
                        l=lineObj(self.drawable,self.drawable.new_gc(self.drawable.get_colormap().alloc_color("green")),[self.connectObject.centreCoords,self.connectObject.getHandlePos(self.connectObject.selected)])#(int(e.x),int(e.y))])
                        l.drawTmp()
                        self.connectObject.currentLine=l#a line ends at the top, ie dataflow ends at the input of an object.
            else:#finish the connection
                tmp=self.getobject(e.x,e.y,remove=0,name="sci")
                if tmp!=None and tmp.name=="sci" and self.connectObject.tagID!=tmp.tagID and tmp.imp==self.connectObject.imp and tmp.object==self.connectObject.object and tmp.rank==self.connectObject.rank:#end the resource sharing connection.. to another object that is valid...
                    l=self.connectObject.currentLine
                    l.eraseTmp()
                    l.pop()
                    tmp.sharedFrom.append((self.connectObject.tagID,l))
                    self.connectObject.sharedTo.append((tmp.tagID,l))
                    l.append(tmp.getHandlePos("main"))
                    self.drawlist.append(l)
                    l.draw()
                    self.selobj={"obj":l}
                    self.connectObject=None
                    self.modified(1)
                    self.daexpose()
                else:
                    print "Invalid share"
                    l=self.connectObject.currentLine
                    l.eraseTmp()
                    self.connectObject=None
        elif self.doing=="group":
            l=boxObj(self.drawable,self.drawable.new_gc(self.drawable.get_colormap().alloc_color("blue")),[(int(e.x-50),int(e.y-50)),(int(e.x+50),int(e.y-50)),(int(e.x+50),int(e.y+50)),(int(e.x-50),int(e.y+50))])
            self.drawlist.append(l)
            self.grouplist.append(l)
            l.draw()
            self.selobj={"obj":l}
            self.connectObject=None
            self.modified(1)
            self.daexpose()
                            

                        
                        
                        
        return False
    def dakey(self,w,a=None):
        print "dakey"
    def damotion(self,w,e=None):
        if self.buttonpressed:
            if self.doing=="select":
                if self.selobj!=None:
                    self.modified(1)
                    if self.selobj["obj"].name=="sci":
                        if self.drawMode!="full":
                            self.drawable.draw_image(self.gc,self.selobj["blit"],0,0,self.selobj["obj"].x,self.selobj["obj"].y,-1,-1)
                        # place the object at new coords...
                        self.selobj["obj"].x=int(e.x-self.selobj["offset"][0])
                        if self.selobj["obj"].x<0: self.selobj["obj"].x=0
                        self.selobj["obj"].y=int(e.y-self.selobj["offset"][1])
                        if self.selobj["obj"].y<0: self.selobj["obj"].y=0
                        # now get the new blit image.
                        self.selobj["blit"]=self.drawable.get_image(*self.selobj["obj"].getRect())#copy the image.
                        # now redraw the object.
                        if self.drawMode!="full":
                            self.selobj["obj"].draw()
                        #self.drawable.draw_rectangle(self.gc,True,self.selobj["obj"].x,self.selobj["obj"].y,self.selobj["obj"].w,self.selobj["obj"].h)
                        self.selobj["coords"]=(e.x,e.y)
                    elif self.selobj["obj"].name in ["line","box"]:
                        self.selobj["obj"].moveSelectedPoint(e.x-self.selobj["coords"][0],e.y-self.selobj["coords"][1])
                        self.selobj["coords"]=(e.x,e.y)
                    if e.x>=0:
                        if e.x+30>self.canvsizex:
                            # make canvas larger.
                            self.canvsizex=e.x+30
                            self.adjustCanvasSize()
                            #print "make canvas larger"
                            self.vpx.set_value(self.canvsizex-self.vpx.page_size)
                            self.daexpose()
                        elif e.x<self.vpx.value+30:
                            #print "move scrollbar left"
                            new=e.x-30
                            if new<0:
                                new=0
                            self.vpx.set_value(new)
                            self.daexpose()
                        elif e.x+30>self.vpx.value+self.vpx.page_size:
                            #print "move scrollbar right"
                            new=e.x+30-self.vpx.page_size
                            if new>self.vpx.upper-self.vpx.page_size:
                                new=self.vpx.upper-self.vpx.page_size
                            self.vpx.set_value(new)
                            self.daexpose()
                    else:#make canv larger in -ve direction.
                        self.canvsizex+=30
                        self.shiftAll(30,0)
                        self.adjustCanvasSize()
                        self.vpx.set_value(0)
                        self.daexpose()
                    if e.y>=0:
                        if e.y+30>self.canvsizey:
                            # make canvas larger.
                            self.canvsizey=e.y+30
                            self.adjustCanvasSize()
                            self.vpy.set_value(self.canvsizey-self.vpy.page_size)
                            self.daexpose()
                        elif e.y<self.vpy.value+30:
                            new=e.y-30
                            if new<0:
                                new=0
                            self.vpy.set_value(new)
                            self.daexpose()
                        elif e.y+30>self.vpy.value+self.vpy.page_size:
                            new=e.y+30-self.vpy.page_size
                            if new>self.vpy.upper-self.vpy.page_size:
                                new=self.vpy.upper-self.vpy.page_size
                            self.vpy.set_value(new)
                            self.daexpose()
                    else:#make canv larger in -ve direction.
                        self.canvsizey+=30
                        self.shiftAll(0,30)
                        self.adjustCanvasSize()
                        self.vpy.set_value(0)
                        self.daexpose()
                    if self.drawMode=="full":
                        self.daexpose()
            elif self.doing=="connect":
                if self.connectObject!=None:
                    l=self.connectObject.currentLine
                    l.eraseTmp()
                    l.pop()
                    l.append((int(e.x),int(e.y)))
                    l.drawTmp()
            elif self.doing=="share":
                if self.connectObject!=None:
                    l=self.connectObject.currentLine
                    l.eraseTmp()
                    l.pop()
                    l.append((int(e.x),int(e.y)))
                    l.drawTmp()
            elif self.doing=="group":
                print "DOING group"
        else:
            if self.doing=="connect":
                if self.connectObject!=None:
                    l=self.connectObject.currentLine
                    l.eraseTmp()
                    l.pop()
                    l.append((int(e.x),int(e.y)))
                    l.drawTmp()
            elif self.doing=="share":
                if self.connectObject!=None:
                    l=self.connectObject.currentLine
                    l.eraseTmp()
                    l.pop()
                    l.append((int(e.x),int(e.y)))
                    l.drawTmp()

                
    def daenter(self,w,a=None):
        #print "daenter"
        self.da.window.set_cursor(self.cursors[self.doing])#change cursor
        self.drawable.focus()
    def daleave(self,w,a=None):
        #print "daleave"
        #self.buttonpressed=0
        pass
    def dabuttonreleased(self,w,e=None):
        #print "buttonrelease",w,e
        self.buttonpressed=0
        if self.doing=="object":#place an object
            #print "obj"
            #if e.x==self.buttonpresspos[0] and e.y==self.buttonpresspos[1]:

            #ow,oh=self.selobj["obj"].draw()#drawobject(e.x,e.y)
            self.modified(1)
            self.newObj.setCentre(e.x,e.y)
            self.drawlist.append(self.newObj)
            ow,oh=self.newObj.draw()
            rd=0
            if self.canvsizex<e.x+ow:
                self.canvsizex=e.x+ow
                rd=1
            if self.canvsizey<e.y+oh:
                self.canvsizey=e.y+oh
                rd=1
            if rd:
                self.adjustCanvasSize()
            self.doing="object"
            self.doParamGui(self.newObj)
            self.newObj=self.createNewSimObj()
            self.selobj={"obj":self.drawlist[-1]}
            #self.updateSelectedObject()
            self.daexpose()
            #self.newObj=simObj(self.drawable,self.gc,"sci",colour="red",pango_context=self.da.create_pango_context())
        elif self.doing=="delete":
            if e.x==self.buttonpresspos[0] and e.y==self.buttonpresspos[1]:
                obj=self.getobject(e.x,e.y)
                if obj!=None:
                    self.modified(1)
                    obj.remove()
                    self.daexpose(None)
        elif self.doing=="select":
            if self.selobj!=None and self.selobj["obj"]!=None:
                if self.drawMode!="full":
                    self.drawlist.append(self.selobj["obj"])
                #self.selobj=None
    def dadragbegin(self,w,e=None):
        print "dragbegin"
    def eventReleaseTxt(self,wid,d=None):
        """Detatch the finalisation text from the gui..."""
        #txtw=self.gladetree.get_widget("frameFinalTxt")
        txtw=wid#.get_parent().get_parent()
        txtw.get_parent().remove(txtw)#detatch from parent...
        #create new window to put it in...
        w=gtk.Window()
        w.connect("delete-event",self.closeFinalTxtWindow,txtw)
        w.set_title("Freed text...")
        w.set_transient_for(self.gladetree.get_widget("window2"))
        w.add(txtw)
        w.show_all()
    def closeFinalTxtWindow(self,w,e,txtw):
        #txtw=self.gladetree.get_widget("frameFinalTxt")
        txtw.get_parent().remove(txtw)
        self.gladetree.get_widget("vbox2").pack_start(txtw,expand=False,fill=False)
        txtw.show_all()
        w.destroy()
    def childDetached(self,h,w,d=None):
        """handle box detatched"""
        w.set_size_request(400,400)
    def createNewSimObj(self):
        d=self.simButtonList[self.selectedSimObj]
        dcpu=self.cpuList[self.selectedCPU]
        rank=(int(dcpu["pe"]),int(dcpu["processor"]))
        imp=d["import"]
        object=d["object"]
        feedback=int(d["givesFeedback"])
        self.currTagID+=1
        os.chdir(self.filepath)
        obj=simObj(self.drawable,self.gc,"sci",self.currTagID,msg=d["shortname"],textcol=d["textcol"],pango_context=self.da.create_pango_context(),pixmap=d["pixmap"],colour=dcpu["colour"],rank=rank,imp=imp,object=object,feedback=feedback)
        os.chdir(self.cwd)
        return obj
    def isInGroup(self,obj):
        gr=None
        x=obj.x+obj.w/2
        y=obj.y+obj.h/2
        for g in self.grouplist:
            if g.isIn(obj):
                gr=g
                break
            #g.calcArea()
            #if x>g.x and x<g.x+g.w and y>g.y and y<g.y+g.h:
            #    gr=g
            #    break
        #print "In group",gr
        return gr
        

        
RED=0
GREEN=1
BLUE=2
coldict={"red":[255,0,0],"blue":[0,0,255],"green":[0,255,0],"orange":[255,255,0],"black":[0,0,0]}
class lineObj:
    """A class to display a polyline."""
    def __init__(self,drawable,gc,coordList):
        self.drawable=drawable
        self.name="line"
        self.gc=gc
        self.selgc=self.drawable.new_gc(self.drawable.get_colormap().alloc_color("green"))
        self.coordList=coordList
        self.x=0
        self.y=0
        self.w=0
        self.h=0
        self.selectedPoint=None
        self.removed=0
        self.pyname=""
        self.idstr=""
    def __repr__(self):
        txt="Line object\n"
        txt+="idstr: %s\n"%str(self.idstr)
        txt+="coordList: %s\n"%str(self.coordList)
        return txt
    def append(self,coord):
        self.coordList.append(coord)
        self.calcArea()
    def pop(self,pos=-1):
        poped=self.coordList.pop(pos)
        self.calcArea()
        return poped
    def calcArea(self):
        xmin,ymin=self.coordList[0]
        xmax,ymax=self.coordList[0]
        for x,y in self.coordList:
            if x<xmin:
                xmin=x
            if x>xmax:
                xmax=x
            if y<ymin:
                ymin=y
            if y>ymax:
                ymax=y
        self.x=xmin
        self.y=ymin
        self.w=xmax-xmin
        self.h=ymax-ymin

    def drawTmp(self):
        self.gc.set_function(gtk.gdk.INVERT)
        self.draw()
        self.gc.set_function(gtk.gdk.COPY)
    def eraseTmp(self):
        self.gc.set_function(gtk.gdk.INVERT)
        self.draw()
        self.gc.set_function(gtk.gdk.COPY)
    def draw(self,rect=None):
        self.drawable.draw_lines(self.gc,self.coordList)
    def drawSelected(self):
        if self.removed==0:
            for x,y in self.coordList:
                self.drawable.draw_rectangle(self.selgc,0, x-2, y-2, 4, 4)
    def atPosition(self,x,y):
        rt=False
        self.selectedPoint=None
        #see if we're at a vertex.
        for i in range(len(self.coordList)):
            xx,yy=self.coordList[i]
            if abs(xx-x)<5 and abs(yy-y)<5:
                self.selectedPoint=i
                rt=True
                break
        if rt==False:#see if we're on the line...
            for i in range(len(self.coordList)-1):
                x1,y1=self.coordList[i]
                x2,y2=self.coordList[i+1]
                if x2==x1:#infinite gradient
                    if abs(x-x2)<2 and y>min(y1,y2) and y<max(y1,y2):
                        rt=True
                        break
                else:
                    m=float(y2-y1)/(x2-x1)
                    c=y1-m*x1
                    ycalc=m*x+c
                    if abs(y-ycalc)<2 and x<=max(x1,x2) and x>=min(x1,x2):
                        rt=True
                        break
        if rt==False:#see if we're on the line in the other direction... (if the line is very steep, the near infinite gradient can cause problems...).
            for i in range(len(self.coordList)-1):
                y1,x1=self.coordList[i]
                y2,x2=self.coordList[i+1]
                if x2==x1:#flat line...
                    if abs(y-x2)<2 and x>min(y1,y2) and x<max(y1,y2):
                        rt=True
                        break
                else:
                    m=float(y2-y1)/(x2-x1)
                    c=y1-m*x1
                    xcalc=m*y+c
                    if abs(x-xcalc)<2 and y<=max(x1,x2) and y>=min(x1,x2):
                        rt=True
                        break
                    
                    
        return rt
    def shift(self,x,y):
        nc=[]
        for coord in self.coordList:
            nc.append((coord[0]+x,coord[1]+y))
        self.coordList=nc
        self.x+=x
        self.y+=y
    def getRect(self):
        return [self.x,self.y,self.w,self.h]
    def remove(self):
        self.removed=1
    def moveSelectedPoint(self,x,y):
        """move a vertex by x,y"""
        if self.selectedPoint!=None:
            self.coordList[self.selectedPoint]=(int(self.coordList[self.selectedPoint][0]+x),int(self.coordList[self.selectedPoint][1]+y))
            self.calcArea()
    def setText(self):
        pass

class boxObj(lineObj):
    """This object extends the lineObj object."""
    
    def __init__(self,drawable,gc,coordList):
        lineObj.__init__(self,drawable,gc,coordList)
        self.name="box"
        self.cpuList="[]"
    def __repr__(self):
        txt="""<groupShare idstr="%s" cpu="%s" coordlist="%s"/>\n"""%(str(self.idstr),str(self.cpuList),str(self.coordList))
        #txt+="idstr: %s\n"%str(self.idstr)
        #txt+="coordList: %s\n"%str(self.coordList)
        return txt
    def append(self,coord):
        raise Exception("Box objects can't be appended")
    def pop(self,pos=-1):
        raise Exception("Box objects can't be popped")
    def drawTmp(self):
        raise Exception("Box objects don't have drawTmp")
    def eraseTmp(self):
        raise Exception("Box objects don't have eraseTmp")
    def draw(self,rect=None):
        #print "draw"
        cl=self.coordList+[self.coordList[0]]
        self.drawable.draw_lines(self.gc,cl)
        #print "draw done"
    def moveSelectedPoint(self,x,y):
        """move a virtex and corresponding sides by x,y"""
        if self.selectedPoint!=None:
            sp=self.selectedPoint
            spp=(sp+1)%4
            spm=(sp-1)%4
            self.coordList[sp]=(int(self.coordList[sp][0]+x),int(self.coordList[sp][1]+y))
            
            if sp%2==0:
                self.coordList[spp]=(int(self.coordList[spp][0]),int(self.coordList[spp][1]+y))
                self.coordList[spm]=(int(self.coordList[spm][0]+x),int(self.coordList[spm][1]))
            else:
                self.coordList[spp]=(int(self.coordList[spp][0]+x),int(self.coordList[spp][1]))
                self.coordList[spm]=(int(self.coordList[spm][0]),int(self.coordList[spm][1]+y))

            self.calcArea()
    def isIn(self,obj):
        """look to see whether an object is within this box"""
        minx=miny=-1
        maxx=maxy=0
        for x,y in self.coordList:
            if x>maxx:
                maxx=x
            if x<minx or minx==-1:
                minx=x
            if y>maxy:
                maxy=y
            if y<miny or miny==-1:
                miny=y
        
        if obj.x>minx and obj.x<maxx and obj.y>miny and obj.y<maxy:#+obj.w/2>self.x and obj.x+obj.w/2<self.x+self.w and obj.y+obj.h/2>self.y and obj.y+obj.h/2<self.y+self.h:
            return True
        else:
            return False

        
class simObj:
    def __init__(self,drawable,gc,name,tagID,x=None,y=None,colour="black",pango_context=None,msg=None,pixmap=None,textcol="white",rank=(1,1),imp=None,object=None,feedback=0):
        self.drawable=drawable
        self.pango_context=pango_context
        self.textlabel=pango.Layout(self.pango_context)
        self.gc=gc
        self.tagID=tagID
        self.selgc=self.drawable.new_gc(self.drawable.get_colormap().alloc_color("green"))
        self.name=name
        self.x=self.y=None
        if x!=None:
            self.x=int(x)
        if y!=None:
            self.y=int(y)
        self.selected="main"
        self.idstr=""
        self.centreCoords=None
        self.colour=colour
        self.textcol=textcol
        self.gcolour=None
        self.colourgc=None
        self.lines=[]#lines which begin here...
        self.endlines=[]#lines which end here...
        self.connectto=[]#tags to which we connect
        self.connectfrom=[]#tags which connect to us
        self.sharedTo=[]#tags, lineobj with which resource sharing is implemented
        self.sharedFrom=[]#tag, lineobj with which resource sharing is implemented
        self.currentLine=None
        self.removed=0
        self.msg=msg
        self.rank=rank
        self.groupShare=0#used if the object is defining a group, and may be set to 1 if this object is to implement resource sharing with others in the group.
        self.imp=imp#eg science.el_dm
        self.object=object#eg dm (appended to imp)
        self.feedback=feedback#eg 1 for a reconstructor, as provides feedback.
        self.precode=""#python code to insert before object created
        self.postcode=""#python for after object creation.
        self.pyname=""
        self.args=""
        self.parentNames=None#only used when opening a xml file...
        if name=="rect":
            self.w=self.h=20
        elif name=="sci":
            self.w=69
            self.h=69
            if self.colour==None:
                self.colour="black"
            self.setText()

        if self.x!=None:self.x-=self.w/2
        if self.y!=None:self.y-=self.h/2
        self.setColour(self.colour)
        if self.textcol!=None:
            self.textgc=self.drawable.get_colormap().alloc_color(self.textcol)
            self.textgc=self.drawable.new_gc(self.textgc)
        self.pixmap=pixmap
        if type(self.pixmap)==type(""):
            self.pixmapfile=self.pixmap
            try:
                self.pixmap,mask=gtk.gdk.pixmap_create_from_xpm(self.drawable,None,self.pixmapfile)
            except:
                self.pixmap=None
            if self.pixmap==None:
                try:
                    self.pixmap,mask=gtk.gdk.pixmap_create_from_xpm(self.drawable,None,self.filepath+self.pixmapfile)
                    self.pixmapfile=self.filepath+self.pixmapfile
                except:
                    self.pixmap=None

                    
            if self.pixmap!=None:
                self.pixmapx,self.pixmapy=self.pixmap.get_size()
                if self.pixmapx>self.w-6:
                    self.pixmapx=self.w-6
                if self.pixmapy>self.h-16:
                    self.pixmapy=self.h-16
            #self.pixmapimg=gtk.Image()
            #self.pixmapimg.set_from_file(self.pixmapfile)
            #print dir(self.pixmapimg)
            #self.pixmap=self.pixmapimg.get_image()[0]
            #print dir(self.pixmap)
            #print self.pixmap.get_size()

    def __repr__(self):
        txt='<simulationObject cpu="%s" import="%s" object="%s" pos="%d,%d" tag="%s" shortname="%s" pixmap="%s" feedback="%d" pyname="%s" groupshare="%d" args="%s" connectto="%s" connectfrom="%s" textcol="%s" idstr="%s">\n'%(str(self.rank),self.imp,self.object,self.x+self.w/2,self.y+self.h/2,str(self.tagID),self.msg,self.pixmapfile,self.feedback,self.pyname,self.groupShare,self.args,self.connectto,self.connectfrom,self.textcol,self.idstr)
        #colour depends on rank.  Text colour depends on pixmap and shortname.
        if len(self.precode.strip())>0:
            txt=txt+"<precode>\n%s\n</precode>\n"%self.precode.strip()
        if len(self.postcode.strip())>0:
            txt=txt+"<postcode>\n%s\n</postcode>\n"%self.postcode.strip()
            
        txt=txt+"<lines>\n[\n"
        for line in self.lines:
            txt=txt+str(line.coordList[1:-1])+",\n"
        txt=txt+"]\n</lines>\n"
        txt=txt+"<endlines>\n[\n"
        for line in self.endlines:
            txt=txt+str(line.coordList[1:-1])+",\n"
        txt=txt+"]\n</endlines>\n"
        l=[]
        for line in self.endlines:
            if len(line.idstr)>0:
                l.append(line.idstr)
            else:
                l.append(line.pyname)
        txt=txt+"<parentNames>\n"+str(l)+"\n</parentNames>\n"
        tlist=[]
        for t,l in self.sharedTo:
            tlist.append(t)
        txt=txt+"<sharedTo>\n"+str(tlist)+"\n</sharedTo>\n"
        #l=[]
        #for t,l in self.sharedFrom:
        #    l.append(t)
        #txt=txt+"<sharedFrom>\n"+str(l)+"\n"</sharedFrom>\n"
        txt=txt+"</simulationObject>\n"
        return txt

    def setColour(self,col):
        self.colour=col
        if self.colour!=None:
            self.gcolour=self.drawable.get_colormap().alloc_color(self.colour)
            self.colourgc=self.drawable.new_gc(self.gcolour)

    def getIDStr(self):
        idstr=self.idstr
        args=None
        if len(self.args)>0:
            try:
                args=eval(self.args)
            except NameError:#maybe just the id string... this is no longer supported.
                print "nameerror - cannot evaluate args - (DEPRECIATION WARNING - single idstr no longer supported here)"
                #idstr=self.args
            except:#just ignore rest.
                print "exception %s"%str(self.args)
                #idstr=self.args
        else:
            idstr=self.idstr
        if type(args)==type({}):
            if args.has_key("idstr"):
                idstr=args["idstr"]
                print "DEPRECIATION WARNING - idstr no longer advised in args"
        elif args!=None and idstr=="":
            idstr=args
        return str(idstr)
    def setText(self):
        idstr=self.getIDStr()
        #self.dispmsg="%s\n%s\n%s"%(self.msg,str(self.rank),idstr[:20])
        msg=self.msg
        self.textlabel.set_text(msg)
        w,h=self.textlabel.get_pixel_size()
        while w>self.w-2:
            msg=msg[:-1]
            self.textlabel.set_text(msg)
            w,h=self.textlabel.get_pixel_size()
        msg+="\n%s\n%s"%(str(self.rank),idstr)
        self.textlabel.set_text(msg)
        w,h=self.textlabel.get_pixel_size()
        while w>self.w-2:
            msg=msg[:-1]
            self.textlabel.set_text(msg)
            w,h=self.textlabel.get_pixel_size()
        
    def setCentre(self,x,y):
        self.x=int(x-self.w/2)
        self.y=int(y-self.h/2)
        if self.x<0:self.x=0
        if self.y<0:self.y=0
    def getHandlePos(self,name="main"):
        """get central coords of the object..."""
        coords=None
        if self.name=="rect":
            coords=(int(self.x+self.w/2),int(self.y+self.h/2))
        elif self.name=="sci":
            if name=="main":
                coords=(int(self.x+self.w/2),int(self.y+self.h/2))
            elif name=="top":
                coords=(int(self.x+self.w/2),int(self.y+2))
            elif name=="bottom":
                coords=(int(self.x+self.w/2),int(self.y+self.h-2))
        return coords
    def getRect(self):
        return [self.x,self.y,self.w,self.h]
    def draw(self,rect=None):
        if self.removed==0:
            if self.name=="rect":
                rt=self.drawRect(rect)
            elif self.name=="sci":
                rt=self.drawSci(rect)
        self.adjustLines()
        return rt
    def drawSelected(self):
        if self.removed==0:
            self.drawable.draw_rectangle(self.selgc,0, self.x, self.y, 4, 4)
            self.drawable.draw_rectangle(self.selgc,0, self.x+self.w-4, self.y, 4, 4)
            self.drawable.draw_rectangle(self.selgc,0, self.x, self.y+self.h-4, 4, 4)
            self.drawable.draw_rectangle(self.selgc,0, self.x+self.w-4, self.y+self.h-4, 4, 4)

        
    def adjustLines(self):
        remlist=[]
        oklines=[]
        okto=[]
        for i in range(len(self.lines)):
            l=self.lines[i]
            l.coordList[0]=self.getHandlePos("bottom")
            if l.removed:
                remlist.append(i)
            else:
                oklines.append(l)
                okto.append(self.connectto[i])
        #for i in remlist:
        #    self.lines.pop(i)
        #    self.connectto.pop(i)
        self.lines=oklines
        self.connectto=okto
        remlist=[]
        oklines=[]
        okfrom=[]
        for i in range(len(self.endlines)):
            l=self.endlines[i]
            if type(l)==type([]):
                print "ERROR - endlines"
            else:
                l.coordList[-1]=self.getHandlePos("top")
                if l.removed:
                    remlist.append(i)
                else:
                    oklines.append(l)
                    okfrom.append(self.connectfrom[i])
        self.endlines=oklines
        self.connectfrom=okfrom

        remlist=[]
        oklines=[]
        for i in range(len(self.sharedTo)):
            t,l=self.sharedTo[i]
            if l.removed:
                remlist.append((t,l))
            else:
                oklines.append((t,l))
        self.sharedTo=oklines
        remlist=[]
        oklines=[]
        for i in range(len(self.sharedFrom)):
            t,l=self.sharedFrom[i]
            if l.removed:
                remlist.append((t,l))
            else:
                oklines.append((t,l))
        self.sharedFrom=oklines


        for i in range(len(self.sharedTo)):
            self.sharedTo[i][1].coordList[0]=self.getHandlePos("main")
        for i in range(len(self.sharedFrom)):
            self.sharedFrom[i][1].coordList[-1]=self.getHandlePos("main")
        #for i in remlist:
        #    self.endlines.pop(i)
        #    self.connectfrom.pop(i)

##         for l in self.lines:
##             #l.eraseTmp()
##             l.coordList[0]=self.getHandlePos("bottom")
##             if l.removed:
##                 remlist.append(l)
##             #l.draw()
##         for l in remlist:
##             self.lines.remove(l)
##         remlist=[]
##         for l in self.endlines:
##             #l.eraseTmp()
##             l.coordList[-1]=self.getHandlePos("top")
##             if l.removed:
##                 remlist.append(l)
##             #l.draw()
##         for l in remlist:
##             self.endlines.remove(l)
    def shift(self,x,y):
        self.x+=x
        self.y+=y
    def remove(self):
        self.removed=1
        for l in self.lines:
            l.remove()
        for l in self.endlines:
            l.remove()
        for t,l in self.sharedTo:
            l.remove()
        for t,l in self.sharedFrom:
            l.remove()
        self.lines=[]
        self.endlines=[]
        self.connectfrom=[]
        self.connectto=[]
        self.sharedFrom=[]
        self.sharedTo=[]
    def atPosition(self,x,y):
        if self.name=="rect":
            if self.x<x and self.x+self.w>x and self.y<y and self.y+self.h>y:
                return True
        elif self.name=="sci":
            if self.x<x and self.x+self.w>x and self.y+5<y and self.y+self.h-10>y:
                self.selected="main"
                return True
            elif self.x+self.w/2-2<x and self.x+self.w/2+3>x and self.y<y and self.y+5>y:
                self.selected="top"
                self.centreCoords=(int(self.x+self.w/2),int(self.y+2))
                return True
            elif self.x+self.w/2-2<x and self.x+self.w/2+3>x and self.y+self.h-5<y and self.y+self.h>y:
                self.selected="bottom"
                self.centreCoords=(int(self.x+self.w/2),int(self.y+self.h-2))
                return True
        return False

    def overlapRects(self,r1,r2):
        """Overlaps rectangles, and returns rectangle where they overlap
        """
        if r2==None:
            return r1
        X,Y,W,H=r1
        x,y,w,h=r2
        
        xn=yn=wn=hn=None
        if (X<x and X+W>=x):#start before and ends in or past redraw area
            xn=x#start position
            if X+W>=x+w:#spans redraw area
                wn=w
            else:#ends in redraw area
                wn=X+W-x
        elif X>=x and X<x+w:#starts in redraw area
            xn=X
            if X+W<x+w:#ends in redraw area
                wn=W
            else:#ends outside redraw area
                wn=x+w-X
        if (Y<y and Y+H>=y):#start before and ends in or past redraw area
            yn=y#start position
            if Y+H>=y+h:#spans redraw area
                hn=h
            else:#ends in redraw area
                hn=Y+H-y
        elif Y>=y and Y<y+h:#starts in redraw area
            yn=Y
            if Y+H<y+h:#ends in redraw area
                hn=H
            else:#ends outside redraw area
                hn=y+h-Y
        return xn,yn,wn,hn
            

    def drawSci(self,rect=None):
        if rect==None:
            self.drawable.draw_rectangle(self.gc,True,self.x,self.y+5,self.w,self.h-10)
            self.drawable.draw_rectangle(self.gc,True,self.x+self.w/2-2,self.y,5,5)
            self.drawable.draw_rectangle(self.gc,True,self.x+self.w/2-2,self.y+self.h-5,5,5)
            #tmp=self.gc.foreground
            #self.gc.foreground=self.gcolour
            self.drawable.draw_rectangle(self.colourgc,True,self.x+1,self.y+6,self.w-2,self.h-12)
            #self.gc.foreground=tmp
            if self.pixmap!=None and type(self.pixmap)!=type(""):
                self.drawable.draw_drawable(self.gc,self.pixmap,0,0,self.x+3,self.y+8,self.pixmapx,self.pixmapy)
            w,h=self.textlabel.get_pixel_size()
            self.drawable.draw_layout(self.textgc,self.x+self.w/2-w/2,self.y+self.h/2-h/2,self.textlabel)
        else:
            x,y,w,h=self.overlapRects((self.x,self.y+5,self.w,self.h-10),rect)
            if x!=None and y!=None:
                self.drawable.draw_rectangle(self.gc,True,x,y,w,h)
            x,y,w,h=self.overlapRects((self.x+self.w/2-2,self.y,5,5),rect)
            if x!=None and y!=None:
                self.drawable.draw_rectangle(self.gc,True,x,y,w,h)
            x,y,w,h=self.overlapRects((self.x+self.w/2-2,self.y+self.h-5,5,5),rect)
            if x!=None and y!=None:
                self.drawable.draw_rectangle(self.gc,True,x,y,w,h)
            x,y,w,h=self.overlapRects((self.x+1,self.y+6,self.w-2,self.h-12),rect)
            if x!=None and y!=None:
                self.drawable.draw_rectangle(self.colourgc,True,x,y,w,h)
                if self.pixmap!=None and type(self.pixmap)!=type(""):
                    self.drawable.draw_drawable(self.gc,self.pixmap,0,0,self.x+3,self.y+8,self.pixmapx,self.pixmapy)

            ww,hh=self.textlabel.get_pixel_size()
            x,y,w,h=self.overlapRects((self.x+self.w/2-ww/2,self.y+self.h/2-hh/2,ww,hh),rect)
            if x!=None and y!=None:
                self.drawable.draw_layout(self.textgc,self.x+self.w/2-ww/2,self.y+self.h/2-hh/2,self.textlabel)
        return self.w,self.h
    def drawRect(self,rect=None):
        if rect==None:
            self.drawable.draw_rectangle(self.gc,True,self.x,self.y,self.w,self.h)
            # self.drawlist.append(simObj("rect",x-10,y-10,20,20))
        else: #draw any part of the object that lies within rect.
            xn,yn,wn,hn=self.overlapRects((self.x,self.y,self.w,self.h),rect)
##             x,y,w,h=rect
##             xn=yn=wn=hn=None
##             if (self.x<x and self.x+self.w>=x):#start before and ends in or past redraw area
##                 xn=x#start position
##                 if self.x+self.w>=x1:#spans redraw area
##                     wn=w
##                 else:#ends in redraw area
##                     wn=self.x+self.w-x
##             elif self.x>=x and self.x<x+w:#starts in redraw area
##                 xn=self.x
##                 if self.x+self.w<x+w:#ends in redraw area
##                     wn=self.w
##                 else:#ends outside redraw area
##                     wn=x+w-self.x
##             if (self.y<y and self.y+self.h>=y):#start before and ends in or past redraw area
##                 yn=y#start position
##                 if self.y+self.h>=y1:#spans redraw area
##                     hn=h
##                 else:#ends in redraw area
##                     hn=self.y+self.h-y
##             elif self.y>=y and self.y<y+h:#starts in redraw area
##                 yn=self.y
##                 if self.y+self.h<y+h:#ends in redraw area
##                     hn=self.h
##                 else:#ends outside redraw area
##                     hn=y+h-self.y
            if xn!=None and yn!=None:
                self.drawable.draw_rectangle(self.gc,True,xn,yn,wn,hn)

            
        return (self.w,self.h)#return object width, height.
        


                



def run():
    d=simsetup()
    gtk.main()

if __name__=="__main__":
    run()
