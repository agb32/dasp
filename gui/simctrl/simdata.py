import xml.parsers.expat,types,string
## import Numeric,time,os,numpy # commented out by UB on 2012Jun07 to remove Numeric
import time,os,numpy # NEW by UB on 2012Jun07 to remove Numeric
try:
    import gist
except:
    print "simdata - cannot import gist - ignoring"
from gui.textbox.textbox import textbox
from gui.pylab import mypylab
#import util.pyfits,numarray
class simdata:
    """A class for reading a simulation data xml file.  This file is used
    to determine what plots etc are available for the gui."""
    def __init__(self,file=None,xmltxt=None,cancelData=None):
        """Initialise.
        canceldata is the function that can be called to cancel a request
        (eg if a plot window is closed) - typically will be sim.execute
        function.
        """
        self.p=None
        self.cancelData=cancelData
        self.filename=file
        self.reset()
        if file!=None:
            self.open(file)
        elif xmltxt!=None:
            self.open(None,xmltxt)

    def reset(self):
        """Reset (ready for a load)"""
        self.p=xml.parsers.expat.ParserCreate()
        self.p.StartElementHandler=self.start_element
        self.p.EndElementHandler = self.end_element
        self.p.CharacterDataHandler = self.char_data
        self.p.returns_unicode=0
        self.dataList=Plot()#No longer a list, but a set of nodes for the tree...
        self.dataList.title="base"
        self.doingPlot=None
        self.infeedbackmsg=0
        self.inpost=0
        self.incmd=0
        self.nodeList=[self.dataList]
    def open(self,file,txt=None):
        """open and read file"""
        if file!=None:
            txt=open(file).read()
            self.filename=file
        if txt!=None:
            self.p.Parse(txt)
    def start_element(self,name,attrs):
        """start of XML tag"""
        if name=="plot":
            plot=Plot()
            self.doingPlot=plot
            self.nodeList[-1].childList.append(plot)
            if attrs.has_key("title"):
                plot.title=attrs["title"]
            if attrs.has_key("cmd"):
                plot.cmd=attrs["cmd"]
            if attrs.has_key("ret"):
                plot.ret=attrs["ret"]
            if attrs.has_key("type"):
                if attrs["type"]=="gist":
                    plot.gisttype=1
                elif attrs["type"]=="pylab":
                    plot.pylabtype=1
                elif attrs["type"]=="save":
                    plot.savetype=1
                elif attrs["type"]=="text":
                    plot.texttype=1
                elif attrs["type"]=="feedback":
                    plot.feedbacktype=1
                #plot.type=attrs["type"]
            if attrs.has_key("gisttype"):
                plot.gisttype=int(attrs["gisttype"])
            if attrs.has_key("savetype"):
                plot.savetype=int(attrs["savetype"])
            if attrs.has_key("texttype"):
                plot.texttype=int(attrs["texttype"])
            if attrs.has_key("pylabtype"):
                plot.pylabtype=int(attrs["pylabtype"])
            if attrs.has_key("feedbacktype"):
                plot.feedbacktype=int(attrs["feedbacktype"])
            if attrs.has_key("preprocess"):
                plot.preprocess=attrs["preprocess"]
            if attrs.has_key("post"):
                plot.post=attrs["post"]
            if attrs.has_key("dim"):
                plot.dim=int(attrs["dim"])
            if attrs.has_key("xaxis"):
                plot.xaxis=attrs["xaxis"]
            if attrs.has_key("when"):
                plot.when=attrs["when"]
            if plot.gisttype:#now investigate params specific for gist.
                if attrs.has_key("window"):
                    plot.info["window"]=int(attrs["window"])
                else:
                    plot.info["window"]=0
                if attrs.has_key("palette"):
                    plot.info["palette"]=attrs["palette"]
                else:
                    plot.info["palette"]="gray.gp"
                if attrs.has_key("gistdpi"):
                    plot.info["gistdpi"]=int(attrs["gistdpi"])
                else:
                    plot.info["gistdpi"]=75
            if plot.texttype:#specific params for text
                if attrs.has_key("wintype"):
                    plot.info["wintype"]=attrs["wintype"]
                else:
                    plot.info["wintype"]="ownwindow"
                if attrs.has_key("textreplace"):
                    plot.info["textreplace"]=int(attrs["textreplace"])
                else:
                    plot.info["textreplace"]=0
            if plot.savetype:
                plot.info["filetype"]="fits"
                plot.info["filename"]="tmp.fits"
                plot.info["filereplace"]=0
                if attrs.has_key("filetype"):
                    plot.info["filetype"]=attrs["filetype"]
                if attrs.has_key("filename"):
                    plot.info["filename"]=attrs["filename"]
                if attrs.has_key("replace"):
                    plot.info["filereplace"]=int(attrs["filereplace"])
            if plot.feedbacktype:
                plot.info["feedbackmsg"]="msg=feedback"
                if attrs.has_key("feedbackmsg"):
                    plot.info["feedbackmsg"]=attrs["feedback"]
                
            if plot.title=="":
                plot.title=plot.cmd
            plot.cancelData=self.cancelData
            #self.dataList.append(plot)
        elif name=="feedbackmsg":
            self.infeedbackmsg=1
        elif name=="cmd":
            if self.doingPlot!=None:
                self.incmd=1
        elif name=="post":
            if self.doingPlot!=None:
                self.inpost=1
            else:
                print "ERROR during parsing - post tag must be inside plot tag"
                raise Exception("Error during parsing - post tag must be inside plot tag")
        elif name=="simobj":
            plot=Plot()
            self.nodeList[-1].childList.append(plot)

            if attrs.has_key("name"):
                plot.title=attrs["name"]
            else:
                plot.title="unknown"
            self.nodeList.append(plot)
    def char_data(self,data):
        """middle of XML tag"""
        if self.doingPlot!=None:
            if self.infeedbackmsg==1:
                if self.doingPlot.info.has_key("feedbackmsg"):
                    self.doingPlot.info["feedbackmsg"]+=data
                else:
                    self.doingPlot.info["feedbackmsg"]=data
            elif self.inpost==1:
                self.doingPlot.post+=data
            elif self.incmd==1:
                self.doingPlot.cmd+=data
            else:
                self.doingPlot.preprocess+=data
    def end_element(self,name):
        """end of XML tag"""
        if name=="plot":
            if self.infeedbackmsg==1:
                self.infeedbackmsg=0
                if self.doingPlot.info.has_key("feedbackmsg"):
                    self.doingPlot.info["feedbackmsg"]=string.strip(self.doingPlot.info["feedbackmsg"])
            else:
                self.doingPlot.preprocess=string.strip(self.doingPlot.preprocess)
            self.doingPlot=None
##            if self.dataList[-1].cmd==None:
##                #self.dataList.pop()
##                self.dataList[-1].cmd=""
        elif name=="post":
            if self.inpost:
                self.doingPlot.post=self.doingPlot.post.strip()
                self.inpost=0
        elif name=="cmd":
            if self.incmd:
                self.doingPlot.cmd=self.doingPlot.cmd.strip()
                self.incmd=0
        elif name=="simobj":
            self.nodeList.pop()
class Plot:
    """Class for dealing with data from simulation"""
    def __init__(self):
        """Initialise"""
        self.title=""
        self.cancelData=None#typically, will be a sim.execute function.
        self.cmd=""#None
        self.ret=""#changed from None... 20070112.
        self.tag=None
        self.gisttype=0
        self.savetype=0
        self.texttype=0
        self.pylabtype=0
        self.feedbacktype=0
        self.preprocess=""
        self.post=""
        self.dim=None
        self.xaxis=None
        self.when="cmd"
        self.info={}
        self.repeating=0
        self.textWindow=None
        self.connList=[]
        self.button=None
        self.gistWindow=None
        self.childList=[]
    def gettype(self):
        """returns a string representation of the types"""
        s="save "*self.savetype+"gist "*self.gisttype+"pylab "*self.pylabtype+"text "*self.texttype+"feedback "*self.feedbacktype
        return s[:-1]
    def aslist(self,displayAll):
        if displayAll:
            return [self.title,self.cmd,self.ret,self.gettype(),self.preprocess,self.post,str(self.dim),self.xaxis,self.when,str(self.info)]
        else:
            cmd=self.cmd.replace("\n","\\n").strip()
            if len(cmd)>33:
                cmd=cmd[:30]+"..."
            ret=self.ret.replace("\n","\\n").strip()
            if len(ret)>33:
                ret=ret[:30]+"..."
            typ=self.gettype().replace("\n","\\n").strip()
            if len(typ)>33:
                typ=typ[:30]+"..."
            pre=self.preprocess.replace("\n","\\n").strip()
            if len(pre)>33:
                pre=pre[:30]+"..."
            post=self.post.replace("\n","\\n").strip()
            if len(post)>33:
                post=post[:30]+"..."
            dim=str(self.dim).replace("\n","\\n").strip()
            if len(dim)>33:
                dim=dim[:30]+"..."
            xaxis=str(self.xaxis).replace("\n","\\n").strip()
            if len(xaxis)>33:
                xaxis=xaxis[:30]+"..."
            when=self.when.replace("\n","\\n").strip()
            if len(when)>33:
                when=when[:30]+"..."
            info=str(self.info).replace("\n","\\n").strip()
            if len(info)>33:
                info=info[:30]+"..."
            return [self.title,cmd,ret,typ,pre,post,dim,xaxis,when,info]
            

        
    def __repr__(self):
        s="%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n"%(self.title,str(self.cmd),str(self.ret),self.gettype(),str(self.preprocess),str(self.post),str(self.dim),str(self.xaxis),str(self.when,self.info))
        return s
    def cancel(self):
        """User has decided to cancel feedback from simulation by clicking button"""
        if self.gisttype:
            gist.winkill(self.info["window"])
            self.gistWindow=None
        if self.texttype:
            if self.textWindow!=None:
                if self.textWindow.destroyed==0:
                    self.textWindow.close(None)
                self.textWindow=None
                if self.button!=None:
                    self.button[0]=0
                self.repeating=0
        if self.pylabtype:
            if self.info.has_key("plotwin"):
                self.info["plotwin"].quit()
                if self.info.has_key("plotwin"):
                    del(self.info["plotwin"])
                if self.button!=None:
                    self.button[0]=0
                self.repeating=0
        if self.cancelData!=None:
            self.cancelData(self.cmd,tag=self.tag,action="del",connList=self.connList)

    def deactivate(self):
        """Can be used to deactivate the data..."""
        if self.button!=None:
            self.button[0]=0
        self.repeating=0
        #del(self.info["plotwin"])
    def handle(self,data):
        """Handle data returned from the simulation"""
        if data[1]=="warning" or data[1]=="error":
            print "Error retrieving data",data[1:],data[0]
        elif len(self.gettype())==0:#not doing anything with it...
            pass
        else:
            if self.when[:3]!="rpt" or self.repeating==1:#still expecting data
                if self.ret!=None and data[3].has_key(self.ret):
                    data=data[3][self.ret]
                else:
                    data=None
                ret=self.ret
                if ret==None:
                    ret="None"
                if self.button==None:
                    d={ret:data,"button":None}
                else:
                    d={ret:data,"button":self.button[0]}
                try:
                    exec self.preprocess in d
                except:
                    pass
                if self.button!=None:
                    self.button[0]=d["button"]
                data=d[ret]
                dim=self.dim
                xaxis=self.xaxis
                if dim==None:
##                    if type(data) in [numpy.ndarray,Numeric.ArrayType]: # commented out by UB on 2012Jun07 to remove Numeric
                    if type(data) == numpy.ndarray: # NEW by UB on 2012Jun07 to remove Numeric
                        if len(data.shape)>1:
                            dim=2
                        elif len(data.shape)==1:
                            dim=1
                        else:
                            dim=0
                    else:
                        dim=0
                if dim==1:
                    if type(self.xaxis)==types.NoneType:
                        if len(data.shape)>1:
                            xaxis=data[0]
                            data=data[1:]
                        else:
                            xaxis=numpy.arange(data.shape[0])
                            data=data
                    else:
                        if type(self.xaxis)==types.StringType:
                            xaxis=eval(self.xaxis)
                if self.gisttype:
                    #print "Plotting gist data (if you know more about gist than me, please improve this in gui/simdata.py)"
##                     try:
##                         print data.typecode(),data.shape,"Dims:",dim
##                     except:
##                         pass
##                    if type(data) in [numpy.ndarray,Numeric.ArrayType]: # commented out by UB on 2012Jun07 to remove Numeric
                    if type(data) == numpy.ndarray:  # NEW by UB on 2012Jun07 to remove Numeric
                        if not self.info.has_key("window"):
                            self.info["window"]=0
                        if not self.info.has_key("palette"):
                            self.info["palette"]="gray.gp"
                        if not self.info.has_key("gistdpi"):
                            self.info["gistdpi"]=75
                        if self.gistWindow==None:
                            self.gistWindow=gist.window(self.info["window"],wait=1,dpi=self.info["gistdpi"])
                            gist.animate(0)
                            gist.animate(1)
                            gist.palette(self.info["palette"])
                        else:
                            gist.window(self.gistWindow)
                        #gist.fma()
                        if dim==1:
                            for i in range(data.shape[0]):
                                gist.plg(data[i],xaxis)
                        else:
                            gist.pli(data)
                        gist.fma()
                    else:
                        print "Cannot display type %s with gist"%str(type(data))
                if self.pylabtype:
##                    if type(data)==Numeric.ArrayType or type(data)==numpy.ndarray: # commented out by UB on 2012Jun07 to remove Numeric
                    if type(data)==numpy.ndarray: # NEW by UB on 2012Jun07 to remove Numeric
                        if not self.info.has_key("palette"):
                            self.info["palette"]="gray"
                        if not self.info.has_key("interp"):
                            self.info["interp"]="nearest"
                        if not self.info.has_key("plotwin"):
                            self.info["plotwin"]=mypylab.plot()
                            p=self.info["plotwin"]
                            p.win.set_title(self.title)
                            p.newPalette(self.info["palette"])
                            p.newInterpolation(self.info["interp"])
                            p.deactivatefn=self.cancel#deactivate
                        p=self.info["plotwin"]
                        if dim==1:
                            p.dims=1
                            axis=xaxis
                        else:
                            p.dims=2
                            axis=None
                        if p.active:
                            p.plot(data,axis=axis)
                        else:
                            if self.button!=None:
                                self.button[0]=0
                            #print "Not expecting this data any more... (simdata.handle, type=pylab)"
                            self.repeating=0
                    else:
                        print "Cannot display type %s with pylab"%str(type(data))


                if self.texttype:
                    #self.info["texttype"]=="ownwindow, mainwindow", default own
                    #self.info["replace"]==1 or 0, default 0
                    if not self.info.has_key("wintype"):
                        self.info["wintype"]="ownwindow"
                    if not self.info.has_key("textreplace"):
                        self.info["textreplace"]=0
                    if self.info["wintype"]=="ownwindow":
                        if self.textWindow==None:
                            self.textWindow=textbox(self.title)
                            self.textWindow.closeFunc=self.cancel#deactivate
                        if self.textWindow.destroyed==0:
                            #print "adding text",str(data)
                            self.textWindow.addText(str(data)+"\n",replace=self.info["textreplace"])
                        else:#tell simulation not to send...
                            print "Not expecting this data any more... (simdata.handle, type=text)"
                            self.textWindow=None
                            self.repeating=0
                    else:
                        print str(data)
                if self.savetype:
                    #self.info["filetype"]=="fits, csv, text" default fits
                    #self.info["filename"]
                    #self.info["replace"]==1 or 0, default 0
                    if not self.info.has_key("filetype"):
                        self.info["filetype"]="fits"
                    if not self.info.has_key("filename"):
                        self.info["filename"]="tmp.fits"
                    if not self.info.has_key("filereplace"):
                        self.info["filereplace"]=0
                    if self.info["filetype"]=="fits":
##                        if type(data) in [numpy.ndarray,Numeric.ArrayType]: # commented out by UB on 2012Jun07 to remove Numeric
                        if type(data) == numpy.ndarray: # NEW by UB on 2012Jun07 to remove Numeric
                            print "WARNING - depreciated - use util.FITS instead (code needs updating)"
                            if self.info["filereplace"]:
                                imghdu=util.pyfits.PrimaryHDU(numarray.array(data))
                                imghdu.header.update("DATE",time.asctime())
                                imghdu.header.update("USER",os.environ["USER"])
                                imghdu.header.update("CREATOR","simctrl.py simulation control")
                                imghdu.header.update("TITLE",str(self.title))
                                imghdu.header.update("COMMAND",str(self.cmd))
                                imghdu.header.update("RETURN",str(self.ret))
                                imghdu.header.update("TYPE",self.gettype())
                                imghdu.header.update("PREPROC",str(self.preprocess))
                                imghdu.header.update("DIMS",str(self.dim))
                                imghdu.header.update("XAXIS",str(self.xaxis))
                                imghdu.header.update("WHEN",str(self.when))
                                imghdu.header.update("INFO",str(self.info))
                                hdulist=util.pyfits.HDUList([imghdu])
                                hdulist.writeto(self.info["filename"],clobber=True)
                            else:
                                f=util.pyfits.open(self.info["filename"],mode="update")
                                imghdu=util.pyfits.ImageHDU(numarray.array(data))
                                imghdu.header.update("DATE",time.asctime())
                                imghdu.header.update("USER",os.environ["USER"])
                                imghdu.header.update("CREATOR","simctrl.py simulation control")
                                imghdu.header.update("TITLE",str(self.title))
                                imghdu.header.update("COMMAND",str(self.cmd))
                                imghdu.header.update("RETURN",str(self.ret))
                                imghdu.header.update("TYPE",self.gettype())
                                imghdu.header.update("PREPROC",str(self.preprocess))
                                imghdu.header.update("DIMS",str(self.dim))
                                imghdu.header.update("XAXIS",str(self.xaxis))
                                imghdu.header.update("WHEN",str(self.when))
                                imghdu.header.update("INFO",str(self.info))
                                f.append(imghdu)
                                f.close()
                        else:
                            print "Cannot save fits data of this format:",type(data)
                    elif self.info["filetype"]=="csv":
                        if self.info["filereplace"]:
                            mode="w"
                        else:
                            mode="a"
                        f=open(self.info["filename"],mode)
                        f.write("#Date\t%s\n#User\t%s\n#Creator\tsimctrl.py simulation control\n#Title\t%s\n#Command\t%s\n#Return\t%s\n#Type\t%s\n#Preprocess\t%s\n#Dims\t%s\n#Xaxis\t%s\n#When\t%s\n#Info\t%s\n"%(time.asctime(),os.environ["USER"],str(self.title),str(self.cmd),str(self.ret),self.gettype(),str(self.preprocess),str(self.dim),str(self.xaxis),str(self.when),str(self.info)))
                        if dim==1:
                            try:
                                for i in range(xaxis.shape[0]):
                                    f.write("%g"%float(xaxis[i]))
                                    for j in range(data.shape[0]):
                                        f.write("\t%g"%float(data[j][i]))
                                    f.write("\n")
                                f.write("\n")
                            except:
                                print "Data not in correct 1D format - can't save as csv"
                                f.write(str(data))
                                f.write("\n\n")
                        else:
                            print "Can't save 2D data as csv... using text instead"
                            f.write(str(data))
                            f.write("\n\n")
                        f.close()
                    elif self.info["filetype"]=="text":
                        if self.info["filereplace"]:
                            mode="w"
                        else:
                            mode="a"
                        f=open(self.info["filename"],mode)
                        f.write("#Date\t%s\n#User\t%s\n#Creator\tsimctrl.py simulation control\n#Title\t%s\n#Command\t%s\n#Return\t%s\n#Type\t%s\n#Preprocess\t%s\n#Dims\t%s\n#Xaxis\t%s\n#When\t%s\n#Info\t%s\n"%(time.asctime(),os.environ["USER"],str(self.title),str(self.cmd),str(self.ret),self.gettype(),str(self.preprocess),str(self.dim),str(self.xaxis),str(self.when),str(self.info)))
                        f.write(str(data))
                        f.write("\n\n")
                        f.close()
                    else:
                        print "Unrecognised filetype - not saving"
                if self.feedbacktype:
                    try:
                        d={"feedback":data}
                        exec self.info["feedbackmsg"] in d
                        msg=d["msg"]
                    except:
                        msg="Feedback data:"+str(data)
                    print msg
                exec self.post
                
            else:
                print "Warning: No longer expecting data for",self.cmd
