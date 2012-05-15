"""
show how to add a matplotlib FigureCanvasGTK or FigureCanvasGTKAgg widget and
a toolbar to a gtk.Window
"""

#import Numeric,RandomArray
import numpy,numpy.random
from matplotlib.axes import Subplot
from matplotlib.figure import Figure
from matplotlib.numerix import arange, sin, pi
import thread
# uncomment to select /GTK/GTKAgg/GTKCairo
from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas

# or NavigationToolbarfor classic
from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
#from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

    


from matplotlib.figure import FigureImage
import pylab
import matplotlib.cm as colour
#from matplotlib import interactive
#interactive(True)
import gtk,gobject
from gui.myFileSelection.myFileSelection import myFileSelection



## def mysave(toolbar=None,button=None,c=None):
##         print "mypylabsave"
##         print toolbar.canvas,dir(toolbar.canvas)
##         print type(toolbar.canvas._pixmap)
##         data=toolbar.canvas.get_data()
##         print type(data)
##         fn=myFileSelection("Open an xml parameter file",self.oldfilename,parent=self.gladetree.get_widget("window1")).fname
##         if fn!=None:
##             if fn[-5:] not in [".fits",".FITS"]:
##                 print "Only saving as FITS is supported"
##             else:
##                 util.FITS.Write(self.data,fn)
##                 print "Data saved as %s"%fn

## NavigationToolbar.save_figure=mysave


        
X=numpy.random.random((20,20)).astype("f")
class myToolbar:
    def __init__(self,plotfn=None):
        """plotfn is a function to call to replot..."""
        self.data=None
        if plotfn!=None:
            self.replot=plotfn
        else:
            self.replot=self.dummyreplot
        self.autoscale=1
        self.freeze=0
        self.logx=0
        self.mangleTxt=""
        self.scale=[0,1]
        self.toolbar=gtk.VBox()
        self.hbox=gtk.HBox()
        self.tooltips=gtk.Tooltips()
        self.freezebutton=gtk.CheckButton("Freeze")
        self.freezebutton.set_active(self.freeze)
        self.freezebutton.connect("toggled",self.togglefreeze,None)
        self.tooltips.set_tip(self.freezebutton,"freeze display")
        self.autobutton=gtk.CheckButton("Scaling")
        self.autobutton.set_active(self.autoscale)
        self.autobutton.connect("toggled",self.toggleAuto,None)
        self.tooltips.set_tip(self.autobutton,"autoscale data")
        self.scaleMinEntry=gtk.Entry()
        self.scaleMinEntry.connect("focus-out-event",self.rescale,"min")
        self.scaleMinEntry.connect("activate",self.rescale,"min")
        self.scaleMinEntry.set_width_chars(8)
        self.tooltips.set_tip(self.scaleMinEntry,"Minimum value to clip when not autoscaling")
        self.scaleMaxEntry=gtk.Entry()
        self.scaleMaxEntry.connect("focus-out-event",self.rescale,"max")
        self.scaleMaxEntry.connect("activate",self.rescale,"max")
        self.scaleMaxEntry.set_width_chars(8)
        self.tooltips.set_tip(self.scaleMaxEntry,"Maximum value to clip when not autoscaling")
        self.logxbutton=gtk.CheckButton("Logx")
        self.logxbutton.set_active(self.logx)
        self.logxbutton.connect("toggled",self.togglelogx,None)
        self.tooltips.set_tip(self.logxbutton,"Logaritm of x axis for 1d plots")
        self.dataMangleEntry=gtk.Entry()
        self.dataMangleEntry.connect("focus-out-event",self.dataMangle,None)
        self.dataMangleEntry.connect("activate",self.dataMangle,None)
        self.tooltips.set_tip(self.dataMangleEntry,"Formatting to perform on data prior to plotting, e.g. data=numpy.log(data) (this gets exec'd)")
        self.hbox.pack_start(self.freezebutton)
        self.hbox.pack_start(self.autobutton)
        self.hbox.pack_start(self.scaleMinEntry)
        self.hbox.pack_start(self.scaleMaxEntry)
        self.hbox.pack_start(self.logxbutton)
        self.toolbar.pack_start(self.hbox)
        self.toolbar.pack_start(self.dataMangleEntry)
        self.toolbar.show_all()
    def dummyreplot(self):
        print "Replot data... (doing nowt)"
    def toggleAuto(self,w,data=None):
        self.autoscale=self.autobutton.get_active()
        self.replot()

    def rescale(self,w,e,data=None):
        if data==None:
            data=e
        if data=="min":
            indx=0
        else:
            indx=1
        try:
            self.scale[indx]=float(w.get_text())
        except:
            pass
        self.replot()
                
    def togglefreeze(self,w,data=None):
        self.freeze=self.freezebutton.get_active()
        if not self.freeze:
            self.replot()
    def togglelogx(self,w,data=None):
        self.logx=self.logxbutton.get_active()
        self.replot()
    def dataMangle(self,w,e=None,data=None):
        txt=w.get_text().strip()
        if self.mangleTxt!=txt:
            self.mangleTxt=txt
            self.replot()
    def prepare(self,data,dim=2):
        if self.freeze==0:
            data=data.copy()
            if len(self.mangleTxt)>0:
                d={"data":data,"numpy":numpy}
                try:
                    exec self.mangleTxt in d
                    data=d["data"]#the new data... after mangling.
                except:
                    pass
            if dim==2:#see if dimensions have changed...
                dim=len(data.shape)
            if self.logx==1 and dim==2:
                #take log of the data...
                m=numpy.min(data.ravel())
                if m<0:
                    data+=0.1-m#now ranges from 0 upwards
                elif m==0.:#get rid of any zeros...
                    data+=0.1
                data=numpy.log(data)
            if self.autoscale:
                tmp=data.flat
                #if data.flags.contiguous:
                #    tmp=data.flat
                #else:
                #    tmp=numpy.array(data).flat
                self.scale[0]=numpy.min(tmp)
                self.scale[1]=numpy.max(tmp)
                self.scaleMinEntry.set_text("%.4g"%(self.scale[0]))
                self.scaleMaxEntry.set_text("%.4g"%(self.scale[1]))
            if dim==1:
                data=numpy.where(data<self.scale[0],self.scale[0],data)
                data=numpy.where(data>self.scale[1],self.scale[1],data)
        self.data=data
        return self.freeze,self.logx,data,self.scale
##     def mysave(self,toolbar=None,button=None,c=None):
##         print "mypylabsave"
##         print a,b,c
##         fn=myFileSelection("Open an xml parameter file",self.oldfilename,parent=self.gladetree.get_widget("window1")).fname
##         if fn!=None:
##             if fn[-5:] not in [".fits",".FITS"]:
##                 print "Only saving as FITS is supported"
##             else:
##                 util.FITS.Write(self.data,fn)
##                 print "Data saved as %s"%fn

class plot:
    """Note, currently, this cant be used interactively - because Gtk has to be running...."""
    def __init__(self,window=None,startGtk=0,dims=2):
        self.dims=dims
        self.data=numpy.zeros((10,10),numpy.float32)
        self.data[:]=numpy.arange(10).astype(numpy.float32)
        self.deactivatefn=None#this can be set by the caller, eg to turn off buttons...
        
        self.win = gtk.Window()
        self.win.connect("destroy", self.quit)
        self.win.set_default_size(400,400)
        self.win.set_title("Window")
        self.cmap=colour.gray
        self.vbox = gtk.VBox()
        self.interpolation="nearest"#see pylab documantation for others.
        self.win.add(self.vbox)
        self.vbox.connect("button_press_event",self.buttonPress)
        self.fig = Figure(figsize=(5,4), dpi=50)
        self.ax=self.fig.add_subplot(111)
        self.fig.subplots_adjust(right=0.99,left=0.08,bottom=0.05,top=0.99)
        #self.ax.imshow(self.data,interpolation=self.interpolation)
        #print type(fig),dir(ax),dir(fig)
        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.vbox.pack_start(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self.win)
        self.vbox.pack_start(self.toolbar, False, False)
        self.mytoolbar=myToolbar(plotfn=self.plot)
        #self.toolbar.save_figure=self.mytoolbar.mysave
        self.vbox.pack_start(self.mytoolbar.toolbar,False,False)
        self.win.show_all()
        self.toolbar.hide()
        self.mytoolbar.toolbar.hide()
        self.active=1#will be set to zero once quit or window closed.
        self.toolbarVisible=0
        self.startedGtk=0
        self.update=0
        #self.plot()
        if startGtk==1 and gtk.main_level()==0:
            self.startedGtk=1
            thread.start_new_thread(gtk.main,())
    def quit(self,w=None,data=None):
        if self.deactivatefn!=None:
            d=self.deactivatefn
            self.deactivatefn=None
            d()
        self.active=0
        self.win.hide()
        if self.startedGtk:
            gtk.main_quit()
    def newPalette(self,palette):
        if palette[-3:]==".gp":
            palette=palette[:-3]
        if palette in colour.datad.keys():
            self.cmap=getattr(colour,palette)
        else:
            print "Palette %s not regocnised"%str(palette)
        #self.plot()
    def newInterpolation(self,interp):
        if interp not in ["bicubic","bilinear","blackman100","blackman256","blackman64","nearest","sinc144","sinc64","spline16","spline36"]:
            print "Interpolation %s not recognised"%str(interp)
        else:
            self.interpolation=interp
        #self.plot()
        
    def buttonPress(self,w,e,data=None):
        """If the user right clicks, we show or hide the toolbar..."""
        if e.button==3:
            if self.toolbarVisible:
                self.toolbar.hide()
                self.mytoolbar.toolbar.hide()
                self.toolbarVisible=0
            else:
                self.toolbar.show()
                self.mytoolbar.toolbar.show()
                self.toolbarVisible=1

        return True

    def queuePlot(self,axis):
        """puts a request to plot in the idle loop... (gives the rest of the
        gui a chance to update before plotting)
        """
        #print type(axis),self.data.shape
        #if type(axis)!=type(None):
        #    print axis.shape
        if self.update:
            if hasattr(self.ax.xaxis,"callbacks"):
                self.ax.xaxis.callbacks.callbacks=dict([(s,dict()) for s in self.ax.xaxis.callbacks.signals])#needed to fix a bug!
                self.ax.yaxis.callbacks.callbacks=dict([(s,dict()) for s in self.ax.yaxis.callbacks.signals])#needed to fix a bug!
            
            self.ax.clear()
            freeze,logscale,data,scale=self.mytoolbar.prepare(self.data,dim=self.dims)
            if len(data.shape)==1 or self.dims==1:
                #1D
                if len(data.shape)==1:
                    if freeze==0:
                        if type(axis)==type(None) or axis.shape[0]!=data.shape[0]:
                            axis=numpy.arange(data.shape[0])+1
                        if logscale:
                            try:
                                axis=numpy.log10(axis)
                            except:
                                print "Cannot take log"
                        #self.fig.axis([axis[0],axis[-1],scale[0],scale[1]])
                        #print dir(self.ax)
                        #print dir(self.ax.axis)
                        #print self.ax.get_position()
                        #self.ax.cla()
                        #self.ax.axis([-1,1,-1,1])
                        #self.ax.autoscale_view()
                        try:#older installations don't have this (e.g. cray)
                            self.ax.set_aspect("auto")
                        except:
                            pass
                        self.ax.plot(axis,data)
                else:#use first row of data for the x axis...
                    #axis=data[0]
                    #freeze,logscale,data,scale=self.mytoolbar.prepare(self.data,dim=1)
                    if freeze==0:
                        if type(axis)==type(None) or axis.shape[0]!=data.shape[-1]:
                            axis=numpy.arange(data.shape[0])+1
                        if logscale:
                            try:
                                axis=numpy.log10(axis)
                            except:
                                print "Cannot take log"
                        #self.fig.axis([axis[0],axis[-1],scale[0],scale[1]])
                        for i in range(data.shape[0]):
                            self.ax.plot(axis,data[i])

            else:#2D
                if len(data.shape)!=2:#force to 2d
                    data=numpy.reshape(data,(reduce(lambda x,y:x*y,data.shape[:-1]),data.shape[-1]))
                #freeze,logscale,data,scale=self.mytoolbar.prepare(self.data)
                if freeze==0:
                    self.ax.imshow(data,interpolation=self.interpolation,cmap=self.cmap,vmin=scale[0],vmax=scale[1],origin="lower")
            if freeze==0:
                try:
                    self.ax.draw()
                except:
                    pass
                #self.ax.update()
                self.canvas.draw()
                #self.canvas.queue_draw()
        self.update=0
        return False
    
    def plot(self,data=None,copy=0,axis=None):
        """Plot new data... axis may be specified if 1d...
        """
        if self.active==0:
            self.active=1
            self.win.show()
        #if type(data)==numpy.ndarray:
        #    data=Numeric.array(data)
        if type(data)!=type(None):
            if copy:
                self.data=data.copy().astype("d")
            else:
                if data.dtype.char=="d":
                    self.data=data
                else:
                    self.data=data.astype("d")
        else:#data==None?
            pass
            #if type(self.data)==numpy.ndarray:
            #    self.data=Numeric.array(self.data)
        self.update=1
        #print "plot"
        #print type(axis)
        #if type(axis)!=type(None):
        #    print axis.shape
        gobject.idle_add(self.queuePlot,axis)
##         self.ax.clear()
##         if len(self.data.shape)==1 or self.dims==1:
##             #1D
##             if len(self.data.shape)==1:
##                 freeze,logscale,data,scale=self.mytoolbar.prepare(self.data,dim=1)
##                 if freeze==0:
##                     if type(axis)==type(None) or axis.shape[0]!=data.shape[0]:
##                         axis=Numeric.arange(data.shape[0])+1
##                     if logscale:
##                         try:
##                             axis=Numeric.log10(axis)
##                         except:
##                             print "Cannot take log"
##                     #self.fig.axis([axis[0],axis[-1],scale[0],scale[1]])
##                     self.ax.plot(axis,data)
##             else:#use first row of data for the x axis...
##                 #axis=data[0]
##                 freeze,logscale,data,scale=self.mytoolbar.prepare(self.data,dim=1)
##                 if freeze==0:
##                     if type(axis)==type(None) or axis.shape[0]!=data.shape[-1]:
##                         axis=Numeric.arange(data.shape[0])+1
##                     if logscale:
##                         try:
##                             axis=Numeric.log10(axis)
##                         except:
##                             print "Cannot take log"
##                     #self.fig.axis([axis[0],axis[-1],scale[0],scale[1]])
##                     for i in range(data.shape[0]):
##                         self.ax.plot(axis,data[i])

##         else:#2D
##             if len(self.data.shape)!=2:#force to 2d
##                 self.data=Numeric.reshape(self.data,(reduce(lambda x,y:x*y,self.data.shape[:-1]),self.data.shape[-1]))
##             freeze,logscale,data,scale=self.mytoolbar.prepare(self.data)
##             if freeze==0:
##                 self.ax.imshow(data,interpolation=self.interpolation,cmap=self.cmap,vmin=scale[0],vmax=scale[1])
##         if freeze==0:
##             try:
##                 self.ax.draw()
##             except:
##                 pass
##             #self.ax.update()
##             self.canvas.draw()
##             #self.canvas.queue_draw()
        return True
def randomisePlot(w,p=None):
    d=numpy.random.random((20,20)).astype("f")
    p.dims=3-p.dims
    if p!=None:
        p.plot(d)
    else:
        print "No plot widget"
    return True

def randomise(w,data=None):
    print "Randomising"
    X[:]=numpy.random.random((20,20)).astype("f")
    #ax = fig.add_subplot(111)
    #print dir(ax),ax.figBottom
    #dr=dir(ax)
    #for d in dr:
    #    print d,"       ",getattr(ax,d)
    ax.clear()
    ax.imshow(X,interpolation="nearest",cmap=colour.gray)
    
    #ax.redraw_in_frame()
    ax.draw()
    #event = gtk.gdk.Event(gtk.gdk.EXPOSE)
    #print dir(event),dir(event.area),type(event.area),type(event.area.x),event.area.x,event.area.height
    #event.time = -1 # assign current time
    #canvas.emit("expose_event",event)
    #self.window.draw_drawable (self.style.fg_gc[self.state],self._pixmap, 0, 0, 0, 0, w, h)
    canvas.queue_draw()
    return True
def buttonPress(w,e,data=None):
    if e.button==3:
        print "Right clicked..."
        print dir(toolbar)
        print dir(toolbar.get_parent())
        print dir(toolbar.window)
        if gtk.Object.flags(toolbar)|gtk.VISIBLE:
            print "visible"
        print toolbar.get_child_visible()
        toolbar.hide()
    if e.button==2:
        print "middle clicked"
        ax.clear()
        ax.plot(arange(X.shape[1]),X[0])
        ax.draw()
        canvas.queue_draw()
    return True

if __name__=="__main__":
    simple=0
    if simple:
        ctrlwin=gtk.Window()
        ctrlwin.connect("destroy",lambda x: gtk.main_quit())
        ctrlwin.set_default_size(400,300)
        ctrlwin.set_title("control window")
        button = gtk.Button("Randomise")
        button.connect("clicked", randomise, None)
        ctrlwin.add(button)
        ctrlwin.show_all()



        win = gtk.Window()
        win.connect("destroy", lambda x: gtk.main_quit())
        win.set_default_size(400,300)
        win.set_title("Embedding in GTK")

        vbox = gtk.VBox()
        win.add(vbox)

        vbox.connect("button_press_event",buttonPress)

        fig = Figure(figsize=(5,4), dpi=50)


        #figorig=pylab.figimage(X)
        #fig=figorig.figure
        #print dir(fig)
        #fig2=FigureImage(X)
        #fig2=fig.figimage(X)
        #ax = fig.figure.add_subplot(111)
        ax = fig.add_subplot(111)
        t = arange(0.0,3.0,0.01)
        s = sin(2*pi*t)
        #print dir(ax)
        #ax.plot(t,s)
        ax.imshow(X,interpolation="nearest")
        print type(fig),dir(ax),dir(fig)

        canvas = FigureCanvas(fig)  # a gtk.DrawingArea
        vbox.pack_start(canvas)
        toolbar = NavigationToolbar(canvas, win)
        vbox.pack_start(toolbar, False, False)


        win.show_all()
        gtk.main()
    else:
        p=plot()
        ctrlwin=gtk.Window()
        ctrlwin.connect("destroy",lambda x: gtk.main_quit())
        ctrlwin.set_default_size(40,30)
        ctrlwin.set_title("control window")
        button = gtk.Button("Randomise")
        button.connect("clicked", randomisePlot, p)
        ctrlwin.add(button)
        ctrlwin.show_all()
        
        gtk.main()
