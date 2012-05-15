import pygtk
import sys
if not sys.modules.has_key("gtk"):
    pygtk.require("2.0")
import gtk, gobject
import gtk.glade as glade
import string
from gui.myFileSelection.myFileSelection import myFileSelection

class textbox:
    """Class to display a box with text in it"""
    def __init__(self,title=None,closeFunc=None):
        """initialise"""
        gladefile=__file__.split("/")
        gladefile[-1]="textbox.glade"
        gladefile=string.join(gladefile,"/")
        self.gladetree=glade.XML(gladefile)
        self.sigdict={"on_buttonSave_clicked":self.save,
                      "on_buttonClose_clicked":self.close,
                      "on_window1_delete_event":self.close
                      }
        self.gladetree.signal_autoconnect(self.sigdict)
        self.textview=self.gladetree.get_widget("textview1")
        self.title=title
        self.destroyed=0
        self.closeFunc=closeFunc
        if title!=None:
            self.gladetree.get_widget("window1").set_title(title)
    def save(self,w,arg2=None):
        """save the text"""
        fn=myFileSelection("Choose a text file name").fname
        if fn!=None:
            f=open(fn,"w")
            tb=self.textview.get_buffer()
            f.write(tb.get_text(tb.get_start_iter(),tb.get_end_iter()))
            f.close()
    def close(self,w,arg2=None):
        """close the window"""
        self.gladetree.get_widget("window1").destroy()
        self.destroyed=1
        if self.closeFunc!=None:
            self.closeFunc()
    def addText(self,txt,replace=0):
        """add (or replace) text in the window"""
        tb=self.textview.get_buffer()
        if replace:
            tb.set_text(txt)
        else:
            tb.insert(tb.get_end_iter(),txt)
        self.textview.scroll_to_mark(tb.get_insert(), 0)
        
if __name__=="__main__":
    t=textbox()
    gtk.main()
    
