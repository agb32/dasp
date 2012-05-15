import pygtk,sys
if not sys.modules.has_key("gtk"):
    pygtk.require("2.0")
import gtk, gobject
import gtk.glade as glade
import string

#import base.readConfig,os,types,sys

class selectionbox:
    """Class to implement a simple selection box"""
    def __init__(self,strlist,txt=None,title=None,multiple=0,parent=None):
        self.choice=None
        try:
            gladefile=__file__.split("/")
            gladefile[-1]="selectionbox.glade"
            gladefile=string.join(gladefile,"/")
        except:
            gladefile="selectionbox.glade"
        self.gladetree=glade.XML(gladefile)
        self.sigdict={
                      "on_window1_delete_event":self.on_quit,
                      "on_window1_destroy_event":self.on_quit,
                      "on_buttonCancel_clicked":self.on_quit,
                      "on_buttonOK_clicked":self.on_ok
                      }
        self.gladetree.signal_autoconnect(self.sigdict)
        self.strlist=strlist
        self.multiple=multiple
        if title==None:
            self.title="Selection box"
        else:
            self.title=title
        self.gladetree.get_widget("window1").set_title(self.title)
        if parent!=None:
            self.gladetree.get_widget("window1").set_transient_for(parent)
        if txt!=None:
            self.gladetree.get_widget("labelMsg").set_text(txt)
        self.treeView=self.gladetree.get_widget("treeviewModuleList")
        self.treeStore=gtk.TreeStore(str)
        if self.multiple:
            self.treeView.get_selection().set_mode(gtk.SELECTION_MULTIPLE)
        self.updateList()
        gtk.main()
        
    def updateList(self):
        """Update the selection list"""
        self.treeStore.clear()
        for item in self.strlist:
            self.treeStore.append(None,(item,))
        self.treeView.set_model(self.treeStore)
        for col in self.treeView.get_columns():
            self.treeView.remove_column(col)
        rend=gtk.CellRendererText()
        #rend.set_property("editable",False)
        treecol=gtk.TreeViewColumn(self.title,rend,text=0)
        self.treeView.append_column(treecol)
    def on_quit(self,w=None,arg2=None):
        """quit"""
        self.gladetree.get_widget("window1").destroy()
        gtk.main_quit()
    def on_ok(self,w=None,arg2=None):
        """Selected an option"""
        pathlist=[]
        self.treeView.get_selection().selected_foreach(lambda a,b,c,d: d.append(b), pathlist)
        self.choice=None
        if self.multiple:
            self.choice=[]
            for p in pathlist:
                self.choice.append(p[0])
        else:
            if len(pathlist)>0:
                self.choice=pathlist[0][0]
        self.on_quit()
if __name__=="__main__":
    l=["a","b","c","d"]
    s=selectionbox(l,txt="Choose a letter",title="Alphabet")
    gtk.main()
    print "I got %s which is letter:"%(str(s.choice))
    if s.choice!=None:
        print l[s.choice]
