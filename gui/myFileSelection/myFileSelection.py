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
import sys
import gobject
if not sys.modules.has_key("gtk"):
    pygtk.require("2.0")
import gtk, gobject
class myFileSelection:
    """Class to get a filename from user"""
    def __init__(self,title="Open an XML parameter file",file=None,complete="*.xml",parent=None):
        """Initialise file selection dialog"""
        self.fname=None
        self.fs=gtk.FileSelection(title)
        if file!=None:
            self.fs.set_filename(file)
        self.fs.connect("destroy", self.fsdestroy)
        self.fs.ok_button.connect("clicked", self.file_ok_sel)
        self.fs.cancel_button.connect("clicked",self.cancel)
        if file!=None and len(file)>0:
            self.fs.complete(file)
        elif len(complete)>0:
            self.fs.complete(complete)
        self.fs.set_modal(1)
        if parent!=None:
            self.fs.set_transient_for(parent)
        self.fs.show()
        gobject.idle_add(self.move)
        if file!=None and len(file)>0:
            self.fs.complete(file)
        elif len(complete)>0:
            self.fs.complete(complete)
        gtk.main()


    def move(self,a=None):
	x,y=self.fs.get_position()
        m=0
        if x<0:
            x=0
            m=1
        if y<0:
            y=0
            m=1
        if m:
            self.fs.move(x,y)
        return False
        
    def file_ok_sel(self, w):
        """ok pressed"""
        self.fname=self.fs.get_filename()
        #print self.fname
        self.fs.destroy()
        gtk.main_quit()
    def cancel(self,w):
        """cancel pressed"""
        self.fs.destroy()
        gtk.main_quit()
    def fsdestroy(self,w):
        """Quit"""
        gtk.main_quit()
