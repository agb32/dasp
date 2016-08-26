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
import pygtk,sys
if not sys.modules.has_key("gtk"):
    pygtk.require("2.0")
import gtk, gobject

class myDialog:
    def __init__(self,msg="OK or cancel?",title="OK or cancel",buttons=(gtk.STOCK_OK,gtk.RESPONSE_ACCEPT,gtk.STOCK_CANCEL,gtk.RESPONSE_REJECT),parent=None):
        self.buttons=buttons
        self.d=gtk.Dialog(title,parent,gtk.DIALOG_MODAL,buttons)
        self.d.set_position(gtk.WIN_POS_MOUSE)
        l=gtk.Label(msg)
        l.set_justify(gtk.JUSTIFY_CENTER)
        self.d.vbox.pack_start(l, True, True, 0)
        l.show()
        self.d.connect("close",self.close)
        self.d.connect("response",self.response)
        self.resp="cancel"
        self.d.show()
        gtk.main()
    def close(self,d,arg1=None):
        print "Closing (cancelling dialog)"
        self.d.destroy()
        gtk.main_quit()
    def response(self,d,id,arg1=None):
        for i in range(1,len(self.buttons),2):
            if id==self.buttons[i]:
                self.resp=self.buttons[i-1]
                break
##         if id==gtk.RESPONSE_ACCEPT:
##             self.resp="ok"
##         else:
##             self.resp="cancel"
        if self.resp==gtk.STOCK_OK:
            self.resp="ok"
        elif self.resp==gtk.STOCK_CANCEL:
            self.resp="cancel"
        self.d.destroy()
        gtk.main_quit()
