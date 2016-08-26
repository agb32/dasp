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
import Numeric,types,base.aobase
class testmod(base.aobase.aobase):
    """A class that can be used to test the framework.
    @cvar parent: Parent object
    @type parent: Instance
    @cvar msg: A message to be printed when generateNext is called
    @type msg: String
    @cvar data: Data returned
    @type data: Numeric array.
    """
    def __init__(self,parent,msg,data=Numeric.zeros((10,10),"i")):
        """Initialise the test object
        @param msg: Message to be printed.
        @type msg: String
        @param parent: Parent object
        @type parent: None or Instance
        @param data: Data to be returned
        @type data: Numeric array or None
        """
        base.aobase.aobase.__init__(self,parent,None)
        self.msg=msg
        self.parent=parent
        self.data=data
    def needData(self,msg=None):
        """Do we need data?
        @param msg: A fwdMsg
        @type msg: None or Instance
        """
        if self.parent!=None:
            try:
                print "%s: calling parent (%s)"%(self.msg,self.parent.objID)
            except:
                print "%s: calling parent"%self.msg
            return True
        else:
            print "%s: has no parent to call"%self.msg
            return False
    def generateNext(self,msg=None):
        """Generate the next data
        @param msg: A fwdMsg
        @type msg: None or Instance
        """
        print "%s: In generate next"%self.msg
        return self.data
