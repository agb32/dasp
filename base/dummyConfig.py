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
class config:
    """A dummy config, that can be used for testing and debugging"""
    def __init__(self,values):
        """values is a dictionary of values, with keys being the name of the variable to be searched for"""
        self.val=values#a dictionary of values
    def getVal(self,name,default="nodefault"):
        if self.val.has_key(name):
            val=self.val[name]
        elif default!="nodefault":
            val=default
        else:
            raise Exception("Variable %s not found in dummy config"%name)
        return val
    def setSearchOrder(self,val1=None,val2=None):
        pass
