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
#import Numeric
import numpy
class fixedList:
    """Class to implement a list in an array, ie having a fixed length,
    determined on initialisation.  If the array gets full,
    it wraps around and starts overwriting.
    @cvar arr: Array holding the data
    @type arr: Numeric array
    @cvar listlen: Max length of list
    @type listlen: Int
    @cvar wrapped: Whether has overflowed yet
    @type wrapped: Int
    @cvar cnt: Current position in list
    @type cnt: Int
    """
    def __init__(self,listlen,typecode):
        """Initialise.
        @param listlen: Max length of array
        @type listlen: Int
        @param typecode: Type of array, e.g. numpy.int32 etc
        @type ptypecode: Character
        """
        self.arr=numpy.zeros((listlen,),typecode)
        self.listlen=listlen
        self.wrapped=0
        self.cnt=0
    def append(self,data):
        """Append data to the last position in the list.
        @param data: The data to append
        @type data: Defined by typecode
        """
        self.arr[self.cnt]=data
        self.cnt+=1
        if self.cnt==self.listlen:
            self.cnt=0
            self.wrapped=1
    def __repr__(self):
        return str(self.aslist())
    def aslist(self):
        """Return the object as a list."""
        if self.wrapped==0:
            lst=list(self.arr[:self.cnt])
        else:
            lst=list(self.arr[self.cnt:,])+list(self.arr[:self.cnt])
        return lst
