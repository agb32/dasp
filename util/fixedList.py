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
