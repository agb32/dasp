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
import Numeric
import cmod.utils, mmap, os, stat,unittest, os.path
"""Code to create a memory mapped array.  This is probably better than my c code because it may allow the array to be unmapped if required.
I think this still isn't 100% fool proof - if you del the array and then close the mmap file and del the mmap object, the ref count of the mmap object will still be 1 because arrayfrombuffer increments it but doesn't dec it.  So there will always be a small memory leak...
But I could be wrong!
Actually, I think this doesn't work - don't bother using it.
My c code now allows deallocation etc anyway.
"""


class mmapArray:
    def __init__(self,file, typecode='l', shape=None):
        """mmap a file as an array, either given a file object or a filename.

        Returns the array, the mmap object (so you can call flush
        explicitly to ensure your changes have gotten copied out to disk),
        and the open file.

        Note, this won't work for files > 2GB.

        """
        self.typecode=typecode
        self.file=file
        
        if shape==None:#shape is equal to file size
            if type(file)==type(""):#is the filename
                if not os.path.exists(file):#doesn't exist...
                    raise Exception("Cannot create mmap file of unknown size")
                shape=(os.stat(file).st_size/self.sizeof(typecode),)
            else:#is the file object...
                pos=f.tell()
                f.seek(0,2)#move to end of file...
                shape=(f.tell()/self.sizeof(typecode),)
                f.seek(pos)
            
        if type(file)==type(""):
            self.f=open(file,"a+")
        else:
            self.f=file
        self.shape=shape
        self.size=reduce(lambda x,y:x*y,shape)*self.sizeof(typecode)
        self.f.truncate(self.size)
        #if not hasattr(file, 'fileno'): file = open(file, 'rb+')
        fn = self.f.fileno()
        self.mem = mmap.mmap(fn, self.size)
        self.arr = cmod.utils.arrayfrombuffer(self.mem, typecode)
        self.arr.shape = shape
        #return arr, mem, file#use m.close() when finished - but dont try to access a after that...

    def sizeof(self,typecode):
        a=Numeric.zeros((1,),typecode)
        s=a.itemsize()
        return s
    def __del__(self):
        """Destroy the mmap object etc"""
        #first del arr.
        del(self.arr)
        #then unmap the memory.
        self.mem.close()
        

class maparraytest(unittest.TestCase):
    def setUp(self):
        self.x = Numeric.arange(10) / 2.0
        open('tmp.foo', 'wb').write(self.x.tostring())
    def tearDown(self):
        os.unlink('tmp.foo')
    def testx(self):
        self.failUnless(Numeric.alltrue(self.x ==
                                        Numeric.array([0, 0.5, 1, 1.5, 2,
                                                       2.5, 3, 3.5, 4, 4.5])))
    def testread(self):
        y, _, _ = maparray('tmp.foo', 'd')
        self.failUnless(Numeric.alltrue(self.x == y))
    def testwrite(self):
        y, mapobj, _ = maparray('tmp.foo', 'd')
        y[1] = 37
        mapobj.flush()
        z, _, _ = maparray('tmp.foo', 'd')
        self.failUnless(Numeric.alltrue(z ==
                                        Numeric.array([0, 37, 1, 1.5, 2,
                                                       2.5, 3, 3.5, 4, 4.5])))
        self.failUnless(Numeric.alltrue(z == y))
    def testimplicitflush(self):
        y, _, _ = maparray('tmp.foo', 'd')
        y[1] = 37
        del y, _  # trigger implicit flush as refcounts go to 0
        z, _, _ = maparray('tmp.foo', 'd')
        self.failUnless(Numeric.alltrue(z ==
                                        Numeric.array([0, 37, 1, 1.5, 2,
                                                       2.5, 3, 3.5, 4, 4.5])))
    def testnewlines(self):
        # If we forgot to open the file in binary mode, some bytes would
        # get mangled on some systems.
        x = Numeric.arange(258).astype('b')  # bytes
        x[-2:,] = (13, 10)  # CRLF
        open('tmp.foo', 'wb').write(x.tostring())
        y, _, _ = maparray('tmp.foo', 'b')
        self.failUnless(Numeric.alltrue(x == y))
    def testshape(self):
        xx = [[1, 3], [5, 7]]
        open('tmp.foo', 'wb').write(Numeric.array(xx).tostring())
        self.assertEqual(maparray('tmp.foo', shape=(2, 2))[0].tolist(), xx)
    # There are a bunch of features that are sort of hard to test, because
    # they're performance-oriented and require big data files or because
    # they involve transcending system limitations that may not be present.
    # For example:
    # - mapping an array should return instantly, even for large files
    # - reading the contents of a mapped array shouldn't take much more time
    #   than reading in the corresponding file with ordinary read() calls
    # - it should be possible to have arrays much larger than RAM
    #   (this is obviously hard to test if you have more than a moby of RAM)
    #   and operate on them reasonably efficiently.
    # - reading or writing a small part of the array should be fast
