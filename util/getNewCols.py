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
import math
class getNewCols:
    """A class to determine how many columns (or rows) should be added to the infinite phase screen in a given iteration"""
    def __init__(self,stepSize):
        self.stepSize=stepSize#the step size (floating point)
        self.extraCols=0#current array size above standard.
        self.interpPosition=0.
        self.nadd=0
    def next(self):
        """Returns as a tuple, the number of colums to remove, the number of columns to add, and the interpolation position to use."""
        newpos=self.interpPosition+self.stepSize
        nremove=math.floor(newpos)
        self.nadd=int(math.ceil(newpos)-self.extraCols)
        if newpos%1==0:
            self.nadd+=1
        oldExtraCols=self.extraCols
        self.extraCols+=self.nadd-nremove
        self.interpPosition=newpos-nremove
        return nremove,self.nadd,self.interpPosition


    def nextOld(self):
        """Returns as a tuple, the number of colums to remove, the number of columns to add, and the interpolation position to use.
        This was wrong - try with stepSize=0.75 to see why...
        """
        newpos=self.interpPosition+self.stepSize
        nremove=math.floor(newpos)
        self.nadd=int(math.ceil(newpos)-self.extraCols)
        oldExtraCols=self.extraCols
        self.extraCols+=self.nadd-nremove
        self.interpPosition=newpos-nremove
        return nremove,self.nadd,self.interpPosition
    


if __name__=="__main__":
    import sys
    if len(sys.argv)!=2:
        print "Usage: %s stepsize(floating point)"%sys.argv[0]
        sys.exit(0)
    stepSize=float(eval(sys.argv[1]))
    g=getNewCols(stepSize)
    print "Phasescreen info for a step size of %g"%stepSize
    print "Iteration, Number to remove, number to add, interpolation position"
    for i in range(10):
        print i,g.next()
