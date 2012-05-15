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
