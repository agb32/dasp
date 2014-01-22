# science/temporalFilter.py
#
# aosim filter operating in the time domain
#
# Author:   Nazim Bharmal
# Date:     Jan/2014
# Version:  not for production
# Platform: Darwin
# Python:   2.6.7
# Numpy:    1.5.1

# This code based on science/template.py

import template
numpy=template.numpy

class aosimTemporalFilter(template.aosimNessuno):
    """aosim temporal filter based on a 1D filter
    [todo, epydoc compatible documentation]
    """
    def _doInitializationCode(self):
           # \/ load configuration variables
        self.active=self.config.getVal(
              "active",
              default=1, raiseerror=0 ) # whether active or not
        self.filterArr=self.config.getVal(
              "filterArr",
              default=None, raiseerror=0 ) # define the filter explicitly...
        self.filterArrAutoNorm=bool(self.config.getVal(
              "filterArrAutoNormalize",
              default=1, raiseerror=0 )) # ...and choose if to normalize or not
        self.filterArrAutoNormTol=float(self.config.getVal(
              "filterArrAutoNormalizeTolerance",
              default=1e-2, raiseerror=0 )) # ...and to within what precision
        filterLen=int(self.config.getVal(
              "filterLen",
              default=0, raiseerror=0 )) # ...or make a filter based on length

        if self.filterArr!=None:
            if filterLen>0:
                raise ValueError("Choose either filterArr *or* filterLen")
            try:
                self.filterArr=numpy.array(self.filterArr,astype=numpy.float32)
            except:
                raise ValueError(
                    "Could not cast the filterArr variable as an array")
            if len(self.filterArr.shape)!=1 or self.filterArr.shape[0]<2:
                raise ValueError(
                    "Filter arr is required to be a 1D array of length >2")
            if abs(self.filterArr.sum()-1)>=self.filterArrAutoNormTol:
                if self.filterArrAutoNorm:
                    self.filterArr/=self.filterArr.sum()
                else:
                    print("WARNING: the filter is *not* normalized to within "+
                            "+/-{0:f}".format(filterArrAutoNormTol))
            print("INFORMATION(**temporalFilter**) using spec'd array, len="
                  +str(len(self.filterArr)))
        else:
            if filterLen<2:
                raise ValueError("FilterLen should be >=2")
            self.filterArr=numpy.empty(filterLen)
            self.filterArr[:]=filterLen**-1.0 # normalize
            print("INFORMATION(**temporalFilter**) making box-car, n="
                  +str(filterLen))

        # here we do not know how big the input data is, so do not yet
        # define the memory that stores previous data
        self.previousData=None
   
    def _doSomeProcessing(self,ip):
        # apply the relevant code
        if self.previousData==None:
            self.previousData=numpy.zeros(
                   [ len(self.filterArr) ]+
                   list(numpy.array(ip).shape),
                ip.dtype ) 
            self.filterArr=self.filterArr.reshape([-1,1])
        self.previousData[1:]=self.previousData[:-1] # shift data
        self.previousData[0]=ip
        return (self.previousData*self.filterArr).sum(axis=0)
