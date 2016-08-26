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
# science/temporalFilter.py
#
# aosim filter operating in the time domain
# Has the ability to operate with sequences and invalid inputs
#
# Author:   Nazim Bharmal
# Date:     Nov/2015
# Version:  not for production
# Platform: Linux
# Python:   2.7.3
# Numpy:    1.6.2

# This code based on science/template.py

import template
numpy=template.numpy

class aosimSequentialFilter(template.aosimNessuno):
    """aosim sequential filter with capability to apply a 
    filter to each member of the sequence.
    This version requires dataValid to always be true at each iteration

    filterM has dimensions 
      {N_seq, N_data}
    where
      N_seq is the number of sequences of filters to use,
      N_data is the number of input data the filter acts over
    so as examples,
      filterM = [[0.5,0.5]] : {1,2}
       averages over two input data at each iteration
      filterM = [[1],[0]] : {2,1}
       alternates between passing through (x1) the input data and nulling it 
       (x0)
    
    input has dimensions
     {N_seq,N_ip}
    where
     N_seq is as above,
     N_ip is the number of data in the input

    If input has dimension {N_ip} then it is converted to a sequence of length
    1.
    An inherent memory with dimension {N_data,N_seq,N_ip} is rolled once per
    iteration so that the action of a filter with N_data equal to 1/N_data is
    equivalent to a rolling average.

    The meaning of sequence only has an interpretation for 2D inputs, N_seq!=1.
    For each iteration, the sequence number is increased and the appropriate
    row from filterM and from input is chosen, and a sum over the axis for
    N_data is computed and returned.

    [todo, epydoc compatible documentation]
    """
    def _doInitializationCode(self):
           # \/ load configuration variables
        self.active=self.config.getVal(
              "active",
              default=1, raiseerror=0 ) # whether active or not
        self.filterM=self.config.getVal(
              "filterM",
              default=None, raiseerror=0 ) # define the filter explicitly...
        filterLen=int(self.config.getVal(
              "filterLen",
              default=0, raiseerror=0 )) # ...or make a filter based on length

        if not self.filterM is None:
            assert filterLen==0, ValueError(
                    "Choose either filterM *or* filterLen"
                )
            try:
                self.filterM=numpy.array(self.filterM,dtype=numpy.float32)
            except:
                raise ValueError(
                    "Could not cast filterM to an array")
            assert self.filterM.shape[0]>0, ValueError("filterM.shape[0]>0")
            #
            if len(self.filterM.shape)==1:
                # to support older code, can specify a 1D array which we
                # convert to a 2D matrix to support sequences
                self.filterM.resize( [1,-1] )
            assert len(self.filterM.shape) is 2, ValueError( 
                    "filterM should be 2D"
                )
        else:
            assert filterLen>0, ValueError("filterLen<=0")
            self.filterM = numpy.ones([filterLen,1])

        # here we do not know how big the input data is, so do not yet
        # define the memory that stores previous data, and then
        # we work out how many sequences there are
        self.previousData = None
        self.numSequences = self.filterM.shape[0]
        self.sequenceNumber = 0

    def generateNext(self,msg=None):
        """Called per iteration.
        It has to be able to handle parent.dataValid==False by 
        passing on the value to _doSomeProcessing and then
        accepting the self.dataValid value directly from that function's
        o/p.
        """
        if self.generate is 1:
            if self.newDataWaiting:#this is always 1
                self.inputInvalid = 0 if self.parent.dataValid else 1
                self.nValid += 0 if self.inputInvalid else 1
                if self.active:
                    # \/ apply the function and accept data validity as
                    #    one return value
                    self.outputData, self.dataValid = self._doSomeProcessing(
                            self.parent.outputData, self.inputInvalid
                        )
                else:
                    # \/ if not active then pass data along & choose data
                    #    validity depending on input invalidity
                    self.dataValid = not self.inputInvalid
                    self.outputData = self.parent.outputData

    def _doSomeProcessing(self, ip, inputInvalid):
        #
        # \/ create an appropriate input
        if not inputInvalid:
           thisIp=numpy.array(ip)
           assert len(ip.shape)<3, "Can only accept 1D or 2D input"
           thisIp=ip if len(ip.shape)==2 else ip.reshape([1,-1])
           if thisIp.shape[0]!=self.filterM.shape[0]:
               if self.filterM.shape[0]==1:
                   # if 1D then can extend filterM to accomodate
                   self.filterM.resize( [thisIp.shape[0],-1] )
                   self.filterM[1:] = self.filterM[0]
                   print("WARNING(**sequentialFilter**) "+
                           "filterM has been extended to match the data"
                       )
               else:
                   raise ValueError( "Input.shape[0]!=filter array.shape[0]: "+
                           "{:d} vs. {:d}".format(
                                   ip.shape[0], self.filterM.shape[0]
                               )
                       )
           #
           # \/ store the input
           if self.previousData==None:
               self.previousData=numpy.zeros(
                      [self.filterM.shape[1]] + list(thisIp.shape),
                      thisIp.dtype
                   )
           else:
              # \/ shift...
              self.previousData = numpy.roll( self.previousData, -1, axis=0 )
           self.previousData[-1] = ip  # < ...& replace
        else:
            # the input was invalid, so we just use the previous data as is
            # i.e. the rolling average is not moved but the sequence can be
            pass
        # \/ compute the filtered result
        result = numpy.sum(
                  self.previousData[:, self.sequenceNumber ]*
                  self.filterM[ self.sequenceNumber ].reshape([-1,1]),
                  axis=0
            )
            # wrap the sequence selection if required
        self.sequenceNumber = (self.sequenceNumber+1)%self.numSequences
        # \/ always return the result and dataValid=1
        return (result, 1)

class aosimTemporalFilter(aosimSequentialFilter):
    """aosim temporal filter based on a 1D filter
    [todo, epydoc compatible documentation]
    """
    def _doInitializationCode(self):
        aosimSequentialFilter._doInitializationCode(self)
        #
        # \/ load additional configuration variables
        self.filterMAutoNorm=bool(self.config.getVal(
              "filterMAutoNormalize",
              default=1, raiseerror=0 )) # ...and choose if to normalize or not
        self.filterMAutoNormTol=float(self.config.getVal(
              "filterMAutoNormalizeTolerance",
              default=1e-2, raiseerror=0 )) # ...and to within what precision
        filterLen=int(self.config.getVal(
              "filterLen",
              default=0, raiseerror=0 )) # ...or make a filter based on length
        #
        if not self.filterM is None:
            if abs(self.filterM.sum()-1)>=self.filterMAutoNormTol:
                if self.filterMAutoNorm:
                    self.filterM/=self.filterM.sum()
                else:
                    print("WARNING: the filter is *not* normalized to within "+
                            "+/-{0:f}".format(filterMAutoNormTol))
            print("INFORMATION(**temporalFilter**) using spec'd array, len="
                  +str(len(self.filterM)))
        else:
            self.filterM[:]=filterLen**-1.0 # normalize
            print("INFORMATION(**temporalFilter**) making box-car, n="
                  +str(filterLen))
