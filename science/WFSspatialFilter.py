# science/WFSspatialFilter.py
#
# aosim spatial filter implementation
#
# Author:   Nazim Ali Bharmal
# Date:     June/2013
# Version:  not for production
# Platform: Darwin
# Python:   2.6.7
# Numpy:    1.5.1

# Lines marked with '#?' are those that may be important but don't appear
# to be necessary.
#
# This code based on the example in 'dummiesGuide.pdf', 'wfscent.py',
# 'science.py'.

# Assumes a 'phaseonly' value for <atmosPhaseType> and then assumes the next
# module (WFS) will accept a 'phaseamp' value for the same variable i.e.
# it changes the data format, so the parameters file will require changes
# in the variable value from default for the modules that follow this one.

import base.aobase
import collections
import numpy
import time

class spatialFilter(object):
   '''Create a persistant object which can apply a spatial filter.
      The basic idea is that upon initialisation, the mask can be constructed
      using specified filter frequency cutoff and shape. It then need not be
      updated.
      Within application, take the input phase and with knowledge of the
      aperture, construct the complex amplitude.  Then apply a DFT to the
      complex amplitude, apply the mask, and apply a IDFT.  Split into phase
      and amplitude and return.
      NOTE: `size' is defined by the extent of the mask relative to sampling
      and this sampling can be fractional; the number of pixels within the mask
      may exceed the number expected, but the frequency cut-off will
      nonetheless be correct and this is the effect of interest in this
      code.'''

   filterType=None
   filterMasks=None
   filterSize=None # circle: diameter, square: width
   apertureMask=None
   haveUnwrapper=False
   n=None
   _cds = lambda self,n : numpy.fft.fftshift((numpy.arange(n)-n//2)*n**-1.0)
    
   def __init__(self,filterType,filterSize,apertureMask,fpScaling=[0],
         paddingExpansion=1, unwrap=True):
      '''
      [document using epydoc compatible format]
      @str filterType The type of the filter, currently either
      a fully transparent circle or square in a fully opaque holder.
      @float filterSize fraction of Nyquist for filter configuration
         => <=0.5.
      @numpy.array apertureMask obv.
      (opt) @iter(int) fpScaling how many pixels to adjust the
         aperture by in order to simulate polychromatism by changing the scale
         of the focal plane, default is to have a one shift of zero (i.e.
         monochromatic).
      (opt) @int paddingExpansion The upper limit in powers of 2 which the 
         padding will increase the apertureMask by, default is 1 (canonical
         no.)
      (opt) @bool unwrap Whether to use a defined unwrapping function, default
         is yes.
      '''
      self.applyFunction=self._applyPython
      self.knownFilterTypes=(
            ('circle',     type(0.0),  self._makeFilterCircle),
            ('square',     type(0.0),  self._makeFilterSquare),
            ('rectangle',  type((0,)), self._makeFilterRectangle),
         )
      if type(filterType)!=type(""):
         raise TypeError("spatialFilter/Filter type must be a string")
      self.filterType=filterType # string
      if filterSize>0.5: raise ValueError("spatialFilter/0<filterSize<=0.5")
      self.filterSize=filterSize # [spatial frequency] 
         # \/ cast, util.tel.Pupil -> numpy.array
      self.apertureMask=numpy.array(apertureMask)
         # \/ relative power (n.b. amplitude==1 hence no square)
      self.aperturePower=(self.apertureMask!=0).sum() 
      self.n=self.apertureMask.shape
         # \/ find the nearest power of 2, rounded up, plus paddingExpansion to
         #    derive how big the padded apertureMask is
      self.n_=numpy.power(2,(paddingExpansion+
          numpy.ceil(numpy.log(numpy.array(self.n))/numpy.log(2))
            ).astype(numpy.int16))
      if not isinstance(fpScaling, collections.Iterable):
         raise TypeError(
               "spatialFilteR: fpScaling should be iterable integers")
      self.fpScaling=fpScaling
      self.foundZeroShift=False
      for thisFPS in self.fpScaling:
         if thisFPS%2!=0 or type(thisFPS)!=type(0):
            raise TypeError("spatialFilter/Wavefront shift ({0:d}) must be "+
                  "even integers".format(thisFPS))
         for i in (0,1):
               # *** assumption made that apertureMask is sized /exactly/ to
               #  hold the aperture and so it is the minimum acceptable size.
            if self.n_[i]-thisFPS<self.n[i]:
               raise ValueError("spatialFilter/Wavefront shift ({0:d}) is "+
                     "too large (by {1:d}) for the aperture size".format(
                     thisFPS,self.n[i]-(self.n_[i]-thisFPS)))
         if thisFPS==0: self.foundZeroShift=True
      if not self.foundZeroShift: 
         print("spatialFilter/(WARNING) "+
               "No zero shift found, it is sensible to include this value")
         # \/ relative (to the zero change focal plane scaling) wavelength
         #  shifts. Typically this zero change is 'wfs_lam' in the params
         #  file.
      self.relativeWaveShifts=[
            1-thisFPS*self.n_[0]**-1.0 for thisFPS in self.fpScaling ]
      for thisFilterType in self.knownFilterTypes:
         if (filterType==thisFilterType[0]
               and type(filterSize)==thisFilterType[1]):
            self._makePolychromaticFilters( thisFilterType[2], filterSize )
      if self.filterMask==None:
         raise ValueError("Could not make a filter mask, currently know about:"+
               string.join( 
                  [ string.join([str(y) for y in x],",")
                     for x in self.knownFilterTypes ])
                     )
         # \/ ATTN USER: Set this to false if you do not have an unwrapping
         #   routine.
      self.haveUnwrapper=True
      if unwrap:
         if not self.haveUnwrapper:
            raise RuntimeError("Cannot ask for unwrapping without supplying"+
               " an unwrapping routine")
         # ATTN USER: Insert your own phase-unwrapping code here.
         self.unwrapperParametersReturn={    
            'params':('phase','mask'),       # Darwin 9.8, 2012/10
            'return':('phase',)              # Darwin 9.8, 2012/10
           }
         import punwrap                      # Darwin 9.8, 2012/10
         self.unwrapper=punwrap.unwrap2D     # Darwin 9.8, 2012/10
      else:
         self.unwrapper=None

   def apply(self,phase,passThrough):
      return self.applyFunction(phase,passThrough)

   def _applyC(self,phase,passThrough):
      '''Apply the algorithm with a partial C-implementation'''
      raise RuntimeError("spatialFilter/_applyC is not implemented")

   def _applyPython(self,phase,passThrough):
      '''Apply the algorithm with a Python/numpy implementation'''
         # nb numpy.fft which should be sufficiently quick
      metapupil=numpy.empty(
            [len(self.fpScaling)]+list(self.n), numpy.complex128)
      phaseampPupil=numpy.empty( [2]+list(self.n), numpy.float32)
         # \/ the following need to be as small as the smallest focal plane
      metafocalPlane=numpy.empty(list(self.n_-max(self.fpScaling)), numpy.float32)
      zeroShiftfocalPlane=numpy.zeros(list(self.n_-max(self.fpScaling)),
            numpy.float32)
      if not self.foundZeroShift:
         # pattern with a cross
         # TODO
         pass
      for fpSI,thisfpS in enumerate(self.fpScaling):
         pupil=self.apertureMask*numpy.exp(
               1.0j*phase*self.relativeWaveShifts[fpSI])
         newpupil=numpy.zeros(self.n_-thisfpS,numpy.complex128)
         newpupil[:self.n[0],:self.n[1]]=pupil         
# (old)          newpupil[:self.n[0]/2,:self.n[1]/2]=pupil[-self.n[0]/2:,-self.n[1]/2:]
# (old)          newpupil[-self.n[0]/2:,:self.n[1]/2]=pupil[:self.n[0]/2,-self.n[1]/2:]
# (old)          newpupil[:self.n[0]/2,-self.n[1]/2:]=pupil[-self.n[0]/2:,:self.n[1]/2]
# (old)          newpupil[-self.n[0]/2:,-self.n[1]/2:]=pupil[:self.n[0]/2,:self.n[1]/2]
         thisfocalPlane=numpy.fft.fft2(newpupil) # always make a focal plane
         if not passThrough: thisfocalPlane*=self.filterMask[fpSI] 
         metapupil[fpSI]=numpy.fft.ifft2(thisfocalPlane)[
                  :self.n[0],:self.n[1] ]
# (old)          metapupil[fpSI]=numpy.fft.fftshift( numpy.fft.ifft2(thisfocalPlane)
# (old)                )[ self.n_[0]/2-self.n[0]/2:self.n_[0]/2+self.n[0]/2,
# (old)                   self.n_[1]/2-self.n[1]/2:self.n_[1]/2+self.n[1]/2 ]
            # \/ record the focal plane
         thisfocalPlaneIntensity=abs(numpy.fft.fftshift(thisfocalPlane)[
                  (max(self.fpScaling)-thisfpS)/2:
                     self.n_[0]-(thisfpS+max(self.fpScaling))/2,
                  (max(self.fpScaling)-thisfpS)/2:
                     self.n_[1]-(thisfpS+max(self.fpScaling))/2,
                   ])**2.0
         if thisfpS==0:
            zeroShiftfocalPlane=thisfocalPlaneIntensity
         metafocalPlane+=thisfocalPlaneIntensity

      # \/ process up the polychromatic pupil to form a monochromatic
      #  approximation. The key fact here is that the ultimate interest
      #  is in the focal plane intensity which is strictly modelled by,
      #  FPI[i,j]=\sigma_k{||F(pup[i,j]_k)_k||^2}
      #  where i,j are the sub-aperture indices and k is the wavelength
      #  index.
      #  If loosely the following approximation is made,
      #  FPI[i,j]\simeq{||F( \sigma_k{pup[i,j]_k} )||^2}
      #  then the lack of interference between different wavelength is ignored,
      #  so it is clearly wrong. If these cross-terms produce a /similar/
      #  result to the self-terms (hypothesis) which will not be far wrong in
      #  the corrected case then the correction factor becomes a constant.
      #  Assuming k runs to 1,N then the extra fractional power is N so the
      #  constant in amplitude is (N)^{-1/2}. Is it important to apply a
      #  constant?
      #  Unknown, but for consistency it is done so.  Note that much of the
      #  issue here is due to the lack of chromatic support in the aosim code.
      #  (begins)
      metapupil=metapupil.sum(axis=0)*len(self.fpScaling)**-0.5
      #  (ends)
      phaseampPupil[0]=numpy.arctan2(
            metapupil.imag,metapupil.real ).astype(numpy.float32)
      phaseampPupil[1]=abs(metapupil)
      if self.unwrapper!=None: # presume an unwrapper exists
         # \/ unwrap the phase, it isn't clear if this step is reqd.
         if self.unwrapperParametersReturn['params']==('phase','mask')\
               and self.unwrapperParametersReturn['return']==('phase',):
            phaseampPupil[0]=self.unwrapper(
                  phaseampPupil[0], self.apertureMask )
         else:
            raise RuntimeError("Do not understand the phase unwrapper")
      relativePower=(len(self.fpScaling)*self.aperturePower)**-1.0\
            *(abs(phaseampPupil[1])**2.0).sum()
      metafocalPlane*=len(self.fpScaling)**-1.0
      return phaseampPupil,(metafocalPlane,zeroShiftfocalPlane),relativePower
#"""  <begins> To slice out the middle of the focal plane:
#"""      [
#"""            int(self.n_[0]/2-2*self.filterSize*self.n_[0]):
#"""             int(self.n_[0]/2+2*self.filterSize*self.n_[0]),
#"""            int(self.n_[1]/2-2*self.filterSize*self.n_[1]):
#"""             int(self.n_[1]/2+2*self.filterSize*self.n_[1]) ]
#"""  <ends>

   def _makePolychromaticFilters( self, filterFn, filterSize ):
      '''[private] create all the filters for all the wavelengths
      '''
      self.filterMask=[ filterFn(self.n_-thisFPS, filterSize)
            for thisFPS in self.fpScaling ]
   def _makeFilterCircle(self,l,s):
      '''[private] create a circular mask of radius s
      '''
      return numpy.where( numpy.add.outer(
            self._cds(l[0])**2.0, self._cds(l[1])**2.0 ) <= s**2.0, 1, 0 )
   def _makeFilterRectangle(self,l,s):
      '''[private] create a rectangle mask of width 2s[0] and height of 2s[1]
      '''
      return numpy.multiply.outer(
            numpy.where(abs(self._cds(l[0]))<=s[1],1,0),
            numpy.where(abs(self._cds(l[1]))<=s[0],1,0) )
   def _makeFilterSquare(self,l,s):
      '''[private] create a square mask of dimension 2s X 2s
      '''
      return self._makeFilterRectangle(l,[s,s])


class aosimSpatialFilter(base.aobase.aobase):
    """aosim spatial filter
    [todo, epydoc compatible documentation]
    """
    def __init__(self,parent,config,args={},forGUISetup=0,passThrough=False,debug=None,idstr=None):
        """ [todo, epydoc compatible documentation]
        """
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.control['close_dm']=0 # assume open-loop to start with
        #  General configuration
        self.npup=self.config.getVal("npup")    # Pupil phase array size
        self.pupil=self.config.getVal("pupil")  # Pupil mask
        atmosPhaseType=self.config.getVal("atmosPhaseType",default="phaseonly")
        if atmosPhaseType!="phaseonly":
           raise RuntimeError("ERROR: The spatial filter expects the atmosPhaseType to be 'phaseonly'")
#?  <begins> The following two variables don't seem to be necessary because we
#?    do not have polychromatic support.
#?        self.wfsLam=self.config.getVal("wfslam",default=sourceLam) # WFS wavelength
#?        self.telDiam=self.config.getVal("telDiam") # Pupil physical diameter
#?  <ends>
        #
        # Specific configuration
        #
           # \/ open the filter when the loop is open?
        self.filterPadding=self.config.getVal('filterPadding',default=1)
        if type(self.filterPadding)!=type(0) or self.filterPadding<1:
           raise ValueError("Require that filterPadding is an integer greater "+
              +"than zero (default=1)")
        self.filterAdaptable=self.config.getVal('filterAdaptable',default=0)
        self.filterType=self.config.getVal('filterType',default='square')
        if type(self.filterType)!=type(""):
            raise ValueError("Require that filterType is a string, "+
                  "suggest a value of 'square'")
        if self.pupil.shape[0]!=self.pupil.shape[1]:
           self.filterSize=self.config.getVal('filterSize')
        else:
              # filter out frequencies beyond a sub-aperture
           minimalSampling=self.config.getVal("wfs_nsubx")*self.npup**-1.0
           self.filterSize=float(self.config.getVal(
               'filterSize',default=minimalSampling))
           if self.filterSize<0 or self.filterSize>0.5:
               raise ValueError(
                     "Require that the filter size is bounded, >0, <=0.5")
        self.unwrap=self.config.getVal('filterUnwrapping',default=False)
        self.filterfpScaling=self.config.getVal(
            'filterFocalPlaneScaling',default=[0])
        if not isinstance(self.filterfpScaling, collections.Iterable):
            raise ValueError("Require that filterfpScaling be iterable, "+
                  "should contain integers")
        self.passThrough=passThrough
        self.spatialFilter=spatialFilter(
              self.filterType,
              self.filterSize,
              self.pupil,
              self.filterfpScaling,
              self.filterPadding,
              self.unwrap)
   
    def generateNext(self,msg=None):
        """[todo, epydoc compatible documentation]
        """
        if self.generate==1:
            if self.newDataWaiting:#this is always 1
                if self.parent.dataValid==1:
                    #self.inputData=self.parent.outputData
                    self.inputInvalid=0
                else:
                    self.inputInvalid=1
                    self.dataValid=0
            if self.inputInvalid==0: # there was an input, so we can integrate...
                self.dataValid=1
                # if the loop is open, we probably don't want to filter
                passThrough=self.passThrough
                if (not passThrough and not self.control['close_dm'] and
                      self.filterAdaptable):
                   # if open loop then do not apply filter
                   passThrough=1
                (self.outputData,self.intermediateFocus, self.throughput)=\
                     self.spatialFilter.apply(
                         self.parent.outputData, passThrough)

    def plottable(self,objname="$OBJ"):
        """[todo] write epydoc compatible description"""
        op=""
        op+=('<plot title="sf polychromatic focus" '+
             'cmd="data=%s.intermediateFocus[0]" ret="data" type="pylab" '+
             'when="rpt" palette="jet"/>\n')%(objname)
        op+=('<plot title="sf monochromatic focus" '+
             'cmd="data=%s.intermediateFocus[1]" ret="data" type="pylab" '+
             'when="rpt" palette="jet"/>\n')%(objname)
        op+=('<plot title="sf output amp" '+
             'cmd="data=%s.outputData[1]" ret="data" type="pylab" '+
             'when="rpt" palette="jet"/>\n')%(objname)
        op+=('<plot title="sf output phs" '+
             'cmd="data=%s.outputData[0]" ret="data" type="pylab" '+
             'when="rpt" palette="gray"/>\n')%(objname)
        op+=('<plot title="sf throughput (%%)" '+
             'cmd="data=%s.throughput*100" ret="data" when="rpt" texttype="1" '+
             'wintype="ownwindow" textreplace="1"/>\n')%(objname)
        return(op)

#? 
#?     def getParams(self):
#?         """parameters required for this module"""
#?         return {"myvalue":99.8,"myothervalue":99.7}
#? 
#?     def getInputType(self):
#?         """Returns the input needed by this module.
#?         Should return an instance of dataType, or a list of dataType objects,
#?         or None if no input is required."""
#?         return base.dataType.dataType("My required input",Numeric.ArrayType,(10,10),Numeric.Int32)
#? 
#?     def getOutputType(self):
#?         """Returns the output given by this module."""
#?         return base.dataType.dataType("my output",Numeric.ArrayType,(20,20),Numeric.Float32)
#? 
#?     def getInitArgs(self):
#?         """return a dictionary of parameters that can be passed to init.
#?         Dictionary key is the parameter name, and dictionary value is a tuple
#?         of (description, default value).  If there is no default value, then
#?         the dictionary value is a tuple of (description,).
#?         """
#?         return {"parent":("Predecessor object",),"config":("Configuration object",),
#?                 "args":("Optional arguments",{}),"forGUISetup":("initialisation type",0)}
