# science/hwr.py
#
# aosim HWR reconstructor
#
# Author:   Nazim Ali Bharmal
# Date:     Sept/2013
# Version:  not for production
# Platform: Darwin
# Python:   2.6.7
# Numpy:    1.5.1

# Lines marked with '#??' are those that may be important but don't appear
# to be necessary.
#
# This code based on 'tomoRecon.py' by inheritance and replication of members
# where neccessary.


import abbot
import abbot.hwr
import base.aobase
import numpy
import os.path
import Scientific.MPI
import scipy.sparse
import scipy.linalg
from scipy.sparse.linalg import cg as spcg
import sys
import time
import tomoRecon
import util.FITS


class recon(tomoRecon.recon):
   """Clarified reconstructor for HWR.
   """
    
   def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
      global util
      newReconType="hwr"
      tomoRecon.recon.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
      
         # \/ Get HWR-specific configuration values
      self.hwrMaxLen=self.config.getVal( "hwrMaxLen",
            default=self.config.getVal("wfs_nsubx"),raiseerror=0) 
      self.hwrBoundary=self.config.getVal( "hwrBoundary",
            default=None,raiseerror=0) 
      if type(self.hwrBoundary)==type(0):
         self.hwrBoundary=[self.hwrBoundary]*2
      self.hwrOverlap=self.config.getVal( "hwrOverlap",
            default=0.1,raiseerror=0) 
         # \/ n.b. Ought this not be supported by the util.dmInfo object
         #  but then that seems to be redundant code now?
      self.decayFactor=self.config.getVal( "decayFactor",
            default=1,raiseerror=0) 
      self.hwrSparse=self.config.getVal( "hwrSparse",
            default=False,raiseerror=0) 
      self.hwrGradientMgmnt=self.config.getVal( "hwrGradientMgmnt",
            default=False,raiseerror=0) 
      self.hwrVerbose=self.config.getVal( "hwrVerbose",
            default=False,raiseerror=0)
      self.hwrArchive=self.config.getVal( "hwrArchive",
            default=False,raiseerror=0)
      self.hwrLoadPrevious=self.config.getVal( "hwrLoadPrevious",
            default=0,raiseerror=0)
      self.hwrRecomputeWFIMandWFMM=self.config.getVal(
            "hwrRecomputeWFIMandWFMM", default=0,raiseerror=0)
      self.hwrWFIMfname=self.config.getVal( "hwrWFIMfname",
            default="hwrWFIM.fits", raiseerror=0 )
      self.hwrWFMMrcon=self.config.getVal( "hwrWFMMrcon",
            default=0.1,raiseerror=0)
      self.hwrWFMMblockReduction=self.config.getVal( "hwrWFMMblockReduction",
            default=None,raiseerror=0)
      self.hwrWFMMfname=self.config.getVal( "hwrWFMMfname",
            default="hwrWFMM.fits", raiseerror=0 )
      if self.hwrSparse: # only bother if using sparse implementation
         self.hwrCGtol=self.config.getVal( "hwrCGtol",
               default=1e-7,raiseerror=0) 
      self.pokeIgnore=self.config.getVal( "pokeIgnore",
            default=1e-10,raiseerror=0)
      
         # \/ Run a few checks to make sure we're in a supported configuration
      if len(self.dmList)!=1:
         raise Exception("ERROR: HWR: Can only support one DM")
      for dm in self.dmList:
         if not dm.zonalDM:
            raise Exception("ERROR: HWR: Can only support zonal DMs")
      if len(self.ngsList)!=1:
         raise Exception("ERROR: HWR: Can only support one NGS")
      if len(self.lgsList)!=0:
         raise Exception("ERROR: HWR: Cannot support any LGS")
      if self.hwrSparse:
         print("INFORMATION(**HWR**): Forcing sparsePmxType to become csr")
         self.sparsePmxType="csr"
         self.ndata=int(
            self.config.getVal("spmxNdata",default=self.nmodes*self.ncents*0.1))
      else:
         pass # hmm, perhaps something ought to go here but its after lunch and I can't quite recall what that something should be...

      nsubx_tmp=self.config.getVal("wfs_nsubx")
      #
      subapMap=self.pupil.getSubapFlag(nsubx_tmp, self.minarea) 
      self.gradOp=abbot.gradientOperator.gradientOperatorType1(
            subapMap )
#??      #
#??      #  \/ This maps neatly to the actual DM actuator locations but
#??      #  doesn't lead to a sufficient mapping: better to specify the
#??      #  same actual sub-apertures so that the wavefront is more accurately
#??      #  estimated on the grid, and then transferred to the actuators.
#??      DMactMap=numpy.array( self.dmList[0].dmflag )
#??      self.gradOp=abbot.gradientOperator.gradientOperatorType1(
#??            pupilMask=DMactMap )
# (redundant?) :      self.gradM=self.gradOp.returnOp().T # *** DEBUG ***
      print(
"INFORMATION(**HWR**): gradOp.(numberPhases={0:d},numberSubaps={1:d})".format(
            self.gradOp.numberPhases, self.gradOp.numberSubaps))
      t1=time.time()
      self.smmtnsDef,self.smmtnsDefStrts,self.smmtnsMap,self.offsetEstM=\
            abbot.hwr.prepHWR( self.gradOp,
               self.hwrMaxLen,
               self.hwrBoundary,
               self.hwrOverlap,
               self.hwrSparse
               )
      print("INFORMATION(**HWR**): took {0:d}s to setup HWR".format(
            int(time.time()-t1)) )
      
      if self.hwrWFMMblockReduction:
         # \/ calculate which indices should be included in the WFIM
         # first, do it for the wavefront grid
         self.WFIMbrIndicesWFGrid=[
            numpy.flatnonzero( (self.gradOp.illuminatedCornersIdx//
                  (self.gradOp.n_[1]*self.hwrWFMMblockReduction)) ==i)
               for i in range( int(numpy.ceil(
                    self.gradOp.n_[0]*self.hwrWFMMblockReduction**-1.0)) ) ]
         # now, do it for for every DM mode
         print("INFORMATION(**HWR**): block sizes, "+str(
               [ len(x) for x in self.WFIMbrIndicesWFGrid ] ))
         dmIdx=numpy.array(self.dmList[0].dmflag).ravel().nonzero()[0]
         self.WFIMbrIndicesModes=[
            numpy.flatnonzero( (dmIdx//
                  (self.gradOp.n_[1]*self.hwrWFMMblockReduction))==i)
               for i in range( int(numpy.ceil(
                    self.gradOp.n_[0]*self.hwrWFMMblockReduction**-1.0)) ) ]


   def calcTakeRef(self):
      self.control["takeRef"]=0
      self.control["zero_dm"]=1
      self.control["close_dm"]=0
      self.takingRef=1

   def calcLoadPreviousData(self):
      status=0
      try:
         if self.hwrSparse:
            self.spmx=util.FITS.loadSparse(self.pmxFilename)
            status=1
            self.WFIM=util.FITS.loadSparse(self.hwrWFIMfname)
            self.WFMM=( util.FITS.loadSparse(self.hwrWFMMfname,0),
                        util.FITS.loadSparse(self.hwrWFMMfname,1), )
            status=2
            if (self.WFMM[0].shape!=(self.nmodes,self.nmodes) or
                  self.WFMM[1].shape!=(self.nmodes,self.gradOp.numberPhases)):
               raise ValueError("Wrong shapes of WFMM (sparse)")
               status=1
         else:
            errmsg=""
            self.spmx=util.FITS.Read(self.pmxFilename)[1]
            status=1
            self.WFIM=util.FITS.Read(self.hwrWFIMfname)[1]
            self.WFMM=util.FITS.Read(self.hwrWFMMfname)[1]
            status=2
            if self.WFMM.shape!=(self.nmodes,self.gradOp.numberPhases):
               errmsg+="Wrong shape of WFMM (dense)"
               status=1
            if self.spmx.shape!=(self.nmodes,self.ncents):
               errmsg+="Wrong shape of PMX (dense)---was a sparse one loaded?"
               status=0
            if errmsg:
               raise ValueError(errmsg)
            else:
               status=2
         # now check other loaded matrices for shape
         if (self.spmx.shape!=(self.nmodes,self.ncents) or
               self.WFIM.shape!=(self.nmodes,self.gradOp.numberPhases)):
            status=0 # drop back to nowt
            raise ValueError("Wrong shapes of spmx, WFIM, or WFMM")
      except:
         print("WARNING:(**HWR**): Failure to load previous data,")
         print("WARNING:(**HWR**):  sys.exc_info()[0]={0:s}".format(
               str(sys.exc_info()[0])) )
         print("WARNING:(**HWR**):  sys.exc_info()[1]={0:s}".format(
               str(sys.exc_info()[1])) )
      return status

   def calcPrepareForPoke(self):
      self.control["poke"]=0
      self.pokeStartClock=time.clock()
      self.pokeStartTime=time.time()
      self.poking=1
      self.control["close_dm"]=0
      self.pokingDMNo=0#the current DM being poked
      self.pokingActNo=0#the current actuator(s) being poked
      self.pokingDMNoLast=0#the current DM being read by centroids
      self.pokingActNoLast=0#the current actuator(s) being read by cents
      dm=self.dmList[0]
      if dm.pokeSpacing!=None:
          self.pokeActMap=numpy.zeros((dm.nact,dm.nact),numpy.int32)
      if not self.hwrSparse:
          print("INFORMATION:(**HWR**): Dense pmx, type={0:s}, "+
               "[{1:d},{2:d}],[{3:d},{4:d}]".format(
                  str(self.reconDtype),self.nmodes,
                  self.gradOp.numberPhases,self.ncents,
                  self.gradOp.numberSubaps*2))
          self.spmx=numpy.zeros(
               (self.nmodes,self.ncents),numpy.float64)
      else: 
          print("INFORMATION:(**HWR**): Sparse poke matrix, forcing csc")
          self.spmxData=[]
          self.spmxRowind=[]
          self.spmxIndptr=[]
      print("INFORMATION:(**HWR**): Poking for {0:d} iterations".format(
            self.npokes+1))

   def calcSetPoke(self):
      if self.dmModeType=="poke":
          self.outputData[:,]=0.
          if self.poking<=self.npokes-self.totalLowOrderModalModes:
              dm=self.dmList[0]
              if dm.pokeSpacing!=None:
                  # poking several actuators at once
                  raise Exception("ERROR:(**HWR**): Not supported")
   #!!!                        self.pokeActMap[:]=0
   #!!!                        for i in range(self.pokingActNo/dm.pokeSpacing,
   #!!!                              dm.nact,dm.pokeSpacing):
   #!!!                           numpy.put(self.pokeActMap[i],
   #!!!                                 range(self.pokingActNo%dm.pokeSpacing,
   #!!!                                 dm.nact,dm.pokeSpacing),1)
   #!!!                        self.pokeActMapFifo.append(self.pokeActMap*dm.dmflag)
   #!!!                        self.outputData[
   #!!!                                 self.nactsCumList[self.pokingDMNo]:
   #!!!                                 self.nactsCumList[self.pokingDMNo+1]]=\
   #!!!                              numpy.take(self.pokeActMap.ravel(),
   #!!!                                         numpy.nonzero(dm.dmflag.ravel())
   #!!!                                        )[0]*self.pokeval
              else:
                  self.outputData[ self.nactsCumList[self.pokingDMNo]
                                  +self.pokingActNo 
                                 ]=self.pokeval
              self.pokingActNo+=1

   def calcFillInInteractionMatrix(self):
      # begin to fill the poke matrix
      # find out which DM we've just poked, and whether it was poked
      # individually or several actuators at once:
      dm=self.dmList[0]
      if dm.pokeSpacing!=None:
         # poking several actuators at once
         raise Exception("ERROR:(**HWR**): Not supported")
#!!!               pokenumber=None
#!!!               #Decide which actuator has produced which centroid values (if any), and then fill the poke matrix.
#!!!               self.fillPokemx(dm,self.pokingDMNoLast)
      else: # poked 1 at once so all centroids belong to this actuator
         pokenumber=self.nactsCumList[self.pokingDMNoLast]\
               +self.pokingActNoLast
      self.pokingActNoLast+=1

      if pokenumber!=None:
         this_spmx_data=self.inputData/self.pokeval
         if self.hwrSparse:
            for i,val in enumerate(this_spmx_data):
               if numpy.fabs(val)>self.pokeIgnore:
                  self.spmxData.append( val )
                  self.spmxRowind.append( pokenumber )
                  self.spmxIndptr.append( i )
         else:
            self.spmx[pokenumber]=this_spmx_data

   def calcCompletePoke(self):
      if self.hwrVerbose:
         print("INFORMATION(**HWR**)(verbose): Poking took: "+
               "{0:g}s/{1:g}s in CPU/wall-clock time".format(
                  time.clock()-self.pokeStartClock,
                  time.time()-self.pokeStartTime)
               )
      self.outputData[:,]=0.
      if self.hwrSparse:
         self.spmx=scipy.sparse.csc_matrix(
            (self.spmxData,
             (self.spmxRowind,
             self.spmxIndptr)
            ),(self.nmodes,self.ncents) )
      if self.pmxFilename!=None:
         print("INFORMATION:(**HWR**): Saving poke matrix to '{0:s}'".format(
               self.pmxFilename))
         if self.hwrSparse:
            util.FITS.saveSparse(self.spmx,self.pmxFilename)
         else:
            util.FITS.Write(self.spmx,self.pmxFilename)
      
      self.poking=0


   def calcComputeWFIMandWFMM(self):
      '''From spmx, compute the WFIM (wavefront-interaction matrix) and then
      compute the WFMM (wavefront mapping matrix).
      '''

         # \/ prep
      if self.hwrSparse:
         self.WFIMData=[]
         self.WFIMRowind=[]
         self.WFIMIndptr=[]
      else:
         self.WFIM=numpy.zeros(
               (self.nmodes,self.gradOp.numberPhases), numpy.float64 )
         # \/ do computations
      self.calcComputeWFIM()
      self.calcComputeWFMM()
         # \/ save, if asked to, the WFIM and WFMM  
      if self.hwrArchive:
         print(
              "INFORMATION:(**HWR**): Writing WFIM & WFMM to file {0:s}".format(
               self.hwrWFIMfname))
         if self.hwrSparse:
            util.FITS.saveSparse(self.WFIM,self.hwrWFIMfname)
               # \/ save the 2 matrices separately
            util.FITS.saveSparse(self.WFMM[0],self.hwrWFMMfname)
            util.FITS.saveSparse(self.WFMM[1],self.hwrWFMMfname,writeMode="a")
         else:
            util.FITS.Write(self.WFIM,self.hwrWFIMfname)
            util.FITS.Write(self.WFMM,self.hwrWFMMfname)

   def calcComputeWFIM(self):
      '''Compute the WFIM based on spmx'''
      if self.hwrVerbose:
         print("INFORMATION:(**HWR**): Calculating WFIM, wavefront-interaction matrix")

      # \/ loop over rows and convert
      for actNo,thisIp in enumerate(self.spmx):
         if self.hwrSparse:
            thisIp=numpy.array(thisIp.todense()).ravel()
         HWRcalc= abbot.hwr.doHWRGeneral( thisIp,
               self.smmtnsDef,self.gradOp,self.offsetEstM,
               self.smmtnsDefStrts,self.smmtnsMap,
               doWaffleReduction=0, doPistonReduction=1,
               sparse=self.hwrSparse )[1]
         if self.hwrWFMMblockReduction:
            thisi=None
            for i,tWFIMbrIndicesModes in enumerate(self.WFIMbrIndicesModes):
               # search through indices to find the right one
               if actNo in tWFIMbrIndicesModes:
                  thisi=i
                  continue
            if type(thisi)==type(None):
               print(actNo)
               raise Exception("ERROR:(**HWR**): Block-reduction index location"
                     " failed, too few defined blocks?")
            validModes=self.WFIMbrIndicesModes[thisi]
         if self.hwrSparse:
            for i,val in enumerate(HWRcalc):
               if self.hwrWFMMblockReduction and (i not in validModes):
                  continue
               self.WFIMData.append( val )
               self.WFIMRowind.append( actNo )
               self.WFIMIndptr.append( i )
         else:
            if self.hwrWFMMblockReduction:
               self.WFIM[actNo][validModes]=HWRcalc[validModes]
            else:
               self.WFIM[actNo]=HWRcalc

      if self.hwrSparse:
         self.WFIM=scipy.sparse.csc_matrix(
               (self.WFIMData, (self.WFIMRowind, self.WFIMIndptr) ),
               (self.nmodes, self.gradOp.numberPhases) )
          
   def calcComputeWFMM(self):
      '''Compute the WFMM based on WFIM'''
      if self.hwrVerbose:
         print("INFORMATION(**HWR**)(verbose): Calculating WFMM, wavefront-mapping matrix")
      #
      # NOTA BENE:
      # The order of the WFIM matrix means that the transposes
      # are inverted from that you might expect in the following creation
      # of the WFMM.
      #
      if not self.hwrSparse:
         self.WFMM=scipy.linalg.inv( numpy.dot(self.WFIM,self.WFIM.T)+
               numpy.identity( self.nmodes )*self.hwrWFMMrcon ).dot(
               self.WFIM )
         self.reconmx=self.WFMM
      else:
         # \/ store the output ready for a CG loop, gets turned into
         # CSR format, because of transposes.
         self.WFMM=( self.WFIM.dot( self.WFIM.T )
                  +scipy.sparse.identity(self.nmodes)*self.hwrWFMMrcon,
               self.WFIM ) 

   def calcClosedLoop(self):
      """Apply reconstruction, and mapping"""
      try:
         self.integratedV=abbot.hwr.doHWRGeneral( self.inputData,
                  self.smmtnsDef,self.gradOp,self.offsetEstM,
                  self.smmtnsDefStrts,self.smmtnsMap,
                  doWaffleReduction=0, doPistonReduction=0,
	          doGradientMgmnt=self.hwrGradientMgmnt,
                  sparse=self.hwrSparse )[1]
      except:
         import traceback,pickle
         print("ERROR:(**HWR**): failed in HWR integration,")
         print("ERROR:(**HWR**):  {0:s}".format(
               traceback.format_exc() ))
         print("ERROR:(**HWR**): Will write /tmp/hwrTmp.pickle ...")
         # The following line could over-write the previous dump, perhaps
         # instead the behaviour should be to throw an exception if it already
         # exists.
         fhandle=open("/tmp/hwrTmp.pickle","w")
         p=pickle.Pickler(fhandle)
         p.dump({ 'self.hwrMaxLen': self.hwrMaxLen,
                  'self.hwrBoundary': self.hwrBoundary,
                  'self.hwrOverlap': self.hwrOverlap,
                  'self.hwrSparse': self.hwrSparse,
                  'self.hwrSparse': self.hwrSparse,
                  'subapMap': self.gradOp.subapMask,
                  'self.inputData': self.inputData})
         raise Exception("ERROR: HWR: ... written")
      if not self.hwrSparse:
         self.outputV=numpy.dot( self.WFMM, self.integratedV )
      else:
         outputV=spcg( self.WFMM[0], self.WFMM[1].dot(self.integratedV),
               tol=self.hwrCGtol )
         if outputV[1]==0:
            self.outputV=outputV[0]
         else:
            raise Exception("ERROR: HWR: Could not converge on mapping")
      self.outputData=self.decayFactor*self.outputData+self.gains*-self.outputV

   def calc(self):
      #
      # \/ This section is the preliminaries; references, prepare the
      #  interaction matrix, and compute the wavefront mapping matrix/matrices.
      #
      if self.control["takeRef"]:
         self.calcTakeRef()
      if self.hwrLoadPrevious:
         success=self.calcLoadPreviousData()
         self.control["poke"]=0
         if success==0:
            print("WARNING:(**HWR**): Loading was not successful, must poke")
            self.control["poke"]=1
         else:
            self.poking=0 # stop poking
            if success==1 or self.hwrRecomputeWFIMandWFMM:
               print("INFORMATION:(**HWR**): Loaded spmx, recomputed WFIM and WFMM.")
               # compute the loaded WFIM and WFMM
               self.calcComputeWFIMandWFMM()
            else:
               print("INFORMATION:(**HWR**): Loaded spmx, WFIM, and WFMM.")
            self.control["close_dm"]=1 # close the loop
         self.hwrLoadPrevious=0
      if self.control["poke"]: 
         self.calcPrepareForPoke()
      if self.control["zero_dm"]:
         self.control["zero_dm"]=0
         self.outputData[:,]=0.
         # \/ don't these two lines skip a stage?
      if self.takingRef==1: self.takingRef=2
      if self.takingRef==2: # the reference centroids now complete
         self.takingRef=0
         self.refCentroids=self.inputData.copy()
         if type(self.saveRefCentroids)==type(""):
             util.FITS.Write(self.refCentroids,self.saveRefCentroids)
             if self.hwrVerbose:
                print("INFORMATION(**HWR**): Saving reference centroids, {0:s}".format(
                     self.saveRefCentroids))
         else:
             if self.hwrVerbose:
                print("INFORMATION(**HWR**): Not saving reference centroids")
      if self.control["subtractRef"] and type(self.refCentroids)!=type(None):
         self.inputData-=self.refCentroids
      if self.poking>0 and self.poking<=self.npokes:
         self.calcSetPoke()
      if self.poking>1 and self.poking<=self.npokes+1:
         self.calcFillInInteractionMatrix() 
      if self.poking==self.npokes+1: # finished poking
         self.calcCompletePoke()
         self.calcComputeWFIMandWFMM()
         if self.abortAfterPoke:
             print("INFORMATION(**HWR**): Finished poking - aborting simulation")
             if Scientific.MPI.world.size==1:
                 Scientific.MPI.world.abort()
             else:
                 Scientific.MPI.world.abort(0)
      if self.poking>0: # still poking
         self.poking+=1

      #
      # \/ This section is when the loop is closed.
      #
      if self.control["close_dm"]:
         self.calcClosedLoop()


   def plottable(self,objname="$OBJ"):
      """Return a XML string which contains the commands to be sent
      over a socket to obtain certain data for plotting.  The $OBJ symbol
      will be replaced by the instance name of the object - e.g.
      if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""

      thisid=(self.idstr==None or self.idstr==[None]
                )*"({0:s})".format(self.idstr)
      op=""
      op+=('<plot title="%s Output data%s" '+
          'cmd="data=%s.outputData" ret="data" type="pylab" '+
          'when="rpt" palette="jet"/>\n')%(self.objID,thisid,objname)
      op+=('<plot title="%s Integrated vector %s" '+
          'cmd="data=%s.integratedV" ret="data" type="pylab" '+
          'when="rpt" palette="jet"/>\n')%(self.objID,thisid,objname)
      op+=('<plot title="%s Mapped integration %s" '+
          'cmd="data=%s.outputV" ret="data" type="pylab" '+
          'when="rpt" palette="jet"/>\n')%(self.objID,thisid,objname)
      op+=('<plot title="%s Use reference centroids%s" '+
          'ret="data" when="cmd" texttype="1" wintype="mainwindow">\n'+
          '<cmd>\ndata=%s.control["subtractRef"]\n'+
          '%s.control["subtractRef"]=1-data\n</cmd>\n'+
          'button=1-data\n</plot>\n')%(self.objID,thisid,objname,objname)

# (redundant?) :      op+=('<plot title="%s gradM%s" '+
# (redundant?) :          'cmd="data=%s.gradM" ret="data" type="pylab" '+
# (redundant?) :          'when="cmd"/>\n')%(self.objID,thisid,objname)
      if self.hwrSparse:
         thiscmd='"data=numpy.array(%s.spmx.todense())"'%(objname)
      else:
         thiscmd='"data=%s.spmx"'%(objname)
      op+=('<plot title="%s pmx%s" '+
           'cmd=%s ret="data" type="pylab" '+
           'when="cmd"/>\n')%(self.objID,thisid,thiscmd)
      if self.hwrSparse:
         thiscmd='"data=numpy.array(%s.WFIM.todense())"'%(objname)
      else:
         thiscmd='"data=%s.WFIM"'%(objname)
      op+=('<plot title="%s WFIM%s" '+
           'cmd=%s ret="data" type="pylab" '+
           'when="cmd"/>\n')%(self.objID,thisid,thiscmd)
      op+=('<plot title="%s inputData%s" '+
          'cmd="data=%s.inputData" ret="data" type="pylab" '+
          'when="rpt"/>\n')%(self.objID,thisid,objname)

      return op 
