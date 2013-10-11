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


import os.path
import numpy
import base.aobase
import tomoRecon
import util.FITS
import cmod.svd
#?? import cmod.utils
#?? import util.dot as quick

import scipy.sparse,scipy.linalg
from scipy.sparse.linalg import cg as spcg
#?? import util.spmatrix
import time,types
import Scientific.MPI

import util.gradientOperator
import abbot.hwr


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
      self.hwrSparse=self.config.getVal( "hwrSparse",
            default=False,raiseerror=0) 
      self.hwrVerbose=self.config.getVal( "hwrVerbose",
            default=False,raiseerror=0)
      self.hwrArchive=self.config.getVal( "hwrArchive",
            default=False,raiseerror=0)
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
         raise Exception("HWR: Can only support one DM")
      for dm in self.dmList:
         if not dm.zonalDM: raise Exception("HWR: Can only support zonal DMs")
      if len(self.ngsList)!=1:
         raise Exception("HWR: Can only support one NGS")
      if len(self.lgsList)!=0:
         raise Exception("HWR: Cannot support any LGS")
      if self.hwrSparse:
         print("HWR: Forcing sparsePmxType to become csr")
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
      self.gradM=self.gradOp.returnOp().T # *** DEBUG ***
      print("HWR: gradM shape is [{0[0]:d},{0[1]:d}]".format(self.gradM.shape))
      self.smmtnsDef,self.smmtnsDefStrts,self.smmtnsMap,self.offsetEstM=\
            abbot.hwr.prepHWR( self.gradOp,
               self.hwrMaxLen,
               self.hwrBoundary,
               self.hwrOverlap,
               self.hwrSparse
               )
      
      if self.hwrWFMMblockReduction:
         # \/ calculate which indices should be included in the WFIM
         # first, do it for the wavefront grid
         self.WFIMbrIndicesWFGrid=[
            numpy.flatnonzero( (self.gradOp.illuminatedCornersIdx//
                  (self.gradOp.n_[1]*self.hwrWFMMblockReduction)) ==i)
               for i in range(
                    self.gradOp.n_[0]//self.hwrWFMMblockReduction)]
         # now, do it for for every DM mode
         dmIdx=numpy.array(self.dmList[0].dmflag).ravel().nonzero()[0]
         self.WFIMbrIndicesModes=[
            numpy.flatnonzero( (dmIdx//
                  (self.gradOp.n_[1]*self.hwrWFMMblockReduction))==i)
               for i in range(
                    self.gradOp.n_[0]//self.hwrWFMMblockReduction)]


   def calcTakeRef(self):
      self.control["takeRef"]=0
      self.control["zero_dm"]=1
      self.control["close_dm"]=0
      self.takingRef=1

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
          print("HWR: Dense pmx, type={0:s}, [{1:d},{2:d}],[{3:d},{4:d}]".format(
               str(self.reconDtype),self.nmodes,
               self.gradOp.numberPhases,self.ncents,
               self.gradOp.numberSubaps*2))
          self.spmx=numpy.zeros(
               (self.nmodes,self.ncents),numpy.float64)
          self.WFIM=numpy.zeros(
               (self.nmodes,self.gradOp.numberPhases),numpy.float64)
      else: 
          print("HWR: Sparse poke matrix, forcing csc")
#??                self.spmxData=numpy.zeros((self.ndata,),numpy.float32)
#??                self.spmxRowind=numpy.zeros((self.ndata,),numpy.int32)
#??                self.spmxIndptr=numpy.zeros((self.ncents+1,),numpy.int32)
          self.spmxData=[]
          self.spmxRowind=[]
          self.spmxIndptr=[]
#??                self.spmx=cmod.svd.sparseMatrixCreate(
#??                     self.ndata,self.nmodes,self.ncents,
#??                     self.spmxData,self.spmxRowind,self.spmxIndptr)
#??                self.WFIMData=numpy.zeros((self.ndata,),numpy.float32)
#??                self.WFIMRowind=numpy.zeros((self.ndata,),numpy.int32)
#??                self.WFIMIndptr=numpy.zeros((self.gradOp.numberPhases+1,),
#??                     numpy.int32)
          self.WFIMData=[]
          self.WFIMRowind=[]
          self.WFIMIndptr=[]
#??                self.WFIM=cmod.svd.sparseMatrixCreate(
#??                     self.ndata,self.nmodes,self.gradOp.numberPhases,
#??                     self.WFIMData,self.WFIMRowind,self.WFIMIndptr)
      print("HWR: Poking for {0:d} iterations".format(self.npokes+1))

   def calcSetPoke(self):
      if self.dmModeType=="poke":
          self.outputData[:,]=0.
          if self.poking<=self.npokes-self.totalLowOrderModalModes:
              dm=self.dmList[0]
              if dm.pokeSpacing!=None:
                  # poking several actuators at once
                  raise Exception("HWR: Not supported")
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
         raise Exception("HWR: Not supported")
#!!!               pokenumber=None
#!!!               #Decide which actuator has produced which centroid values (if any), and then fill the poke matrix.
#!!!               self.fillPokemx(dm,self.pokingDMNoLast)
      else: # poked 1 at once so all centroids belong to this actuator
         pokenumber=self.nactsCumList[self.pokingDMNoLast]\
               +self.pokingActNoLast
      self.pokingActNoLast+=1

      if pokenumber!=None:
         if self.hwrWFMMblockReduction:
            thisi=None
            for i in xrange(int(self.gradOp.n_[0]//self.hwrWFMMblockReduction)):
               # search through indices to find the right one
               if pokenumber in self.WFIMbrIndicesWFGrid[i]:
                  thisi=i
                  continue
            if type(thisi)==type(None):
               print(pokenumber)
               raise Exception("HWR: Block-reduction index location failed, too few defined blocks?")
            validModes=self.WFIMbrIndicesModes[thisi]
         this_spmx_data=self.inputData/self.pokeval
         HWRcalc= abbot.hwr.doHWRGeneral( this_spmx_data,
               self.smmtnsDef,self.gradOp,self.offsetEstM,
               self.smmtnsDefStrts,self.smmtnsMap,
               doWaffleReduction=0, doPistonReduction=0,
               sparse=self.hwrSparse )[1]
         if self.hwrSparse:
            for i,val in enumerate(this_spmx_data):
               if numpy.fabs(val)>self.pokeIgnore:
                  self.spmxData.append( val )
                  self.spmxRowind.append( pokenumber )
                  self.spmxIndptr.append( i )
#??                        cmod.svd.sparseMatrixInsert(
#??                              self.spmx, pokenumber,i, val )
            for i,val in enumerate(HWRcalc):
               if self.hwrWFMMblockReduction and (i not in validModes):
                  continue
               self.WFIMData.append( val )
               self.WFIMRowind.append( pokenumber )
               self.WFIMIndptr.append( i )
#??                     cmod.svd.sparseMatrixInsert(
#??                           self.WFIM, pokenumber,i, val )
         else:
            self.spmx[pokenumber]=this_spmx_data
            if self.hwrWFMMblockReduction:
               self.WFIM[pokenumber][validModes]=HWRcalc[validModes]
            else:
               self.WFIM[pokenumber]=HWRcalc

   def calcCompletePoke(self):
      self.outputData[:,]=0.
      if self.hwrSparse:
#??               converted_spmx=scipy.sparse.csc_matrix(
         self.spmx=scipy.sparse.csc_matrix(
            (self.spmxData,
             (self.spmxRowind,
             self.spmxIndptr)
            ),(self.nmodes,self.ncents) )
#??               cmod.svd.sparseMatrixFree(self.spmx)
#??               self.spmx=converted_spmx.T # switch to prefered format
         self.WFIM=scipy.sparse.csc_matrix(
            (self.WFIMData,
             (self.WFIMRowind,
             self.WFIMIndptr)
            ),(self.nmodes,self.gradOp.numberPhases) )
      print("HWR: Finished poking")
      if self.hwrArchive:
         print("HWR: Writing WFIM to file {0:s}".format(self.hwrWFIMfname))
         if self.hwrSparse:
            util.FITS.Write(self.WFIM.todense(),self.hwrWFIMfname)
         else:
            util.FITS.Write(self.WFIM,self.hwrWFIMfname)
          
      if not self.computeControl:
         raise Exception("HWR: You must set the computeControl to true")
      print("HWR: Calculating wavefront-interaction matrix")
            
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
         if self.hwrArchive:
            print("HWR: Writing WFMM to file '{0:s}'".format(self.hwrWFMMfname))
            util.FITS.Write(self.WFMM,self.hwrWFMMfname)
            if self.pmxFilename!=None:
               print("HWR: Saving poke matrix to '{0:s}'".format(self.pmxFilename))
               util.FITS.Write(self.spmx,self.pmxFilename)
      else:
         # \/ store the output ready for a CG loop, gets turned into
         # CSR format, because of transposes.
         self.WFMM=(
               self.WFIM.dot( self.WFIM.T )
                  +scipy.sparse.identity(self.nmodes)
                  *self.hwrWFMMrcon,
               self.WFIM ) 
         self.reconmx=None
         if self.hwrArchive:
            print("HWR: Writing WFMM to file '{0:s}'".format(self.hwrWFMMfname))
            util.FITS.Write(self.WFMM[0].todense(),self.hwrWFMMfname)
            util.FITS.Write(self.WFMM[1].todense(),self.hwrWFMMfname,writeMode="a")
            if self.pmxFilename!=None:
               print("HWR: Saving poke matrix to '{0:s}'".format(self.pmxFilename))
               util.FITS.Write(self.spmx.todense(),self.pmxFilename)
      if self.hwrVerbose:
         print("HWR: Poking took: "+
               "{0:g}s/{1:g}s in CPU/wall-clock time".format(
                  time.clock()-self.pokeStartClock,
                  time.time()-self.pokeStartTime)
               )
      self.poking=0
      if self.abortAfterPoke:
          print("HWR: Finished poking - aborting simulation")
          if Scientific.MPI.world.size==1:
              Scientific.MPI.world.abort()
          else:
              Scientific.MPI.world.abort(0)

   def calcClosedLoop(self):
      """Apply reconstruction, and mapping"""
      integratedV=abbot.hwr.doHWRGeneral( self.inputData,
               self.smmtnsDef,self.gradOp,self.offsetEstM,
               self.smmtnsDefStrts,self.smmtnsMap,
               doWaffleReduction=0, doPistonReduction=0,
               sparse=self.hwrSparse )[1]
      if not self.hwrSparse:
         outputV=numpy.dot( self.WFMM, integratedV )
      else:
         outputV=spcg( self.WFMM[0], self.WFMM[1].dot(integratedV),
               tol=self.hwrCGtol )
         if outputV[1]==0:
            outputV=outputV[0]
         else:
            raise Exception("HWR: Could not converge on mapping")
      self.outputData+=self.gains*-outputV

   def calc(self):
      #
      # \/ This section is the preliminaries; references, prepare the
      #  interaction matrix, and compute the wavefront mapping matrix/matrices.
      #
      if self.control["takeRef"]:
         self.calcTakeRef()
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
                print("HWR: Saving reference centroids, {0:s}".format(
                     self.saveRefCentroids))
         else:
             if self.hwrVerbose:
                print("HWR: Not saving reference centroids")
      if self.control["subtractRef"] and type(self.refCentroids)!=type(None):
         self.inputData-=self.refCentroids
      if self.poking>0 and self.poking<=self.npokes:
         self.calcSetPoke()
      if self.poking>1 and self.poking<=self.npokes+1:
         self.calcFillInInteractionMatrix() 
      if self.poking==self.npokes+1:#finished poking
         self.calcCompletePoke()
      if self.poking>0:
         self.poking+=1

      #
      # \/ This section is when the loop is closed.
      #
      if self.control["close_dm"]:
         self.calcClosedLoop()

#!!!    def fillPokemx(self,dm,dmindx):
#!!!        """Here, when we've been poking multiple actuators at once, we need to decide which centroids \
#!!!           belong to which actuator, and then place them into the poke matrix.
#!!!        """
#!!!        #First compute the poked actuator coords for this DM.
#!!!        if not hasattr(dm,"coords"):
#!!!            dm.computeCoords(self.telDiam)
#!!!        actuators=self.pokeActMapFifo.pop(0)
#!!!        nonzero=numpy.nonzero(actuators.ravel())[0]
#!!!        if nonzero.size==0:
#!!!            print "Got no actuators used for this poke"
#!!!            return
#!!!        pokedCoordsx=numpy.take(dm.coords.ravel(),nonzero*2)#coordinates of poked ones only.
#!!!        pokedCoordsy=numpy.take(dm.coords.ravel(),nonzero*2+1)#coordinates of poked ones only.
#!!!        cpos=0
#!!!
#!!!        cnt=0
#!!!        gsList=self.ngsList+self.lgsList
#!!!        for i in range(len(self.wfsIDList)):#for each wfs...
#!!!            gs=gsList[i]
#!!!            gs.computeCoords(self.telDiam,dm.height)
#!!!            key=self.wfsIDList[i]
#!!!            ns=self.ncentList[i]
#!!!            for subap in self.centIndex[cnt:cnt+ns]/2:#for each subap of this wfs
#!!!                x=subap%gs.nsubx
#!!!                y=subap/gs.nsubx
#!!!                dist=(pokedCoordsx-gs.coords[y,x,0])**2+(pokedCoordsy-gs.coords[y,x,1])**2
#!!!                m=numpy.argmin(dist)#this actuator is closest.
#!!!                # find the coords for this actuator, and then the number that it corresponds to in the poke matrix.
#!!!                indx=nonzero[m]
#!!!                pos=dm.dmflag.ravel()[:indx].sum()+self.nactsCumList[dmindx]
#!!!                #if pos==196+65 and (cpos in [145,146]):
#!!!                #    print (m,indx,pos,x,y,subap,cpos,dist,nonzero)
#!!!                #    print self.nactsCumList[dmindx]
#!!!                #    print dm.dmflag.shape
#!!!                #    print pokedCoordsx
#!!!                #    print pokedCoordsy
#!!!                #    print gs.coords[y,x]
#!!!                #    #util.FITS.Write(dm.dmflag,"dmflag.fits")
#!!!                #    if cpos==146:
#!!!                #        util.FITS.Write(actuators,"acts.fits")
#!!!                # pos is the actuator number in the outputData for this DM
#!!!                # This is where the centroid should go in the poke matrix.
#!!!                if self.reconType in ["svd","pinv","MAP","reg","regularised","regBig","regSmall","pcg"]:
#!!!                    self.spmx[pos,cpos]=self.inputData[cpos]/self.pokeval
#!!!                    self.spmx[pos,cpos+self.ncents/2]=self.inputData[cpos+self.ncents/2]/self.pokeval
#!!!                elif self.reconType in ["spmx","spmxSVD","spmxGI"]:
#!!!                    if self.sparsePmxType=="lil":
#!!!                        val=self.inputData[cpos]
#!!!                        if numpy.fabs(val)>self.pokeIgnore:
#!!!                            self.spmx[pos,cpos]=val/self.pokeval
#!!!                        val=self.inputData[cpos+self.ncents/2]
#!!!                        if numpy.fabs(val)>self.pokeIgnore:
#!!!                            self.spmx[pos,cpos+self.ncents/2]=val/self.pokeval
#!!!                    elif self.sparsePmxType=="csc":
#!!!                        val=self.inputData[cpos]
#!!!                        if numpy.fabs(val)>self.pokeIgnore:
#!!!                            cmod.svd.sparseMatrixInsert(self.spmx,pos,cpos,val/self.pokeval)
#!!!                        val=self.inputData[cpos+self.ncents/2]
#!!!                        if numpy.fabs(val)>self.pokeIgnore:
#!!!                            cmod.svd.sparseMatrixInsert(self.spmx,pos,cpos+self.ncents/2,val/self.pokeval)
#!!!                    else:
#!!!                        val=self.inputData[cpos]
#!!!                        val2=self.inputData[cpos+self.ncents/2]
#!!!                        if numpy.fabs(val)>self.pokeIgnore:
#!!!                            if self.transposeDensePmx==0:
#!!!                                self.spmx[pos,cpos]=val/self.pokeval
#!!!                            else:
#!!!                                self.spmx[cpos,pos]=val/self.pokeval
#!!!                            self.spmxValidCnt+=1
#!!!                        if numpy.fabs(val2)>self.pokeIgnore:
#!!!                            if self.transposeDensePmx==0:
#!!!                                self.spmx[pos,cpos+self.ncents/2]=val2/self.pokeval
#!!!                            else:
#!!!                                self.spmx[cpos+self.ncents/2,pos]=val2/self.pokeval
#!!!                            self.spmxValidCnt+=1
#!!!                cpos+=1
#!!!            cnt+=ns

#!!!    def doCompressedDot(self,data,rmx=None,bits=None,shape=None,work=None):
#!!!        """Here, rmx is a 1D array in compressed float format, with bits bits per element.
#!!!        """
#!!!        if rmx==None:
#!!!            rmx=self.reconmx
#!!!        if bits==None:
#!!!            bits=self.compressedBits
#!!!        if shape==None:
#!!!            shape=self.compressedShape
#!!!        if work==None:
#!!!            work=self.compressedWork
#!!!        expMin=self.compressedExpMin
#!!!        expMax=self.compressedExpMax
#!!!        if shape[0]==data.shape[0]:
#!!!            ncents,nacts=shape
#!!!            if work==None:
#!!!                mem=getMem("MemFree:")
#!!!                nrows=min(mem/4/nacts,ncents)
#!!!                print "Compressed rmx decompressing to %d rowsa at a time"%nrows
#!!!                work=numpy.zeros((nrows,nacts),numpy.float32)
#!!!                self.compressedWork=work
#!!!            nrows=work.shape[0]
#!!!            tmp=numpy.zeros((nacts,),numpy.float32)
#!!!            r=work.ravel()
#!!!            nsteps=(ncents+nrows-1)/nrows
#!!!            for i in xrange(nsteps):
#!!!                #uncompress into work
#!!!                #if (i+1)*nrows>=ncents:#this is the last one...
#!!!                if i==nsteps-1:#this is the last one...
#!!!                    end=nrows-(nsteps*nrows-ncents)
#!!!                else:
#!!!                    end=nrows
#!!!                #print "running %d"%i
#!!!                cmod.utils.uncompressFloatArrayAllThreaded(rmx,r[:end*nacts],bits,expMin,expMax,i*nrows*nacts,8)
#!!!                #print "done"
#!!!                tmp+=quick.dot(data[i*nrows:i*nrows+end],work[:end])
#!!!        else:
#!!!            nacts,ncents=shape
#!!!            if work==None:
#!!!                mem=getMem("MemFree:")
#!!!                nrows=min(mem/4/ncents,nacts)
#!!!                print "Compressed rmx decompressing to %d rowsb at a time"%nrows
#!!!                work=numpy.zeros((nrows,ncents),numpy.float32)
#!!!                self.compressedWork=work
#!!!            nrows=work.shape[0]
#!!!            tmp=numpy.empty((nacts,),numpy.float32)
#!!!            r=work.ravel()
#!!!            nsteps=(nacts+nrows-1)/nrows
#!!!            for i in xrange(nsteps):
#!!!                if i==nsteps-1:#last one...
#!!!                    end=nrows-(nsteps*nrows-nacts)
#!!!                else:
#!!!                    end=nrows
#!!!                #print "running %d"%i
#!!!                cmod.utils.uncompressFloatArrayAllThreaded(rmx,r[:end*ncents],bits,expMin,expMax,i*nrows*ncents,8)
#!!!                #print "done"
#!!!                tmp[i*nrows:i*nrows+end]=quick.dot(work[:end],data)
#!!!        return tmp

#!!!    def savecsc(self,csc,filename,hdr=None):
#!!!        if type(hdr)==type(""):
#!!!            hdr=[hdr]
#!!!        elif type(hdr)==type(None):
#!!!            hdr=[]
#!!!        hdr.append("SHAPE   = %s"%str(csc.shape))
#!!!        util.FITS.Write(csc.data[:csc.indptr[-1]],filename,extraHeader=hdr)
#!!!        util.FITS.Write(csc.rowind[:csc.indptr[-1]],filename,writeMode="a")
#!!!        util.FITS.Write(csc.indptr,filename,writeMode="a")
#!!!                
#!!!    def loadReconmx(self,reconmxFilename):
#!!!        if reconmxFilename==None:
#!!!            print "Reconmxfilename not specified - using 0."
#!!!            return 0.
#!!!        print "tomoRecon: Loading reconstructor from file: %s"%reconmxFilename
#!!!        if os.path.exists(reconmxFilename):
#!!!            head=util.FITS.ReadHeader(reconmxFilename)["parsed"]
#!!!            if head.has_key("COMPBITS"):
#!!!                print "Reconstructor in compressed FP format"
#!!!                #its a compressed floating point format rmx...
#!!!                reconmx=util.FITS.Read(reconmxFilename)[1]
#!!!                self.compressedBits=int(head["COMPBITS"])
#!!!                self.compressedExpMin=int(head.get("EXPMIN",1))
#!!!                self.compressedExpMax=int(head.get("EXPMAX",255))
#!!!                self.compressedShape=eval(head["SHAPE"])
#!!!            else:
#!!!                # f=util.FITS.Read(reconmxFilename,savespace=1)
#!!!                reconmx=util.FITS.loadSparse(reconmxFilename)
#!!!        else:
#!!!            print "Error/warning - unable to load reconmx %s, using 0. instead"%reconmxFilename
#!!!            reconmx=0.
#!!!            #f=[0.]
#!!!        #if len(f)==2:
#!!!        #    reconmx=f[1].astype("f")
#!!!        #elif len(f)>2:
#!!!        #    if len(f[1].shape)==1:
#!!!        #        reconmx=scipy.sparse.csc_matrix((numpy.array(f[1],copy=0),numpy.array(f[3],copy=0),numpy.array(f[5],copy=0)),eval(f[0]["parsed"]["SHAPE"]))
#!!!        #    else:#f[3,5,etc] are probably phase covariance/noise cov etc...
#!!!        #        reconmx=f[1]
#!!!        return reconmx

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
      op+=('<plot title="%s Use reference centroids%s" '+
          'ret="data" when="cmd" texttype="1" wintype="mainwindow">\n'+
          '<cmd>\ndata=%s.control["subtractRef"]\n'+
          '%s.control["subtractRef"]=1-data\n</cmd>\n'+
          'button=1-data\n</plot>\n')%(self.objID,thisid,objname,objname)

      op+=('<plot title="%s gradM%s" '+
          'cmd="data=%s.gradM" ret="data" type="pylab" '+
          'when="cmd"/>\n')%(self.objID,thisid,objname)
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
          'when="cmd"/>\n')%(self.objID,thisid,objname)

      return op 
