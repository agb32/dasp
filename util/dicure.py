#(deprecated):# What is this?
#(deprecated):# Two classes for dicure (diagonal cure) wavefront reconstruction:
#(deprecated):#
#(deprecated):# - gradientOperatorType1 : utility class to calculate the gradient operator
#(deprecated):#
#(deprecated):# - DiCuRe : main class
#(deprecated):#   - calc() : the method that performs the reconstruction
#(deprecated):# 
#(deprecated):# Author: NAB
#(deprecated):# Included into aosim: UB
#(deprecated):#

raise RuntimeError("DEPRECATED, use hwr.py instead")

#(deprecated):import numpy as np
#(deprecated):import util.dot as quick
#(deprecated):
#(deprecated):# 
#(deprecated):# UB, 2012 Aug 6th:
#(deprecated):# I copied (this part of) class gradientOperatorType1 from NAB's "gradientOperator.py".
#(deprecated):# I only copied what is needed by aosim and nothing else.
#(deprecated):# 
#(deprecated):# Hardy, 1998, p.270 for configurations
#(deprecated):# Type 1 = Fried
#(deprecated):class gradientOperatorType1:
#(deprecated):    '''Using Type 1 geometry, define a gradient operator matrix.'''
#(deprecated):    op=None
#(deprecated):    numberSubaps=None
#(deprecated):    numberPhases=None
#(deprecated):
#(deprecated):    def __init__( self, subapMask=None, pupilMask=None ):
#(deprecated):        if subapMask!=None: self.newSubaperturesGiven(subapMask)
#(deprecated):        if pupilMask!=None: self.newPupilGiven(pupilMask)
#(deprecated):
#(deprecated):    def newSubaperturesGiven(self, subapMask):
#(deprecated):        self.op=None # reset
#(deprecated):        self.n=subapMask.shape
#(deprecated):        self.n_=[ x+1 for x in self.n]
#(deprecated):
#(deprecated):      # now define the illuminated sub-apertures (-1 = not illuminated)
#(deprecated):        self.maskIdx=np.arange(self.n[0]*self.n[1]).reshape(self.n)*subapMask+(subapMask-1)
#(deprecated):         # \/ corners on the array
#(deprecated):        self.cornersIdx=np.arange(self.n_[0]*self.n_[1]).reshape(self.n_) 
#(deprecated):         # \/ per corner, how many subapertures use this corner
#(deprecated):        self.illuminatedCorners=np.zeros(self.n_)
#(deprecated): 
#(deprecated):      # \/ array to specify the first corner for each sub-aperture
#(deprecated):        self.mask2allcorners=\
#(deprecated):            np.flatnonzero(subapMask.ravel())%(self.n[1])\
#(deprecated):            +(np.flatnonzero(subapMask.ravel())//(self.n[1]))*self.n_[1]
#(deprecated):      # \/ figure out which phases then contribute to a sub-aperture gradient
#(deprecated):        for i in range(self.n[0]):
#(deprecated):            for j in range(self.n[1]):
#(deprecated):                if self.maskIdx[i,j]==-1: continue # skip not illuminated
#(deprecated):                self.illuminatedCorners+=(
#(deprecated):                    (abs(self.cornersIdx//self.n_[1]-0.5-self.maskIdx[i,j]//self.n[1])<=0.5)*\
#(deprecated):                        (abs(self.cornersIdx%self.n_[1]-0.5-self.maskIdx[i,j]%self.n[1])<=0.5) )
#(deprecated):        self.illuminatedCornersIdx=np.flatnonzero(self.illuminatedCorners.ravel())
#(deprecated):
#(deprecated):      # \/ some handy numbers
#(deprecated):        self.numberSubaps=int(subapMask.sum()) # =n**2 for all illuminated
#(deprecated):        self.numberPhases=(self.illuminatedCorners!=0).sum() # =(n+1)**2 for all illum.
#(deprecated):
#(deprecated):    def newPupilGiven(self, pupilMask):
#(deprecated):        '''Define the mask given a pupil mask, by creating a sub-aperture mask
#(deprecated):        and then re-creating the pupil'''
#(deprecated):
#(deprecated):        self.subapMask=((pupilMask[:-1,:-1]+pupilMask[1:,:-1]\
#(deprecated):                  +pupilMask[1:,1:]+pupilMask[:-1,1:])==4)
#(deprecated):        self.newSubaperturesGiven(self.subapMask)
#(deprecated):        print "DICURE SUBAPMASK:", self.subapMask
#(deprecated):
#(deprecated):
#(deprecated):#
#(deprecated):# Class to handle the dicure reconstruction of the wavefront:
#(deprecated):#
#(deprecated):class DiCuRe:
#(deprecated):#    def __init__(self, subapMap):
#(deprecated):    def __init__(self, actuatorMap):
#(deprecated):
#(deprecated):#        if subapMap.shape[0] != subapMap.shape[1]:
#(deprecated):#            raise Exception("DiCuRe can not be used with non-square subaperture maps.")
#(deprecated):        if actuatorMap.shape[0] != actuatorMap.shape[1]:
#(deprecated):            raise Exception("DiCuRe can not be used with non-square actuator maps.")
#(deprecated):
#(deprecated):        # Set up the geometry:
#(deprecated):#        self.gInst = gradientOperatorType1( subapMap )
#(deprecated):        self.gInst = gradientOperatorType1( pupilMask = actuatorMap )
#(deprecated):        # Define the chains:
#(deprecated):        self.chainsDef,\
#(deprecated):        self.chainsDefChStrts = self.chainsDefine() # chainsDefine returns 2 values
#(deprecated):        # the actuatorMask:
#(deprecated):        self.actuatorMap = actuatorMap
#(deprecated):        # Get the number of subapertures in one direction:
#(deprecated):#        self.nsubx = subapMap.shape[0]
#(deprecated):        self.nsubx = actuatorMap.shape[0]-1
#(deprecated):        # the number of actuators in one direction:
#(deprecated):        self.nactx = self.nsubx+1
#(deprecated):
#(deprecated):    def calc(self, gradV):
#(deprecated):        rgradV = self.rotateVectors(gradV).ravel()
#(deprecated):        chains = self.chainsIntegrate(rgradV)
#(deprecated):        chainsV,chainsVOffsets = self.chainsVectorize(chains)
#(deprecated):        chainsOvlps = self.chainsOverlaps()
#(deprecated):        A,B = self.chainsDefMatrices(chainsOvlps, chainsV, chainsVOffsets)
#(deprecated):        # ought to be able to remove 'chains' here, but been lazy further down
#(deprecated):        del chainsOvlps,chainsVOffsets
#(deprecated):        invchOScovM = np.identity(A.shape[1])
#(deprecated):        offsetEstM = quick.dot( quick.dot( np.linalg.inv( quick.dot(A.T,A)+invchOScovM*1e-3 ), A.T), -B )
#(deprecated):        offsetEstV = np.dot( offsetEstM, chainsV ) # 2012 Aug 03: currently you can not use quick.dot
#(deprecated):                                                   # here because chainsV is a list and not an array
#(deprecated):        comp = np.zeros(2*self.nactx**2, np.float64)
#(deprecated):
#(deprecated):        # do one way...
#(deprecated):        for x in range(len(chains[0])):
#(deprecated):            toeI=x
#(deprecated):            for i in range((self.chainsDef[0][x][1])):
#(deprecated):                comp[ self.chainsDef[0][x][0][i] ] += (chains[0][x][i]+offsetEstV[toeI])
#(deprecated):                pass
#(deprecated):        # ...then another
#(deprecated):        for x in range(len(chains[1])):
#(deprecated):            toeI=self.chainsDefChStrts[1][2][x]
#(deprecated):            for i in range((self.chainsDef[1][x][1])):
#(deprecated):                comp[ self.nactx**2+self.chainsDef[1][x][0][i] ] += (chains[1][x][i]+offsetEstV[toeI])
#(deprecated):                pass
#(deprecated):        np.set_printoptions(precision = 3)
#(deprecated):#        print "COMP:", comp.reshape((2, self.nactx, self.nactx))
#(deprecated):#        print "ACTUATORS:", self.actuatorMap, "VSOTA:", sum(sum(self.actuatorMap))
#(deprecated):#        print "ILLUMINATED CORNERS:", self.gInst.illuminatedCornersIdx
#(deprecated):#        print "size(ilum corners):", self.gInst.illuminatedCornersIdx.size
#(deprecated):        comp.resize(2, self.nactx**2) # separate chains 1 and chains 2
#(deprecated):        rs = comp.mean(0) # calculate mean over chains 1 and chains 2; rs is nactx**2 long
#(deprecated):#        print "RS:", rs
#(deprecated):        idx = self.actuatorMap.ravel().nonzero() # get the indices of active DM actuators
#(deprecated):#        print "IDX:", idx
#(deprecated):        result  = rs[idx] # get the active-actuators' values from rs to result
#(deprecated):        result -= result.mean() # subtract the mean value from all the elements
#(deprecated):#        print "RESULT:", result
#(deprecated):        return result # that's it
#(deprecated):
#(deprecated):    def rotateVectors(self, grads):
#(deprecated):        return np.array(
#(deprecated):            [ grads[:grads.shape[0]/2]*2**-0.5+grads[grads.shape[0]/2:]*2**-0.5, 
#(deprecated):              grads[:grads.shape[0]/2]*-2**-0.5+grads[grads.shape[0]/2:]*2**-0.5 ])
#(deprecated):
#(deprecated):    def chainsDefine(self, maxLen=None ):
#(deprecated):          # \/ define the chains
#(deprecated):          # Let the 1st method be to zip through the subaperture index into the
#(deprecated):          # pupil mask and find which sub-aperture corners are in the
#(deprecated):          # mask2allcorners and that this element is the bottom-left (thus there
#(deprecated):          # should also be four partners).
#(deprecated):          # Let the 2nd method be similar but mask2allcorners+n_ (top left).
#(deprecated):       N=self.gInst.n_ # alias
#(deprecated):       blank=np.zeros([2]+N)
#(deprecated):       chainsNumber=0
#(deprecated):       chaninsStrtNum=0
#(deprecated):       chainsDef=[],[]
#(deprecated):       chainsStarts=[ ([],[],[]), ([],[],[]), 0]
#(deprecated):       insubap=False
#(deprecated):       for i in range(N[1]+N[0]-1):
#(deprecated):          # top left, upwards right
#(deprecated):          y=int(np.where( N[0]-i-1<0, 0, N[0]-i-1 )) # starting y, then x
#(deprecated):          x=int(np.where( i-N[0]+1<0, 0, i-N[0]+1 ))
#(deprecated):          end=int(np.where( N[0]-y<N[1]-x, N[0]-y, N[1]-x ))
#(deprecated):          if insubap:
#(deprecated):             # terminate chain
#(deprecated):             chainsDef[0].append( newChain )
#(deprecated):          insubap=False
#(deprecated):          for j in range(end):
#(deprecated):             thisIdx=y*N[1]+x+j*(N[0]+1)
#(deprecated):             if insubap:
#(deprecated):                forceTerminate=False
#(deprecated):                if type(None)!=type(maxLen):
#(deprecated):                   if maxLen==newChain[1]:
#(deprecated):                      forceTerminate=True
#(deprecated):                if ( thisIdx in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1]+1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx%N[0]<N[1]-1 ) and not forceTerminate:
#(deprecated):                   # continue chain
#(deprecated):                   newChain[0].append( thisIdx ) ; newChain[1]+=1
#(deprecated):                elif thisIdx in self.gInst.illuminatedCornersIdx:
#(deprecated):                   # must terminate chain but include this point
#(deprecated):                   newChain[0].append( thisIdx ) ; newChain[1]+=1
#(deprecated):                   insubap=False
#(deprecated):                   chainsDef[0].append( newChain )
#(deprecated):                else:
#(deprecated):                   # terminate chain
#(deprecated):                   insubap=False
#(deprecated):                   chainsDef[0].append( newChain )
#(deprecated):             else:
#(deprecated):                if ( thisIdx in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1]+1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx%N[0]<N[1]-1 ):
#(deprecated):                   # start chain
#(deprecated):                   chainsNumber+=1
#(deprecated):                   insubap=True
#(deprecated):                   newChain=[[thisIdx],1, int(chainsNumber)]
#(deprecated):
#(deprecated):                   # \/ will always be a new entry
#(deprecated):                   if thisIdx in chainsStarts[0][1]:
#(deprecated):                      # so this should never happen
#(deprecated):                      raise ValueError("dicure.py: Gone bonkers")
#(deprecated):                   chainsStarts[0][0].append(chainsNumber)
#(deprecated):                   chainsStarts[0][1].append(thisIdx)
#(deprecated):                   chainsStarts[0][2].append(chaninsStrtNum)
#(deprecated):                   chaninsStrtNum+=1
#(deprecated):
#(deprecated):       for i in range(N[1]+N[0]-1):
#(deprecated):          # top right, downwards left
#(deprecated):          y=int(np.where( N[0]-i-1<0, 0, N[0]-i-1 )) # starting y, then x
#(deprecated):          x=int(np.where( N[0]-i-1>0, N[1]-1, N[1]-1+N[0]-i-1 ))
#(deprecated):          end=int(np.where( N[0]-y<=x, N[0]-y, x+1 ))
#(deprecated):          if insubap:
#(deprecated):             # terminate chain
#(deprecated):             chainsDef[1].append( newChain )
#(deprecated):          insubap=False
#(deprecated):          for j in range(end):
#(deprecated):             thisIdx=y*N[1]+x+j*(N[0]-1)
#(deprecated):             if insubap:
#(deprecated):                forceTerminate=False
#(deprecated):                if type(None)!=type(maxLen):
#(deprecated):                   if maxLen==newChain[1]:
#(deprecated):                      forceTerminate=True
#(deprecated):                if ( thisIdx in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1]-1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx-1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx%N[1]>0) and not forceTerminate:
#(deprecated):                   # continue chain
#(deprecated):                   newChain[0].append( thisIdx ) ; newChain[1]+=1
#(deprecated):                elif thisIdx in self.gInst.illuminatedCornersIdx:
#(deprecated):                   # must terminate chain but include this point
#(deprecated):                   newChain[0].append( thisIdx ) ; newChain[1]+=1
#(deprecated):                   insubap=False
#(deprecated):                   chainsDef[1].append( newChain )
#(deprecated):                else:
#(deprecated):                   # terminate chain
#(deprecated):                   insubap=False
#(deprecated):                   chainsDef[1].append( newChain )
#(deprecated):             else:
#(deprecated):                if ( thisIdx in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx+N[1]-1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx-1 in self.gInst.illuminatedCornersIdx and
#(deprecated):                      thisIdx%N[1]>0 ):
#(deprecated):                   # start chain
#(deprecated):                   chainsNumber+=1
#(deprecated):                   insubap=True
#(deprecated):                   newChain=[[thisIdx],1,int(chainsNumber)]
#(deprecated):
#(deprecated):                   # \/ principle code to check overlaps of chains
#(deprecated):                   if thisIdx in chainsStarts[0][1]:
#(deprecated):                      # have an overlap so direct there
#(deprecated):                      chainsStarts[1][2].append(
#(deprecated):                            chainsStarts[0][2][
#(deprecated):                               chainsStarts[0][1].index(thisIdx)
#(deprecated):                            ])
#(deprecated):                   elif thisIdx in chainsStarts[1][1]:
#(deprecated):                      # should never happen
#(deprecated):                      raise ValueError("dicure.py: Gone bonkers")
#(deprecated):                   else:
#(deprecated):                      chainsStarts[1][2].append(chaninsStrtNum)
#(deprecated):                      chaninsStrtNum+=1
#(deprecated):
#(deprecated):                   # \/ and store which chain we are talking about
#(deprecated):                   chainsStarts[1][0].append(chainsNumber)
#(deprecated):                   chainsStarts[1][1].append(thisIdx)
#(deprecated):
#(deprecated):       chainsStarts[-1]=chaninsStrtNum
#(deprecated):       return chainsDef, chainsStarts
#(deprecated):    # END of chainsDefine
#(deprecated):
#(deprecated):    def chainsIntegrate(self, rgradV ):
#(deprecated):       chains=([],[])
#(deprecated):       for i in range(len(self.chainsDef[0])):
#(deprecated):          chains[0].append([0])
#(deprecated):          for j in range(self.chainsDef[0][i][1]-1): # never use first element, defined as zero 
#(deprecated):             tC=self.chainsDef[0][i][0][j] # chain index
#(deprecated):             tgC=(self.gInst.mask2allcorners==tC).nonzero()[0]  # index into gradient vector
#(deprecated):             if type(tgC)==None:
#(deprecated):                raise ValueError("Failed to index gradient")
#(deprecated):             if len(tgC)>1:
#(deprecated):                raise ValueError("Conflict!")
#(deprecated):             chains[0][-1].append(chains[0][-1][-1]+rgradV[tgC]*2**0.5)
#(deprecated):          chains[0][-1]=np.array(chains[0][-1])
#(deprecated):       for i in range(len(self.chainsDef[1])):
#(deprecated):          chains[1].append([0])
#(deprecated):          for j in range(self.chainsDef[1][i][1]-1): # never use first element, defined as zero 
#(deprecated):             tC=self.chainsDef[1][i][0][j] # chain index
#(deprecated):                # \/ index into gradient vector
#(deprecated):             tgC=self.gInst.numberSubaps+(self.gInst.mask2allcorners==(tC-1)).nonzero()[0]
#(deprecated):             if type(tgC)==None:
#(deprecated):                raise ValueError("Failed to index gradient")
#(deprecated):             if len(tgC)>1:
#(deprecated):                raise ValueError("Conflict!")
#(deprecated):             chains[1][-1].append(float(chains[1][-1][-1]+rgradV[tgC]*2**0.5))
#(deprecated):          chains[1][-1]=np.array(chains[1][-1])
#(deprecated):       return chains
#(deprecated):
#(deprecated):    def chainsVectorize(self, chains ):
#(deprecated):       # \/ produce the chain vector and indexing into it
#(deprecated):       chainsV=[]
#(deprecated):       chainsVOffsets=[]
#(deprecated):       for dirn in (0,1):
#(deprecated):          for tchain in chains[dirn]:
#(deprecated):             chainsVOffsets.append( len(chainsV) )
#(deprecated):             chainsV+=tchain.tolist()
#(deprecated):       return chainsV, chainsVOffsets
#(deprecated):
#(deprecated):    def chainsOverlaps(self):
#(deprecated):       # Solution approach is: each chain has an offset, as does the orthogonal
#(deprecated):       # chain at the edge where they meet. The system can be pinned down by
#(deprecated):       # defining the chain values plus their offsets to be equal, or equivalently,
#(deprecated):       # their difference to be zero. Note that some orthogonal chains share
#(deprecated):       # offsets.  This latter fact permits the reconstruction. 
#(deprecated):          # \/ for each point in one set of chains, find matching in the other
#(deprecated):       chainsMatching=[]
#(deprecated):       for i in range(len(self.chainsDef[0])): # index the chains we want to align
#(deprecated):          for k in range(0, self.chainsDef[0][i][1], 1):
#(deprecated):                # /\ can choose which points along a chaing to try and match
#(deprecated):             for j in range(len(self.chainsDef[1])): # index the orthogonal chains
#(deprecated):                for l in range(0, self.chainsDef[1][j][1], 1):
#(deprecated):                   if self.chainsDef[0][i][0][k]==self.chainsDef[1][j][0][l]: # cross-over
#(deprecated):                      chainsMatching.append( (i,j,k,l) )
#(deprecated):       return chainsMatching
#(deprecated):
#(deprecated):    def chainsDefMatrices(self, chainsOvlps, chainsV, chainsVOffsets):
#(deprecated):       # If the problem is Ao+Bc=0, where o is the vector of offsets for the
#(deprecated):       # chains, which are essentially the currently undefined start value for each
#(deprecated):       # chain, and c is the vector of chain values (chainsV here), then A is the
#(deprecated):       # matrix which selects the appropriate offsets and B selects the appropriate
#(deprecated):       # chain differences which ought to be zero thus the equation.  ->
#(deprecated):       # o=A^{+}B(-c) where A^{+} is a pseudo-inverse. Or you can use the
#(deprecated):       # least-squares approach, o=(A^{T}A+{\alpha}I)^{-1}A^{T}B(-c) and standard
#(deprecated):       # matrix algebra works, as long as \alpha!=0 i.e. regularisation is applied.
#(deprecated):
#(deprecated):       # only need to consider unique offsets; some chains will start at the same
#(deprecated):       # point so that would overconstrain the matrix
#(deprecated):       A=np.zeros([ len(chainsOvlps), self.chainsDefChStrts[2] ], np.int32)
#(deprecated):       B=np.zeros([ len(chainsOvlps), len(chainsV) ], np.int32 )
#(deprecated):       for i in range(len(chainsOvlps)):
#(deprecated):          tcO=chainsOvlps[i]
#(deprecated):
#(deprecated):          B[ i, chainsVOffsets[self.chainsDef[0][tcO[0]][-1]-1]+tcO[2] ]=1
#(deprecated):          B[ i, chainsVOffsets[self.chainsDef[1][tcO[1]][-1]-1]+tcO[3] ]=-1
#(deprecated):
#(deprecated):
#(deprecated):          tcCN=[ self.chainsDef[dirn][tcO[dirn]][2] for dirn in (0,1) ]
#(deprecated):          coI=[ self.chainsDefChStrts[dirn][0].index( tcCN[dirn] ) for dirn in (0,1) ]
#(deprecated):
#(deprecated):          A[ i, self.chainsDefChStrts[0][2][coI[0]]]=1
#(deprecated):          A[ i, self.chainsDefChStrts[1][2][coI[1]]]+=-1 # may overwrite above
#(deprecated):
#(deprecated):#(deprecated):       return A,B
