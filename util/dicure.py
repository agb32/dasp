# What is this?
# Two classes for dicure (diagonal cure) wavefront reconstruction:
#
# - gradientOperatorType1 : utility class to calculate the gradient operator
#
# - DiCuRe : main class
#   - calc() : the method that performs the reconstruction
# 
# Author: NAB
# Included into aosim: UB
#

import numpy as np
import util.dot as quick

# 
# UB, 2012 Aug 6th:
# I copied (this part of) class gradientOperatorType1 from NAB's "gradientOperator.py".
# I only copied what is needed by aosim and nothing else.
# 
# Hardy, 1998, p.270 for configurations
# Type 1 = Fried
class gradientOperatorType1:
    '''Using Type 1 geometry, define a gradient operator matrix.'''
    op=None
    numberSubaps=None
    numberPhases=None

    def __init__( self, subapMask=None, pupilMask=None ):
        if subapMask!=None: self.newSubaperturesGiven(subapMask)
        if pupilMask!=None: self.newPupilGiven(pupilMask)

    def newSubaperturesGiven(self, subapMask):
        self.op=None # reset
        self.n=subapMask.shape
        self.n_=[ x+1 for x in self.n]

      # now define the illuminated sub-apertures (-1 = not illuminated)
        self.maskIdx=np.arange(self.n[0]*self.n[1]).reshape(self.n)*subapMask+(subapMask-1)
         # \/ corners on the array
        self.cornersIdx=np.arange(self.n_[0]*self.n_[1]).reshape(self.n_) 
         # \/ per corner, how many subapertures use this corner
        self.illuminatedCorners=np.zeros(self.n_)
 
      # \/ array to specify the first corner for each sub-aperture
        self.mask2allcorners=\
            np.flatnonzero(subapMask.ravel())%(self.n[1])\
            +(np.flatnonzero(subapMask.ravel())//(self.n[1]))*self.n_[1]
      # \/ figure out which phases then contribute to a sub-aperture gradient
        for i in range(self.n[0]):
            for j in range(self.n[1]):
                if self.maskIdx[i,j]==-1: continue # skip not illuminated
                self.illuminatedCorners+=(
                    (abs(self.cornersIdx//self.n_[1]-0.5-self.maskIdx[i,j]//self.n[1])<=0.5)*\
                        (abs(self.cornersIdx%self.n_[1]-0.5-self.maskIdx[i,j]%self.n[1])<=0.5) )
        self.illuminatedCornersIdx=np.flatnonzero(self.illuminatedCorners.ravel())

      # \/ some handy numbers
        self.numberSubaps=int(subapMask.sum()) # =n**2 for all illuminated
        self.numberPhases=(self.illuminatedCorners!=0).sum() # =(n+1)**2 for all illum.

    def newPupilGiven(self, pupilMask):
        '''Define the mask given a pupil mask, by creating a sub-aperture mask
        and then re-creating the pupil'''

        self.subapMask=((pupilMask[:-1,:-1]+pupilMask[1:,:-1]\
                  +pupilMask[1:,1:]+pupilMask[:-1,1:])==4)
        self.newSubaperturesGiven(self.subapMask)


#
# Class to handle the dicure reconstruction of the wavefront:
#
class DiCuRe:
    def __init__(self, subapMap):

        if subapMap.shape[0] != subapMap.shape[1]:
            raise Exception("DiCuRe can not be used with non-square subaperture maps.")

        # Set up the geometry:
        self.gInst = gradientOperatorType1( subapMap )
        # Define the chains:
        self.chainsDef,\
        self.chainsDefChStrts = self.chainsDefine() # chainsDefine returns 2 values
        # Get the number of subapertures in one direction:
        self.nsubx = subapMap.shape[0]
        # the number of actuators in one direction:
        self.nactx = self.nsubx+1

    def calc(self, gradV):
        rgradV = self.rotateVectors(gradV).ravel()
        chains = self.chainsIntegrate(rgradV)
        chainsV,chainsVOffsets = self.chainsVectorize(chains)
        chainsOvlps = self.chainsOverlaps()
        A,B = self.chainsDefMatrices(chainsOvlps, chainsV, chainsVOffsets)
        # ought to be able to remove 'chains' here, but been lazy further down
        del chainsOvlps,chainsVOffsets
        invchOScovM = np.identity(A.shape[1])
        offsetEstM = quick.dot( quick.dot( np.linalg.inv( quick.dot(A.T,A)+invchOScovM*1e-3 ), A.T), -B )
        offsetEstV = np.dot( offsetEstM, chainsV ) # 2012 Aug 03: currently you can not use quick.dot
                                                   # here because chainsV is a list and not an array
        comp = np.zeros(2*self.nactx**2, np.float64)

        # do one way...
        for x in range(len(chains[0])):
            toeI=x
            for i in range((self.chainsDef[0][x][1])):
                comp[ self.chainsDef[0][x][0][i] ] += (chains[0][x][i]+offsetEstV[toeI])
                pass
        # ...then another
        for x in range(len(chains[1])):
            toeI=self.chainsDefChStrts[1][2][x]
            for i in range((self.chainsDef[1][x][1])):
                comp[ self.nactx**2+self.chainsDef[1][x][0][i] ] += (chains[1][x][i]+offsetEstV[toeI])
                pass
        return comp

    def rotateVectors(self, grads):
        return np.array(
            [ grads[:grads.shape[0]/2]*2**-0.5+grads[grads.shape[0]/2:]*2**-0.5, 
              grads[:grads.shape[0]/2]*-2**-0.5+grads[grads.shape[0]/2:]*2**-0.5 ])

    def chainsDefine(self, maxLen=None ):
          # \/ define the chains
          # Let the 1st method be to zip through the subaperture index into the
          # pupil mask and find which sub-aperture corners are in the
          # mask2allcorners and that this element is the bottom-left (thus there
          # should also be four partners).
          # Let the 2nd method be similar but mask2allcorners+n_ (top left).
       N=self.gInst.n_ # alias
       blank=np.zeros([2]+N)
       chainsNumber=0
       chaninsStrtNum=0
       chainsDef=[],[]
       chainsStarts=[ ([],[],[]), ([],[],[]), 0]
       insubap=False
       for i in range(N[1]+N[0]-1):
          # top left, upwards right
          y=int(np.where( N[0]-i-1<0, 0, N[0]-i-1 )) # starting y, then x
          x=int(np.where( i-N[0]+1<0, 0, i-N[0]+1 ))
          end=int(np.where( N[0]-y<N[1]-x, N[0]-y, N[1]-x ))
          if insubap:
             # terminate chain
             chainsDef[0].append( newChain )
          insubap=False
          for j in range(end):
             thisIdx=y*N[1]+x+j*(N[0]+1)
             if insubap:
                forceTerminate=False
                if type(None)!=type(maxLen):
                   if maxLen==newChain[1]:
                      forceTerminate=True
                if ( thisIdx in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1]+1 in self.gInst.illuminatedCornersIdx and
                      thisIdx+1 in self.gInst.illuminatedCornersIdx and
                      thisIdx%N[0]<N[1]-1 ) and not forceTerminate:
                   # continue chain
                   newChain[0].append( thisIdx ) ; newChain[1]+=1
                elif thisIdx in self.gInst.illuminatedCornersIdx:
                   # must terminate chain but include this point
                   newChain[0].append( thisIdx ) ; newChain[1]+=1
                   insubap=False
                   chainsDef[0].append( newChain )
                else:
                   # terminate chain
                   insubap=False
                   chainsDef[0].append( newChain )
             else:
                if ( thisIdx in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1]+1 in self.gInst.illuminatedCornersIdx and
                      thisIdx+1 in self.gInst.illuminatedCornersIdx and
                      thisIdx%N[0]<N[1]-1 ):
                   # start chain
                   chainsNumber+=1
                   insubap=True
                   newChain=[[thisIdx],1, int(chainsNumber)]

                   # \/ will always be a new entry
                   if thisIdx in chainsStarts[0][1]:
                      # so this should never happen
                      raise ValueError("dicure.py: Gone bonkers")
                   chainsStarts[0][0].append(chainsNumber)
                   chainsStarts[0][1].append(thisIdx)
                   chainsStarts[0][2].append(chaninsStrtNum)
                   chaninsStrtNum+=1

       for i in range(N[1]+N[0]-1):
          # top right, downwards left
          y=int(np.where( N[0]-i-1<0, 0, N[0]-i-1 )) # starting y, then x
          x=int(np.where( N[0]-i-1>0, N[1]-1, N[1]-1+N[0]-i-1 ))
          end=int(np.where( N[0]-y<=x, N[0]-y, x+1 ))
          if insubap:
             # terminate chain
             chainsDef[1].append( newChain )
          insubap=False
          for j in range(end):
             thisIdx=y*N[1]+x+j*(N[0]-1)
             if insubap:
                forceTerminate=False
                if type(None)!=type(maxLen):
                   if maxLen==newChain[1]:
                      forceTerminate=True
                if ( thisIdx in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1]-1 in self.gInst.illuminatedCornersIdx and
                      thisIdx-1 in self.gInst.illuminatedCornersIdx and
                      thisIdx%N[1]>0) and not forceTerminate:
                   # continue chain
                   newChain[0].append( thisIdx ) ; newChain[1]+=1
                elif thisIdx in self.gInst.illuminatedCornersIdx:
                   # must terminate chain but include this point
                   newChain[0].append( thisIdx ) ; newChain[1]+=1
                   insubap=False
                   chainsDef[1].append( newChain )
                else:
                   # terminate chain
                   insubap=False
                   chainsDef[1].append( newChain )
             else:
                if ( thisIdx in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1] in self.gInst.illuminatedCornersIdx and
                      thisIdx+N[1]-1 in self.gInst.illuminatedCornersIdx and
                      thisIdx-1 in self.gInst.illuminatedCornersIdx and
                      thisIdx%N[1]>0 ):
                   # start chain
                   chainsNumber+=1
                   insubap=True
                   newChain=[[thisIdx],1,int(chainsNumber)]

                   # \/ principle code to check overlaps of chains
                   if thisIdx in chainsStarts[0][1]:
                      # have an overlap so direct there
                      chainsStarts[1][2].append(
                            chainsStarts[0][2][
                               chainsStarts[0][1].index(thisIdx)
                            ])
                   elif thisIdx in chainsStarts[1][1]:
                      # should never happen
                      raise ValueError("dicure.py: Gone bonkers")
                   else:
                      chainsStarts[1][2].append(chaninsStrtNum)
                      chaninsStrtNum+=1

                   # \/ and store which chain we are talking about
                   chainsStarts[1][0].append(chainsNumber)
                   chainsStarts[1][1].append(thisIdx)

       chainsStarts[-1]=chaninsStrtNum
       return chainsDef, chainsStarts
    # END of chainsDefine

    def chainsIntegrate(self, rgradV ):
       chains=([],[])
       for i in range(len(self.chainsDef[0])):
          chains[0].append([0])
          for j in range(self.chainsDef[0][i][1]-1): # never use first element, defined as zero 
             tC=self.chainsDef[0][i][0][j] # chain index
             tgC=(self.gInst.mask2allcorners==tC).nonzero()[0]  # index into gradient vector
             if type(tgC)==None:
                raise ValueError("Failed to index gradient")
             if len(tgC)>1:
                raise ValueError("Conflict!")
             chains[0][-1].append(chains[0][-1][-1]+rgradV[tgC]*2**0.5)
          chains[0][-1]=np.array(chains[0][-1])
       for i in range(len(self.chainsDef[1])):
          chains[1].append([0])
          for j in range(self.chainsDef[1][i][1]-1): # never use first element, defined as zero 
             tC=self.chainsDef[1][i][0][j] # chain index
                # \/ index into gradient vector
             tgC=self.gInst.numberSubaps+(self.gInst.mask2allcorners==(tC-1)).nonzero()[0]
             if type(tgC)==None:
                raise ValueError("Failed to index gradient")
             if len(tgC)>1:
                raise ValueError("Conflict!")
             chains[1][-1].append(float(chains[1][-1][-1]+rgradV[tgC]*2**0.5))
          chains[1][-1]=np.array(chains[1][-1])
       return chains

    def chainsVectorize(self, chains ):
       # \/ produce the chain vector and indexing into it
       chainsV=[]
       chainsVOffsets=[]
       for dirn in (0,1):
          for tchain in chains[dirn]:
             chainsVOffsets.append( len(chainsV) )
             chainsV+=tchain.tolist()
       return chainsV, chainsVOffsets

    def chainsOverlaps(self):
       # Solution approach is: each chain has an offset, as does the orthogonal
       # chain at the edge where they meet. The system can be pinned down by
       # defining the chain values plus their offsets to be equal, or equivalently,
       # their difference to be zero. Note that some orthogonal chains share
       # offsets.  This latter fact permits the reconstruction. 
          # \/ for each point in one set of chains, find matching in the other
       chainsMatching=[]
       for i in range(len(self.chainsDef[0])): # index the chains we want to align
          for k in range(0, self.chainsDef[0][i][1], 1):
                # /\ can choose which points along a chaing to try and match
             for j in range(len(self.chainsDef[1])): # index the orthogonal chains
                for l in range(0, self.chainsDef[1][j][1], 1):
                   if self.chainsDef[0][i][0][k]==self.chainsDef[1][j][0][l]: # cross-over
                      chainsMatching.append( (i,j,k,l) )
       return chainsMatching

    def chainsDefMatrices(self, chainsOvlps, chainsV, chainsVOffsets):
       # If the problem is Ao+Bc=0, where o is the vector of offsets for the
       # chains, which are essentially the currently undefined start value for each
       # chain, and c is the vector of chain values (chainsV here), then A is the
       # matrix which selects the appropriate offsets and B selects the appropriate
       # chain differences which ought to be zero thus the equation.  ->
       # o=A^{+}B(-c) where A^{+} is a pseudo-inverse. Or you can use the
       # least-squares approach, o=(A^{T}A+{\alpha}I)^{-1}A^{T}B(-c) and standard
       # matrix algebra works, as long as \alpha!=0 i.e. regularisation is applied.

       # only need to consider unique offsets; some chains will start at the same
       # point so that would overconstrain the matrix
       A=np.zeros([ len(chainsOvlps), self.chainsDefChStrts[2] ], np.int32)
       B=np.zeros([ len(chainsOvlps), len(chainsV) ], np.int32 )
       for i in range(len(chainsOvlps)):
          tcO=chainsOvlps[i]

          B[ i, chainsVOffsets[self.chainsDef[0][tcO[0]][-1]-1]+tcO[2] ]=1
          B[ i, chainsVOffsets[self.chainsDef[1][tcO[1]][-1]-1]+tcO[3] ]=-1


          tcCN=[ self.chainsDef[dirn][tcO[dirn]][2] for dirn in (0,1) ]
          coI=[ self.chainsDefChStrts[dirn][0].index( tcCN[dirn] ) for dirn in (0,1) ]

          A[ i, self.chainsDefChStrts[0][2][coI[0]]]=1
          A[ i, self.chainsDefChStrts[1][2][coI[1]]]+=-1 # may overwrite above

       return A,B
