# What is this?
# Generate gradient operator based on pupil mask
#
# 2012/Sept/22 - renamed maskIdx as subapMaskSequential

#For ali:
#g=gradientOperator.gradientOperatorType1(pupilMask=actmask,sparse=1)
#gradientOperator.genericLaplacianCalcOp_NumpyArray(g)
#plot(g.op)
#Then: numpy.dot(g.op,g.op) for each layer...

#import matplotlib.pyplot as pg
import numpy

# Hardy, 1998, p.270 for configurations
# Type 1 = Fried
#      2 = Hudgin
#      3 = centred

class gradientOperatorType1:
   '''Using Type 1 geometry, define a gradient operator matrix.'''
   op=None
   numberSubaps=None
   numberPhases=None
   sparse=False

   def __init__( self, subapMask=None, pupilMask=None, sparse=False ):
      self.sparse=sparse
      if self.sparse:
         self.calcOp=self.calcOp_scipyCSR
      else:
         self.calcOp=self.calcOp_NumpyArray
      if subapMask!=None: self.newSubaperturesGiven(subapMask)
      if pupilMask!=None: self.newPupilGiven(pupilMask)

   def newSubaperturesGiven(self, subapMask):
      self.numberSubaps=int(subapMask.sum()) # =n**2 for all illuminated
      if self.numberSubaps>1000 and not self.sparse:
         raise ValueError("Insist on sparsity for numberSubaps>1000")
      self.op=None # reset
      self.n=subapMask.shape
      self.n_=[ x+1 for x in self.n]

      # now define the illuminated sub-apertures (-1 = not illuminated)
      self.subapMaskSequential=\
       numpy.arange(self.n[0]*self.n[1]).reshape(self.n)*subapMask+(subapMask-1)
         # \/ corners on the array
      self.cornersIdx=numpy.arange(self.n_[0]*self.n_[1]).reshape(self.n_) 
         # \/ per corner, how many subapertures use this corner
      self.illuminatedCorners=numpy.zeros(self.n_)
 
      # \/ array to specify the first corner for each sub-aperture
      self.mask2allcorners=\
         numpy.flatnonzero(subapMask.ravel())%(self.n[1])\
        +(numpy.flatnonzero(subapMask.ravel())//(self.n[1]))*self.n_[1]
      # \/ figure out which phases then contribute to a sub-aperture gradient
      for i in range(self.n[0]):
        for j in range(self.n[1]):
          if self.subapMaskSequential[i,j]==-1: continue # skip not illuminated
          self.illuminatedCorners+=(
            (abs(self.cornersIdx//self.n_[1]-0.5
                  -self.subapMaskSequential[i,j]//self.n[1])<=0.5)*
            (abs(self.cornersIdx%self.n_[1]
                  -0.5-self.subapMaskSequential[i,j]%self.n[1])<=0.5) )
      self.illuminatedCornersIdx=numpy.flatnonzero(
            self.illuminatedCorners.ravel())
      self.subapMaskIdx=numpy.flatnonzero( subapMask.ravel() )

      self.numberPhases=(self.illuminatedCorners!=0).sum() # =(n+1)**2 for all illum.
      self.subapMask=subapMask
      print self.numberSubaps,self.numberPhases

   def newPupilGiven(self, pupilMask):
      '''Define the mask given a pupil mask, by creating a sub-aperture mask
      and then re-creating the pupil'''
      self.newSubaperturesGiven( ((pupilMask[:-1,:-1]+pupilMask[1:,:-1]\
                  +pupilMask[1:,1:]+pupilMask[:-1,1:])==4) )

   def testInit(self, subapMask):
      import matplotlib.pyplot as pg
      # Test
      pg.imshow(subapMask,interpolation='nearest',
         extent=[0,self.n[1],0,self.n[0]])
      pg.plot( numpy.flatnonzero(self.illuminatedCorners.ravel())%self.n_[1],
               numpy.flatnonzero(self.illuminatedCorners.ravel())//self.n_[0],
                  'wo' )
      pg.title("red=illuminated, white=their corners")
      pg.show()

   def returnOp(self):
      if self.numberSubaps==None: return None
      self.calcOp()
      return self.op

   cornerIdxs=lambda self, i : (
      self.mask2allcorners[i], self.mask2allcorners[i]+1,
      self.mask2allcorners[i]+self.n_[1], self.mask2allcorners[i]+self.n_[1]+1 )
      # \/ want for each sub-aperture, the corresponding corners
   r=lambda self,i : [
      numpy.flatnonzero(self.illuminatedCornersIdx==tm)[0] for tm in
      self.cornerIdxs(i) ]
   def calcOp_NumpyArray(self):
      # gradient matrix
      self.op=numpy.zeros(
         [2*self.numberSubaps,self.numberPhases],numpy.float64)
      for i in range(self.numberSubaps):
         r=self.r(i)
            # \/ grad x
         self.op[i,r[0]]=-0.5
         self.op[i,r[1]]=0.5
         self.op[i,r[2]]=-0.5
         self.op[i,r[3]]=0.5
            # \/ self.grad y
         self.op[i+self.numberSubaps,r[0]]=-0.5
         self.op[i+self.numberSubaps,r[1]]=-0.5
         self.op[i+self.numberSubaps,r[2]]=0.5
         self.op[i+self.numberSubaps,r[3]]=0.5

   def calcOp_scipyCSR(self):
      import scipy.sparse # only required if we reach this stage
      row,col=[],[] ; data=[]
      for i in range(self.numberSubaps):
         r=self.r(i)
         for j in range(4): # x,y pairs
            row.append(i)
            col.append(r[j])
            data.append( ((j%2)*2-1)*0.5 )
            #
            row.append(i+self.numberSubaps)
            col.append(r[j])
            data.append( (1-2*(j<2))*0.5 )
         self.op=scipy.sparse.csr_matrix((data,(row,col)), dtype=numpy.float64)

# ------------------------------

class gradientOperatorType2Fried(gradientOperatorType1):
   '''Using Type 2 geometry to define a gradient operator matrix.
    Can rotate the gradient basis to get a Fried-like geometry so this
    is a straightforward approach which can promote waffle.'''

   def calcOp_NumpyArray(self):
      # Gradient matrix
      # Similar to Fried but rotation of phases
      self.op=numpy.zeros(
         [2*self.numberSubaps,self.numberPhases],numpy.float64)
      for i in range(self.numberSubaps):
         # want for each sub-aperture, the corresponding corners
         r=[ numpy.flatnonzero(self.illuminatedCornersIdx==tm)[0] for tm in
            (self.mask2allcorners[i],
             self.mask2allcorners[i]+1,
             self.mask2allcorners[i]+self.n_[1],
             self.mask2allcorners[i]+self.n_[1]+1) ]
            # \/ grad x1->y2
         self.op[i,r[0]]=-1
         self.op[i,r[3]]=1
            # \/ self.grad x2->y1
         self.op[i+self.numberSubaps,r[1]]=1
         self.op[i+self.numberSubaps,r[2]]=-1

class gradientOperatorType2(gradientOperatorType1):
   '''Using Type 2 geometry to define a gradient operator matrix.
    Define explicitly.
    The subaperture mask is now not equal for the x and y directions,
    so the interpretation is that it is extended by one in the x and y
    directions for the respective gradients (meaning they aren't identical).
    The class must therefore provide a method for extending the subaperture
    mask from that given which is some 'mean' value.'''
   subapYMask=None
   subapXMask=None
   nx=None ; ny=None
   numberXSubaps=None
   numberYSubaps=None

   def extendSubapMask(self, subapMask):
      # in the Hudgin geometry, there are more measurements for x in the y
      # direction and vice-versa.
      self.n=subapMask.shape
      
      self.subapXMask=numpy.zeros([self.n[1]+1,self.n[0]],numpy.bool)
      self.subapXMask[:-1]=subapMask
      self.subapXMask[1:]+=subapMask 
      
      self.subapYMask=numpy.zeros([self.n[1],self.n[0]+1],numpy.bool)
      self.subapYMask[:,:-1]=subapMask
      self.subapYMask[:,1:]+=subapMask 

      self.nx=self.subapXMask.shape
      self.ny=self.subapYMask.shape


   def newSubaperturesGiven(self, subapMask):
      self.op=None # reset
      self.extendSubapMask(subapMask)
      self.n=subapMask.shape
      self.n_=[ x+1 for x in self.n]

      # now define the illuminated sub-apertures (-1 = not illuminated)
      self.maskYIdx= numpy.arange(self.ny[0]*self.ny[1]).reshape(self.ny)\
                     *self.subapYMask+(self.subapYMask-1)
      self.maskXIdx= numpy.arange(self.nx[0]*self.nx[1]).reshape(self.nx)\
                     *self.subapXMask+(self.subapXMask-1)
         # \/ corners on the array
      self.cornersIdx=numpy.arange(self.n_[0]*self.n_[1]).reshape(self.n_) 
         # \/ per corner, how many subapertures use this corner
      self.illuminatedCorners=numpy.zeros(self.n_,numpy.complex64)
 
      # \/ array to specify the first corner for each sub-aperture
      self.maskX2allcorners=\
         numpy.flatnonzero(self.subapXMask.ravel())%(self.nx[1])\
        +(numpy.flatnonzero(self.subapXMask.ravel())//(self.nx[0]-1))*self.n_[0]
      self.maskY2allcorners=\
         numpy.flatnonzero(self.subapYMask.ravel())%(self.ny[1])\
        +(numpy.flatnonzero(self.subapYMask.ravel())//(self.ny[0]+1))*self.n_[0]
      # \/ figure out which phases then contribute to a sub-aperture gradient
      for i in range(self.nx[0]):
        for j in range(self.nx[1]):
          maskVal=self.maskXIdx[i,j]
          if maskVal==-1: continue # skip not illuminated
          self.illuminatedCorners+=(
            (self.cornersIdx//self.n_[1]==maskVal//self.nx[1])*\
            (abs(self.cornersIdx%self.n_[0]-0.5-maskVal%(self.nx[0]-1))<=0.5) )
      for i in range(self.ny[0]):
        for j in range(self.ny[1]):
          maskVal=self.maskYIdx[i,j]
          if maskVal==-1: continue # skip not illuminated
          self.illuminatedCorners+=1.0j*(
            (self.cornersIdx%self.n_[0]==maskVal%(self.ny[0]+1))*\
            (abs(self.cornersIdx//self.n_[1]-0.5-maskVal//self.ny[1])<=0.5) )
         # \/ have to compare illumination via this slightly strage method
      if (self.illuminatedCorners.real>0).sum() !=\
         (self.illuminatedCorners.imag>0).sum():
         raise ValueError("Failure of commonality between phase points")
      self.illuminatedCornersIdx=\
         numpy.flatnonzero(self.illuminatedCorners.ravel())

      # \/ some handy numbers
      self.numberXSubaps=int(self.subapXMask.sum())
      self.numberYSubaps=int(self.subapYMask.sum()) 
      self.numberPhases=(self.illuminatedCorners!=0).sum() # =(n+1)**2 for all illum.


   def returnOp(self):
      if self.numberXSubaps==None: return None
      if self.op!=None: return self.op
      self.calcOp()
      return self.op

   def calcOp_NumpyArray(self):
      self.op=numpy.zeros( [self.numberXSubaps+self.numberYSubaps,
         self.numberPhases],numpy.float64)
      for i in range(self.numberXSubaps):
         # want for each sub-aperture, the corresponding corners
         r=[ numpy.flatnonzero(self.illuminatedCornersIdx==tm)[0] for tm in
            (self.maskX2allcorners[i], self.maskX2allcorners[i]+1) ]
         self.op[i,r[0]]=-1
         self.op[i,r[1]]=1

      for i in range(self.numberYSubaps):
         # want for each sub-aperture, the corresponding corners
         r=[ numpy.flatnonzero(self.illuminatedCornersIdx==tm)[0] for tm in
            (self.maskY2allcorners[i], self.maskY2allcorners[i]+self.n_[1]) ]
         self.op[i+self.numberXSubaps,r[0]]=-1
         self.op[i+self.numberXSubaps,r[1]]=1

# ------------------------------


class gradientOperatorType3Avg(gradientOperatorType1):
   '''Using Type 3 geometry and averaging gradients, define a gradient operator
     matrix.'''
   # Convert the Southwell geometry into a Hudgin one via an averaging matrix,
   # so overload the Hudgin operator.
   def __init__(self, subapMask):
      raise NotimplementedError("Type 3 with averaging not yet implemented")


class gradientOperatorType3Centred(gradientOperatorType1):
   '''Using Type 3 geometry and fake sub-apertures, define a gradient operator
     matrix.'''
   # Change from type 1 is that there isn't enough information for corners so
   # fake phases must be added (these can be arbirtrary) but the responsibility
   # for this lies with the constructed phase vector This makes this code far
   # simpler but we ought to also supply a member that suitably constricts and
   # expands a phase vector to the 'right' length.

   def newSubaperturesGiven(self, subapMask):
      self.op=None # reset
      self.n=subapMask.shape
      self.n_=[ x+2 for x in self.n] # pad by an extra row and column

      # now define the illuminated sub-apertures (-1 = not illuminated)
      self.subapMaskSequential=numpy.arange(
            self.n[0]*self.n[1]).reshape(self.n)*subapMask+(subapMask-1)
         # \/ per phase centre, how many subapertures use this corner
      self.allCentres=numpy.zeros(self.n_)
      self.illuminatedCentres=numpy.zeros(self.n_)
      self.illuminatedCentres[1:-1,1:-1]=subapMask
         # \/ zip around the mask and add 'fake' subaperture where
         #  reqd.
      for i in range(self.n[0]):
         for j in range(self.n[1]):
            if subapMask[i,j]:
               for k in (0,2):
                  self.allCentres[i+k,j+1]+=1
                  self.allCentres[i+1,j+k]+=1

      self.illuminatedCentresIdx=\
         numpy.flatnonzero(self.illuminatedCentres.ravel())
      self.allCentresIdx=numpy.flatnonzero(self.allCentres.ravel())

      # \/ now figure out the real and fake subaperture indices
      # now define the illuminated sub-apertures (-1 = not illuminated,0
      # fake)
      self.subapMaskSequential=-numpy.ones(self.n_)
      self.subapMaskSequential.ravel()[self.allCentresIdx]=0
      self.subapMaskSequential.ravel()[self.illuminatedCentresIdx]=1
         # \/ index into the phase vectors, to insert/extract real values
         #  into/from fake values vector
      self.fake2real=numpy.flatnonzero(
         self.subapMaskSequential.ravel()[self.allCentresIdx]==1)
      self.real2fake=numpy.flatnonzero(self.subapMaskSequential.ravel()==1)

      # some handy numbers
      self.numberSubaps=subapMask.sum() # =n**2 for all valid illuminated
      self.numberPhases=(self.allCentres!=0).sum() 
         # =(n+2)**2 for all illum., this is of course the fake version

   def testInit(self, subapMask):
      import matplotlib.pyplot as pg
      # Test
      pg.imshow(subapMask,interpolation='nearest',
         extent=[0.5,self.n[1]+0.5,0.5,self.n[0]+0.5])
      pg.axis([0,self.n_[1],0,self.n_[0]])
      pg.plot( numpy.flatnonzero(self.allCentres.ravel())%self.n_[1],
               numpy.flatnonzero(self.allCentres.ravel())//self.n_[0],
                  'wo' )
      pg.title("red=illuminated, white=their corners")
      pg.waitforbuttonpress()
   
   def constrictPhase(self,phaseV):
      '''Given a fake phase vector, extract the real values'''
      return numpy.take(phaseV,self.fake2real,axis=-1)
   def expandPhase(self,phaseV):
      '''Given a real phase vector, create a fake one'''
      newPhaseV=numpy.zeros(self.numberPhases)
      newPhaseV[self.fake2real]=phaseV
      return newPhaseV
   
   def calcOp_NumpyArray(self):
      # gradient matrix
      self.op=numpy.zeros(
         [2*self.numberSubaps,self.numberPhases],numpy.float64)
      for i in range(self.numberSubaps):
         # want for each /real/ sub-aperture, the centred difference
         r=[ numpy.flatnonzero(self.allCentresIdx==tm)[0] for tm in
            (self.illuminatedCentresIdx[i]-1,
             self.illuminatedCentresIdx[i]+1,
             self.illuminatedCentresIdx[i]-self.n_[1],
             self.illuminatedCentresIdx[i]+self.n_[1]) ]
            # \/ grad x
         self.op[i,r[0]]=-0.5
         self.op[i,r[1]]=0.5
            # \/ self.grad y
         self.op[i+self.numberSubaps,r[2]]=-0.5
         self.op[i+self.numberSubaps,r[3]]=0.5

# ------------------------------


class waffleOperatorType1(gradientOperatorType1):
   '''Define an operator that extracts waffle'''
   def calcOp_NumpyArray(self):
      self.op=((self.illuminatedCornersIdx//self.n_[1]%2)
              +(self.illuminatedCornersIdx%self.n_[1])%2)%2*2-1.0
      self.op-=numpy.mean(self.op) # remove the mean component
      self.op*=numpy.dot(self.op,self.op)**-0.5 # rescale


   # biharmonic stencil is
   # [0  1  0]
   # [1 -4  1]
   # [0  1  0]
   # so for each point, want the four nearest neighbours,
laplacianStencil=[1,1,-4,1,1]

# The approximation to the kolmogorov inverse stencil, based on explicit
# SVD-based inversion is, with the centre weighted at +4 c.f. the laplacian,
#
#([[ -1.8e-03,   9.1e-02,   2.5e-01,   9.1e-02,  -2.9e-03],
#  [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
#  [  2.5e-01,  -1.3e+00,   4.0e+00,  -1.3e+00,   2.5e-01],
#  [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
#  [ -2.9e-03,   8.9e-02,   2.5e-01,   8.9e-02,  -4.2e-03]]) 

kolmogInverseStencil=numpy.array(
      [[ -1.8e-03,   9.1e-02,   2.5e-01,   9.1e-02,  -2.9e-03],
       [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
       [  2.5e-01,  -1.3e+00,   4.0e+00,  -1.3e+00,   2.5e-01],
       [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
       [ -2.9e-03,   8.9e-02,   2.5e-01,   8.9e-02,  -4.2e-03]]).ravel()

def genericKolmogInverse_findLocation(self, i):
   thisIdx=self.illuminatedCornersIdx[i]
   wantedIdx=thisIdx+\
       numpy.add.outer(numpy.arange(-2,3)*self.n_[1],numpy.arange(-2,3)).ravel()
      # \/ check if found indices are in the right place
   rowcolIdx=[ numpy.add.outer(numpy.arange(-2,3)+thisIdx//self.n_[1],
                  numpy.zeros(5)).ravel(),
               numpy.add.outer(numpy.zeros(5),
                  numpy.arange(-2,3)+thisIdx%self.n_[1]).ravel() ]
   valid=[ numpy.flatnonzero(
         (self.illuminatedCornersIdx==wantedIdx[i])*
         ((wantedIdx[i]//self.n_[1])==rowcolIdx[0][i])*
         ((wantedIdx[i]%self.n_[1])==rowcolIdx[1][i])
         )
      for i in range(25) ]
   return valid

def genericLaplacian_findLocation(self, i):
   thisIdx=self.illuminatedCornersIdx[i]
   wantedIdx=[
      ( thisIdx-self.n_[1],(-1,0) ),
      ( thisIdx-1,(0,-1) ),
      ( thisIdx,(0,0) ),
      ( thisIdx+1,(0,1) ),
      ( thisIdx+self.n_[1],(1,0) ) ]
   valid=[ numpy.flatnonzero(
         (self.illuminatedCornersIdx==twi[0])*
         ((twi[0]//self.n_[1]-thisIdx//self.n_[1])==twi[1][0])*
         ((twi[0]%self.n_[1]-thisIdx%self.n_[1])==twi[1][1])
         )
      for twi in wantedIdx ]
   return valid

def genericKolmogInverseCalcOp_NumpyArray(operator):
   self=operator
   self.op=numpy.zeros( [self.numberPhases]*2, numpy.float64 )
   for i in range(self.numberPhases):
      valid=genericKolmogInverse_findLocation(self,i)
      for j in range(25):
         if valid[j].shape[0]==1:
            loc=valid[j][0]
            self.op[i, loc ]=kolmogInverseStencil[j]

def genericLaplacianCalcOp_NumpyArray(operator):
   self=operator
   self.op=numpy.zeros( [self.numberPhases]*2, numpy.float64 )
   # numpy.flatnonzero(gO.illuminatedCornersIdx==3)
   for i in range(self.numberPhases):
      valid=genericLaplacian_findLocation(self,i)
      for j in range(len(laplacianStencil)):
         if valid[j].shape[0]==1:
            loc=valid[j][0]
            self.op[i, loc ]=laplacianStencil[j]

def genericLaplacianCalcOp_scipyCSR(operator):
   import scipy.sparse
   self=operator
   data=[] ; row=[] ; column=[]
   for i in range(self.numberPhases):
      valid=genericLaplacian_findLocation(self,i)
      for j in range(len(laplacianStencil)):
         if valid[j].shape[0]==1:
            row.append(i) ; column.append(valid[j][0])
            data.append(laplacianStencil[j])
   self.op=scipy.sparse.csr_matrix( (data, (row,col)), dtype=numpy.float64 )

#class laplacianOperatorPupil:
#   '''Define an operator that computes the Laplacian, based on a pupil mask'''
#   def __init__(self, pupilMask):
#         # Create relevant data from the given pupilMask
#      self.numberPhases=(pupilMask!=0).sum()
#      self.illuminatedCornersIdx=numpy.flatnonzero( pupilMask.ravel() )
#      self.n_=pupilMask.shape
#      genericLaplacianCalcOp_NumpyArray(self)

class laplacianOperatorType1(gradientOperatorType1):
   '''Define an operator that computes the Laplacian'''
   def calcOp_NumpyArray(self):
      genericLaplacianCalcOp_NumpyArray(self)
   def calcOp_scipyCSR(self):
      genericLaplacianCalcOp_scipyCSR(self)

class laplacianOperatorType2(gradientOperatorType2):
   '''Define an operator that computes the Laplacian, for Type 2'''
   def calcOp_NumpyArray(self):
      genericLaplacianCalcOp_NumpyArray(self)

class laplacianOperatorType3Centred(gradientOperatorType3Centred):
   '''Define an operator that computes the Laplacian, for Type 3'''
   def calcOp_NumpyArray(self):
      genericLaplacianCalcOp_NumpyArray(self)

class kolmogInverseOperatorType1(gradientOperatorType1):
   '''Define an operator that computes an approximation to the -11/3
      power spectrum'''
   def calcOp_NumpyArray(self):
      genericKolmogInverseCalcOp_NumpyArray(self)

# ------------------------------


class localMaskOperatorType1(gradientOperatorType1):
   '''Define an operator that incorporates only certain distances within an
   operator matrix.'''
   def __init__( self, subapMask=None, distances=[0] ):
      self.distances=distances
      gradientOperatorType1.__init__(self,subapMask)

   def calcOp_NumpyArray(self):
      self.op=numpy.zeros( [self.numberPhases]*2, numpy.float64 )
      # distance stencil is,
      #  for every distance,
      #   if dist not 0:
      #      wantedIdx+=[
      # thisIdx-dist*self.n_[1],thisIdx-dist,thisIdx+dist,thisIdx+dist*self.n_[1] ]
      #   else:
      #      wantedIdx+=[thisIdx]
      #  valid=[ twi in self.illuminatedCornersIdx for twi in wantedIdx ]
      for i in range(self.numberPhases):
         thisIdx=self.illuminatedCornersIdx[i]
         self.finalWantedIdx=[]
         for td in self.distances:
            if td!=0:
               wantedIdx=[ thisIdx-td*self.n_[1],thisIdx-td,
                  thisIdx+td,thisIdx+td*self.n_[1] ]
               valid=numpy.array([
                  twi in self.illuminatedCornersIdx for twi in wantedIdx ]
                     ).nonzero()[0]
               wantedIdx=numpy.take(wantedIdx,valid)
            else:
               wantedIdx=[ thisIdx ]
            self.finalWantedIdx+=numpy.searchsorted(
               self.illuminatedCornersIdx, wantedIdx ).tolist()
         for j in self.finalWantedIdx:
            self.op[i, j]+=1


if __name__=="__main__":
   nfft=7 # pixel size
   onePhase=numpy.random.normal(size=[nfft+2]*2 ) # Type 1 & Type 3
   # define pupil mask as sub-apertures
   pupilMask=numpy.ones([nfft]*2,numpy.int32)
   pupilCds=numpy.add.outer(
      (numpy.arange(nfft)-(nfft-1)/2.)**2.0,
      (numpy.arange(nfft)-(nfft-1)/2.)**2.0 )
   pupilMask=(pupilCds>2**2)*(pupilCds<4**2)

   # test, 06/Aug/2012 to confirm UB's issue for CURE
   pupilMask=numpy.array([[0,0,1,1,1,0,0], [0,1,1,1,1,1,0], [1,1,1,1,1,1,1],
            [1,1,1,0,1,1,1], [1,1,1,1,1,1,1], [0,1,1,1,1,1,0], [0,0,1,1,1,0,0]])

   def checkType1():
      # checking code, Type 1
      gO=gradientOperatorType1(pupilMask)
     
         # \/ constrict the phase screen by averaging (makes it
         #  consistent with Type 3
      thisPhase=\
         onePhase[1:,1:]+onePhase[:-1,1:]+onePhase[:-1,:-1]+onePhase[1:,:-1]
      onePhaseV=thisPhase.ravel()[gO.illuminatedCornersIdx] # phases

      gM=gO.returnOp()
      gradV=numpy.dot(gM,onePhaseV)
      
      displayGradsD=numpy.zeros([2*nfft**2],numpy.float64)
      displayGradsD[
       numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1).tolist()+
       (nfft**2+numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1)).tolist()]\
        =gradV
      displayGradsD.resize([2,nfft,nfft])
      displayGradsD=\
         numpy.ma.masked_array( displayGradsD, [pupilMask==0,pupilMask==0] )

      # theoretical gradients 
      calcdGradsD=numpy.ma.masked_array(numpy.array([
        0.5*(thisPhase[1:,1:]-thisPhase[1:,:-1])
         +0.5*(thisPhase[:-1,1:]-thisPhase[:-1,:-1]),
        0.5*(thisPhase[1:,1:]+thisPhase[1:,:-1])
         -0.5*(thisPhase[:-1,1:]+thisPhase[:-1,:-1]) ]
        ), [pupilMask==0,pupilMask==0] )

      return calcdGradsD, displayGradsD
   
   def checkType2():
      # checking code, Type 2
      gO=gradientOperatorType2(pupilMask)
     
         # \/ constrict the phase screen by averaging (makes it
         #  consistent with Type 3
      thisPhase=\
         onePhase[1:,1:]+onePhase[:-1,1:]+onePhase[:-1,:-1]+onePhase[1:,:-1]
      onePhaseV=thisPhase.ravel()[gO.illuminatedCornersIdx] # phases

      gM=gO.returnOp()
      gradV=numpy.dot(gM,onePhaseV)
      
      displayGradsD=[ numpy.zeros([(nfft+1)*nfft],numpy.float64),
                      numpy.zeros([(nfft+1)*nfft],numpy.float64) ]
      displayGradsD[0][numpy.flatnonzero(gO.maskXIdx.ravel()>-1)]=\
         gradV[:gO.numberXSubaps]
      displayGradsD[1][numpy.flatnonzero(gO.maskYIdx.ravel()>-1)]=\
         gradV[gO.numberXSubaps:]
      displayGradsD[0].resize([nfft+1,nfft])
      displayGradsD[1].resize([nfft,nfft+1])
      displayGradsD[0]=\
         numpy.ma.masked_array( displayGradsD[0], gO.subapXMask==0 )
      displayGradsD[1]=\
         numpy.ma.masked_array( displayGradsD[1], gO.subapYMask==0 )

      # theoretical gradients 
      calcdGradsD=[
         numpy.ma.masked_array(thisPhase[:,1:]-thisPhase[:,:-1],
            gO.subapXMask==0 ),\
         numpy.ma.masked_array(thisPhase[1:]-thisPhase[:-1],
            gO.subapYMask==0 )
            ]
      return calcdGradsD, displayGradsD



   def checkType3C():
      # checking code, Type 3, centred
      gO=gradientOperatorType3Centred(pupilMask)
     
      thisPhase=onePhase # the edge values are never used
      onePhaseV=thisPhase.ravel()[gO.allCentresIdx]
      
      gM=gO.returnOp()
      gradV=numpy.dot(gM,onePhaseV)

      displayGradsD=numpy.zeros([2*nfft**2],numpy.float64)
      thisIdx=( (gO.illuminatedCentresIdx//(nfft+2)-1)*nfft
               +(gO.illuminatedCentresIdx%(nfft+2)-1) )
      displayGradsD[thisIdx.tolist()+(thisIdx+nfft**2).tolist()]=gradV
      displayGradsD.resize([2,nfft,nfft])
      displayGradsD=\
         numpy.ma.masked_array( displayGradsD, [pupilMask==0,pupilMask==0] )
      
      # recast phase and gradients onto square matrices
      calcdGradsD=numpy.ma.masked_array( numpy.array([
         0.5*(thisPhase[1:-1,2:]-thisPhase[1:-1,:-2]),
         0.5*(thisPhase[2:,1:-1]-thisPhase[:-2,1:-1])]),
            [pupilMask==0,pupilMask==0] )

      return calcdGradsD, displayGradsD

   def checkType2F():
      # checking code, Type 2 (Fried)
      gO=gradientOperatorType2Fried(pupilMask)
     
         # \/ constrict the phase screen by averaging (makes it
         #  consistent with Type 3
      thisPhase=\
         onePhase[1:,1:]+onePhase[:-1,1:]+onePhase[:-1,:-1]+onePhase[1:,:-1]
      onePhaseV=thisPhase.ravel()[gO.illuminatedCornersIdx] # phases

      gM=gO.returnOp()
      gradV=numpy.dot(gM,onePhaseV)
      
      displayGradsD=numpy.zeros([2*nfft**2],numpy.float64)
      displayGradsD[
       numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1).tolist()+
       (nfft**2+numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1)).tolist()]\
            =gradV
      displayGradsD.resize([2,nfft,nfft])
      displayGradsD=\
         numpy.ma.masked_array( displayGradsD, [pupilMask==0,pupilMask==0] )

      # theoretical gradients 
      calcdGradsD=numpy.ma.masked_array(
         numpy.array([ thisPhase[1:,1:]-thisPhase[:-1,:-1],
                       thisPhase[:-1,1:]-thisPhase[1:,:-1] ]
        ), [pupilMask==0,pupilMask==0] )

      return calcdGradsD, displayGradsD

   def display(title, calcdGradsD, displayGradsD):
      pg.figure()
      pg.subplot(221)
      pg.imshow( displayGradsD[0], interpolation='nearest' )
      pg.title(title+": gradOp, x")
      pg.subplot(222)
      pg.imshow( calcdGradsD[0], interpolation='nearest' )
      pg.title("direct, x")

      pg.subplot(223)
      pg.imshow( displayGradsD[1], interpolation='nearest' )
      pg.xlabel("gradOp, y")
      pg.subplot(224)
      pg.imshow( calcdGradsD[1], interpolation='nearest' )
      pg.title("direct, y")

   for dat in (checkType1,'type 1'),(checkType2F,'type 2, Fried'),\
            (checkType2,'type 2'),(checkType3C,'type 3, centred'):
      theory,opVals=dat[0]()
      display(dat[1], theory,opVals)
   import sys
   print("Waiting for button press..."), ; sys.stdout.flush()
   pg.waitforbuttonpress()
   print("(done)")
