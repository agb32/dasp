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
#(deprecated):      if self.numberSubaps>1000 and not self.sparse:
#(deprecated):         raise ValueError("Insist on sparsity for numberSubaps>1000")
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

#(deprecated):class gradientOperatorType2Fried(gradientOperatorType1):
#(deprecated):   '''Using Type 2 geometry to define a gradient operator matrix.
#(deprecated):    Can rotate the gradient basis to get a Fried-like geometry so this
#(deprecated):    is a straightforward approach which can promote waffle.'''
#(deprecated):
#(deprecated):   def calcOp_NumpyArray(self):
#(deprecated):      # Gradient matrix
#(deprecated):      # Similar to Fried but rotation of phases
#(deprecated):      self.op=numpy.zeros(
#(deprecated):         [2*self.numberSubaps,self.numberPhases],numpy.float64)
#(deprecated):      for i in range(self.numberSubaps):
#(deprecated):         # want for each sub-aperture, the corresponding corners
#(deprecated):         r=[ numpy.flatnonzero(self.illuminatedCornersIdx==tm)[0] for tm in
#(deprecated):            (self.mask2allcorners[i],
#(deprecated):             self.mask2allcorners[i]+1,
#(deprecated):             self.mask2allcorners[i]+self.n_[1],
#(deprecated):             self.mask2allcorners[i]+self.n_[1]+1) ]
#(deprecated):            # \/ grad x1->y2
#(deprecated):         self.op[i,r[0]]=-1
#(deprecated):         self.op[i,r[3]]=1
#(deprecated):            # \/ self.grad x2->y1
#(deprecated):         self.op[i+self.numberSubaps,r[1]]=1
#(deprecated):         self.op[i+self.numberSubaps,r[2]]=-1

#(deprecated):class gradientOperatorType2(gradientOperatorType1):
#(deprecated):   '''Using Type 2 geometry to define a gradient operator matrix.
#(deprecated):    Define explicitly.
#(deprecated):    The subaperture mask is now not equal for the x and y directions,
#(deprecated):    so the interpretation is that it is extended by one in the x and y
#(deprecated):    directions for the respective gradients (meaning they aren't identical).
#(deprecated):    The class must therefore provide a method for extending the subaperture
#(deprecated):    mask from that given which is some 'mean' value.'''
#(deprecated):   subapYMask=None
#(deprecated):   subapXMask=None
#(deprecated):   nx=None ; ny=None
#(deprecated):   numberXSubaps=None
#(deprecated):   numberYSubaps=None
#(deprecated):
#(deprecated):   def extendSubapMask(self, subapMask):
#(deprecated):      # in the Hudgin geometry, there are more measurements for x in the y
#(deprecated):      # direction and vice-versa.
#(deprecated):      self.n=subapMask.shape
#(deprecated):      
#(deprecated):      self.subapXMask=numpy.zeros([self.n[1]+1,self.n[0]],numpy.bool)
#(deprecated):      self.subapXMask[:-1]=subapMask
#(deprecated):      self.subapXMask[1:]+=subapMask 
#(deprecated):      
#(deprecated):      self.subapYMask=numpy.zeros([self.n[1],self.n[0]+1],numpy.bool)
#(deprecated):      self.subapYMask[:,:-1]=subapMask
#(deprecated):      self.subapYMask[:,1:]+=subapMask 
#(deprecated):
#(deprecated):      self.nx=self.subapXMask.shape
#(deprecated):      self.ny=self.subapYMask.shape
#(deprecated):
#(deprecated):
#(deprecated):   def newSubaperturesGiven(self, subapMask):
#(deprecated):      self.op=None # reset
#(deprecated):      self.extendSubapMask(subapMask)
#(deprecated):      self.n=subapMask.shape
#(deprecated):      self.n_=[ x+1 for x in self.n]
#(deprecated):
#(deprecated):      # now define the illuminated sub-apertures (-1 = not illuminated)
#(deprecated):      self.maskYIdx= numpy.arange(self.ny[0]*self.ny[1]).reshape(self.ny)\
#(deprecated):                     *self.subapYMask+(self.subapYMask-1)
#(deprecated):      self.maskXIdx= numpy.arange(self.nx[0]*self.nx[1]).reshape(self.nx)\
#(deprecated):                     *self.subapXMask+(self.subapXMask-1)
#(deprecated):         # \/ corners on the array
#(deprecated):      self.cornersIdx=numpy.arange(self.n_[0]*self.n_[1]).reshape(self.n_) 
#(deprecated):         # \/ per corner, how many subapertures use this corner
#(deprecated):      self.illuminatedCorners=numpy.zeros(self.n_,numpy.complex64)
#(deprecated): 
#(deprecated):      # \/ array to specify the first corner for each sub-aperture
#(deprecated):      self.maskX2allcorners=\
#(deprecated):         numpy.flatnonzero(self.subapXMask.ravel())%(self.nx[1])\
#(deprecated):        +(numpy.flatnonzero(self.subapXMask.ravel())//(self.nx[0]-1))*self.n_[0]
#(deprecated):      self.maskY2allcorners=\
#(deprecated):         numpy.flatnonzero(self.subapYMask.ravel())%(self.ny[1])\
#(deprecated):        +(numpy.flatnonzero(self.subapYMask.ravel())//(self.ny[0]+1))*self.n_[0]
#(deprecated):      # \/ figure out which phases then contribute to a sub-aperture gradient
#(deprecated):      for i in range(self.nx[0]):
#(deprecated):        for j in range(self.nx[1]):
#(deprecated):          maskVal=self.maskXIdx[i,j]
#(deprecated):          if maskVal==-1: continue # skip not illuminated
#(deprecated):          self.illuminatedCorners+=(
#(deprecated):            (self.cornersIdx//self.n_[1]==maskVal//self.nx[1])*\
#(deprecated):            (abs(self.cornersIdx%self.n_[0]-0.5-maskVal%(self.nx[0]-1))<=0.5) )
#(deprecated):      for i in range(self.ny[0]):
#(deprecated):        for j in range(self.ny[1]):
#(deprecated):          maskVal=self.maskYIdx[i,j]
#(deprecated):          if maskVal==-1: continue # skip not illuminated
#(deprecated):          self.illuminatedCorners+=1.0j*(
#(deprecated):            (self.cornersIdx%self.n_[0]==maskVal%(self.ny[0]+1))*\
#(deprecated):            (abs(self.cornersIdx//self.n_[1]-0.5-maskVal//self.ny[1])<=0.5) )
#(deprecated):         # \/ have to compare illumination via this slightly strage method
#(deprecated):      if (self.illuminatedCorners.real>0).sum() !=\
#(deprecated):         (self.illuminatedCorners.imag>0).sum():
#(deprecated):         raise ValueError("Failure of commonality between phase points")
#(deprecated):      self.illuminatedCornersIdx=\
#(deprecated):         numpy.flatnonzero(self.illuminatedCorners.ravel())
#(deprecated):
#(deprecated):      # \/ some handy numbers
#(deprecated):      self.numberXSubaps=int(self.subapXMask.sum())
#(deprecated):      self.numberYSubaps=int(self.subapYMask.sum()) 
#(deprecated):      self.numberPhases=(self.illuminatedCorners!=0).sum() # =(n+1)**2 for all illum.
#(deprecated):
#(deprecated):
#(deprecated):   def returnOp(self):
#(deprecated):      if self.numberXSubaps==None: return None
#(deprecated):      if self.op!=None: return self.op
#(deprecated):      self.calcOp()
#(deprecated):      return self.op
#(deprecated):
#(deprecated):   def calcOp_NumpyArray(self):
#(deprecated):      self.op=numpy.zeros( [self.numberXSubaps+self.numberYSubaps,
#(deprecated):         self.numberPhases],numpy.float64)
#(deprecated):      for i in range(self.numberXSubaps):
#(deprecated):         # want for each sub-aperture, the corresponding corners
#(deprecated):         r=[ numpy.flatnonzero(self.illuminatedCornersIdx==tm)[0] for tm in
#(deprecated):            (self.maskX2allcorners[i], self.maskX2allcorners[i]+1) ]
#(deprecated):         self.op[i,r[0]]=-1
#(deprecated):         self.op[i,r[1]]=1
#(deprecated):
#(deprecated):      for i in range(self.numberYSubaps):
#(deprecated):         # want for each sub-aperture, the corresponding corners
#(deprecated):         r=[ numpy.flatnonzero(self.illuminatedCornersIdx==tm)[0] for tm in
#(deprecated):            (self.maskY2allcorners[i], self.maskY2allcorners[i]+self.n_[1]) ]
#(deprecated):         self.op[i+self.numberXSubaps,r[0]]=-1
#(deprecated):         self.op[i+self.numberXSubaps,r[1]]=1

# ------------------------------


#(deprecated):class gradientOperatorType3Avg(gradientOperatorType1):
#(deprecated):   '''Using Type 3 geometry and averaging gradients, define a gradient operator
#(deprecated):     matrix.'''
#(deprecated):   # Convert the Southwell geometry into a Hudgin one via an averaging matrix,
#(deprecated):   # so overload the Hudgin operator.
#(deprecated):   def __init__(self, subapMask):
#(deprecated):      raise NotimplementedError("Type 3 with averaging not yet implemented")
#(deprecated):

#(deprecated):class gradientOperatorType3Centred(gradientOperatorType1):
#(deprecated):   '''Using Type 3 geometry and fake sub-apertures, define a gradient operator
#(deprecated):     matrix.'''
#(deprecated):   # Change from type 1 is that there isn't enough information for corners so
#(deprecated):   # fake phases must be added (these can be arbirtrary) but the responsibility
#(deprecated):   # for this lies with the constructed phase vector This makes this code far
#(deprecated):   # simpler but we ought to also supply a member that suitably constricts and
#(deprecated):   # expands a phase vector to the 'right' length.
#(deprecated):
#(deprecated):   def newSubaperturesGiven(self, subapMask):
#(deprecated):      self.op=None # reset
#(deprecated):      self.n=subapMask.shape
#(deprecated):      self.n_=[ x+2 for x in self.n] # pad by an extra row and column
#(deprecated):
#(deprecated):      # now define the illuminated sub-apertures (-1 = not illuminated)
#(deprecated):      self.subapMaskSequential=numpy.arange(
#(deprecated):            self.n[0]*self.n[1]).reshape(self.n)*subapMask+(subapMask-1)
#(deprecated):         # \/ per phase centre, how many subapertures use this corner
#(deprecated):      self.allCentres=numpy.zeros(self.n_)
#(deprecated):      self.illuminatedCentres=numpy.zeros(self.n_)
#(deprecated):      self.illuminatedCentres[1:-1,1:-1]=subapMask
#(deprecated):         # \/ zip around the mask and add 'fake' subaperture where
#(deprecated):         #  reqd.
#(deprecated):      for i in range(self.n[0]):
#(deprecated):         for j in range(self.n[1]):
#(deprecated):            if subapMask[i,j]:
#(deprecated):               for k in (0,2):
#(deprecated):                  self.allCentres[i+k,j+1]+=1
#(deprecated):                  self.allCentres[i+1,j+k]+=1
#(deprecated):
#(deprecated):      self.illuminatedCentresIdx=\
#(deprecated):         numpy.flatnonzero(self.illuminatedCentres.ravel())
#(deprecated):      self.allCentresIdx=numpy.flatnonzero(self.allCentres.ravel())
#(deprecated):
#(deprecated):      # \/ now figure out the real and fake subaperture indices
#(deprecated):      # now define the illuminated sub-apertures (-1 = not illuminated,0
#(deprecated):      # fake)
#(deprecated):      self.subapMaskSequential=-numpy.ones(self.n_)
#(deprecated):      self.subapMaskSequential.ravel()[self.allCentresIdx]=0
#(deprecated):      self.subapMaskSequential.ravel()[self.illuminatedCentresIdx]=1
#(deprecated):         # \/ index into the phase vectors, to insert/extract real values
#(deprecated):         #  into/from fake values vector
#(deprecated):      self.fake2real=numpy.flatnonzero(
#(deprecated):         self.subapMaskSequential.ravel()[self.allCentresIdx]==1)
#(deprecated):      self.real2fake=numpy.flatnonzero(self.subapMaskSequential.ravel()==1)
#(deprecated):
#(deprecated):      # some handy numbers
#(deprecated):      self.numberSubaps=subapMask.sum() # =n**2 for all valid illuminated
#(deprecated):      self.numberPhases=(self.allCentres!=0).sum() 
#(deprecated):         # =(n+2)**2 for all illum., this is of course the fake version
#(deprecated):
#(deprecated):   def testInit(self, subapMask):
#(deprecated):      import matplotlib.pyplot as pg
#(deprecated):      # Test
#(deprecated):      pg.imshow(subapMask,interpolation='nearest',
#(deprecated):         extent=[0.5,self.n[1]+0.5,0.5,self.n[0]+0.5])
#(deprecated):      pg.axis([0,self.n_[1],0,self.n_[0]])
#(deprecated):      pg.plot( numpy.flatnonzero(self.allCentres.ravel())%self.n_[1],
#(deprecated):               numpy.flatnonzero(self.allCentres.ravel())//self.n_[0],
#(deprecated):                  'wo' )
#(deprecated):      pg.title("red=illuminated, white=their corners")
#(deprecated):      pg.waitforbuttonpress()
#(deprecated):   
#(deprecated):   def constrictPhase(self,phaseV):
#(deprecated):      '''Given a fake phase vector, extract the real values'''
#(deprecated):      return numpy.take(phaseV,self.fake2real,axis=-1)
#(deprecated):   def expandPhase(self,phaseV):
#(deprecated):      '''Given a real phase vector, create a fake one'''
#(deprecated):      newPhaseV=numpy.zeros(self.numberPhases)
#(deprecated):      newPhaseV[self.fake2real]=phaseV
#(deprecated):      return newPhaseV
#(deprecated):   
#(deprecated):   def calcOp_NumpyArray(self):
#(deprecated):      # gradient matrix
#(deprecated):      self.op=numpy.zeros(
#(deprecated):         [2*self.numberSubaps,self.numberPhases],numpy.float64)
#(deprecated):      for i in range(self.numberSubaps):
#(deprecated):         # want for each /real/ sub-aperture, the centred difference
#(deprecated):         r=[ numpy.flatnonzero(self.allCentresIdx==tm)[0] for tm in
#(deprecated):            (self.illuminatedCentresIdx[i]-1,
#(deprecated):             self.illuminatedCentresIdx[i]+1,
#(deprecated):             self.illuminatedCentresIdx[i]-self.n_[1],
#(deprecated):             self.illuminatedCentresIdx[i]+self.n_[1]) ]
#(deprecated):            # \/ grad x
#(deprecated):         self.op[i,r[0]]=-0.5
#(deprecated):         self.op[i,r[1]]=0.5
#(deprecated):            # \/ self.grad y
#(deprecated):         self.op[i+self.numberSubaps,r[2]]=-0.5
#(deprecated):         self.op[i+self.numberSubaps,r[3]]=0.5

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

#(deprecated):# The approximation to the kolmogorov inverse stencil, based on explicit
#(deprecated):# SVD-based inversion is, with the centre weighted at +4 c.f. the laplacian,
#(deprecated):#
#(deprecated):#([[ -1.8e-03,   9.1e-02,   2.5e-01,   9.1e-02,  -2.9e-03],
#(deprecated):#  [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
#(deprecated):#  [  2.5e-01,  -1.3e+00,   4.0e+00,  -1.3e+00,   2.5e-01],
#(deprecated):#  [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
#(deprecated):#  [ -2.9e-03,   8.9e-02,   2.5e-01,   8.9e-02,  -4.2e-03]]) 
#(deprecated):
#(deprecated):kolmogInverseStencil=numpy.array(
#(deprecated):      [[ -1.8e-03,   9.1e-02,   2.5e-01,   9.1e-02,  -2.9e-03],
#(deprecated):       [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
#(deprecated):       [  2.5e-01,  -1.3e+00,   4.0e+00,  -1.3e+00,   2.5e-01],
#(deprecated):       [  9.1e-02,  -3.8e-02,  -1.3e+00,  -3.8e-02,   8.9e-02],
#(deprecated):       [ -2.9e-03,   8.9e-02,   2.5e-01,   8.9e-02,  -4.2e-03]]).ravel()
#(deprecated):
#(deprecated):def genericKolmogInverse_findLocation(self, i):
#(deprecated):   thisIdx=self.illuminatedCornersIdx[i]
#(deprecated):   wantedIdx=thisIdx+\
#(deprecated):       numpy.add.outer(numpy.arange(-2,3)*self.n_[1],numpy.arange(-2,3)).ravel()
#(deprecated):      # \/ check if found indices are in the right place
#(deprecated):   rowcolIdx=[ numpy.add.outer(numpy.arange(-2,3)+thisIdx//self.n_[1],
#(deprecated):                  numpy.zeros(5)).ravel(),
#(deprecated):               numpy.add.outer(numpy.zeros(5),
#(deprecated):                  numpy.arange(-2,3)+thisIdx%self.n_[1]).ravel() ]
#(deprecated):   valid=[ numpy.flatnonzero(
#(deprecated):         (self.illuminatedCornersIdx==wantedIdx[i])*
#(deprecated):         ((wantedIdx[i]//self.n_[1])==rowcolIdx[0][i])*
#(deprecated):         ((wantedIdx[i]%self.n_[1])==rowcolIdx[1][i])
#(deprecated):         )
#(deprecated):      for i in range(25) ]
#(deprecated):   return valid

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

#(deprecated):def genericKolmogInverseCalcOp_NumpyArray(operator):
#(deprecated):   self=operator
#(deprecated):   self.op=numpy.zeros( [self.numberPhases]*2, numpy.float64 )
#(deprecated):   for i in range(self.numberPhases):
#(deprecated):      valid=genericKolmogInverse_findLocation(self,i)
#(deprecated):      for j in range(25):
#(deprecated):         if valid[j].shape[0]==1:
#(deprecated):            loc=valid[j][0]
#(deprecated):            self.op[i, loc ]=kolmogInverseStencil[j]

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

class laplacianOperatorType1(gradientOperatorType1):
   '''Define an operator that computes the Laplacian'''
   def calcOp_NumpyArray(self):
      genericLaplacianCalcOp_NumpyArray(self)
   def calcOp_scipyCSR(self):
      genericLaplacianCalcOp_scipyCSR(self)

#(deprecated):class laplacianOperatorType2(gradientOperatorType2):
#(deprecated):   '''Define an operator that computes the Laplacian, for Type 2'''
#(deprecated):   def calcOp_NumpyArray(self):
#(deprecated):      genericLaplacianCalcOp_NumpyArray(self)

#(deprecated):class laplacianOperatorType3Centred(gradientOperatorType3Centred):
#(deprecated):   '''Define an operator that computes the Laplacian, for Type 3'''
#(deprecated):   def calcOp_NumpyArray(self):
#(deprecated):      genericLaplacianCalcOp_NumpyArray(self)

#(deprecated):class kolmogInverseOperatorType1(gradientOperatorType1):
#(deprecated):   '''Define an operator that computes an approximation to the -11/3
#(deprecated):      power spectrum'''
#(deprecated):   def calcOp_NumpyArray(self):
#(deprecated):      genericKolmogInverseCalcOp_NumpyArray(self)

# ------------------------------


#(deprecated):class localMaskOperatorType1(gradientOperatorType1):
#(deprecated):   '''Define an operator that incorporates only certain distances within an
#(deprecated):   operator matrix.'''
#(deprecated):   def __init__( self, subapMask=None, distances=[0] ):
#(deprecated):      self.distances=distances
#(deprecated):      gradientOperatorType1.__init__(self,subapMask)
#(deprecated):
#(deprecated):   def calcOp_NumpyArray(self):
#(deprecated):      self.op=numpy.zeros( [self.numberPhases]*2, numpy.float64 )
#(deprecated):      # distance stencil is,
#(deprecated):      #  for every distance,
#(deprecated):      #   if dist not 0:
#(deprecated):      #      wantedIdx+=[
#(deprecated):      # thisIdx-dist*self.n_[1],thisIdx-dist,thisIdx+dist,thisIdx+dist*self.n_[1] ]
#(deprecated):      #   else:
#(deprecated):      #      wantedIdx+=[thisIdx]
#(deprecated):      #  valid=[ twi in self.illuminatedCornersIdx for twi in wantedIdx ]
#(deprecated):      for i in range(self.numberPhases):
#(deprecated):         thisIdx=self.illuminatedCornersIdx[i]
#(deprecated):         self.finalWantedIdx=[]
#(deprecated):         for td in self.distances:
#(deprecated):            if td!=0:
#(deprecated):               wantedIdx=[ thisIdx-td*self.n_[1],thisIdx-td,
#(deprecated):                  thisIdx+td,thisIdx+td*self.n_[1] ]
#(deprecated):               valid=numpy.array([
#(deprecated):                  twi in self.illuminatedCornersIdx for twi in wantedIdx ]
#(deprecated):                     ).nonzero()[0]
#(deprecated):               wantedIdx=numpy.take(wantedIdx,valid)
#(deprecated):            else:
#(deprecated):               wantedIdx=[ thisIdx ]
#(deprecated):            self.finalWantedIdx+=numpy.searchsorted(
#(deprecated):               self.illuminatedCornersIdx, wantedIdx ).tolist()
#(deprecated):         for j in self.finalWantedIdx:
#(deprecated):            self.op[i, j]+=1


#(deprecated):if __name__=="__main__":
#(deprecated):   nfft=7 # pixel size
#(deprecated):   onePhase=numpy.random.normal(size=[nfft+2]*2 ) # Type 1 & Type 3
#(deprecated):   # define pupil mask as sub-apertures
#(deprecated):   pupilMask=numpy.ones([nfft]*2,numpy.int32)
#(deprecated):   pupilCds=numpy.add.outer(
#(deprecated):      (numpy.arange(nfft)-(nfft-1)/2.)**2.0,
#(deprecated):      (numpy.arange(nfft)-(nfft-1)/2.)**2.0 )
#(deprecated):   pupilMask=(pupilCds>2**2)*(pupilCds<4**2)
#(deprecated):
#(deprecated):   # test, 06/Aug/2012 to confirm UB's issue for CURE
#(deprecated):   pupilMask=numpy.array([[0,0,1,1,1,0,0], [0,1,1,1,1,1,0], [1,1,1,1,1,1,1],
#(deprecated):            [1,1,1,0,1,1,1], [1,1,1,1,1,1,1], [0,1,1,1,1,1,0], [0,0,1,1,1,0,0]])
#(deprecated):
#(deprecated):   def checkType1():
#(deprecated):      # checking code, Type 1
#(deprecated):      gO=gradientOperatorType1(pupilMask)
#(deprecated):     
#(deprecated):         # \/ constrict the phase screen by averaging (makes it
#(deprecated):         #  consistent with Type 3
#(deprecated):      thisPhase=\
#(deprecated):         onePhase[1:,1:]+onePhase[:-1,1:]+onePhase[:-1,:-1]+onePhase[1:,:-1]
#(deprecated):      onePhaseV=thisPhase.ravel()[gO.illuminatedCornersIdx] # phases
#(deprecated):
#(deprecated):      gM=gO.returnOp()
#(deprecated):      gradV=numpy.dot(gM,onePhaseV)
#(deprecated):      
#(deprecated):      displayGradsD=numpy.zeros([2*nfft**2],numpy.float64)
#(deprecated):      displayGradsD[
#(deprecated):       numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1).tolist()+
#(deprecated):       (nfft**2+numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1)).tolist()]\
#(deprecated):        =gradV
#(deprecated):      displayGradsD.resize([2,nfft,nfft])
#(deprecated):      displayGradsD=\
#(deprecated):         numpy.ma.masked_array( displayGradsD, [pupilMask==0,pupilMask==0] )
#(deprecated):
#(deprecated):      # theoretical gradients 
#(deprecated):      calcdGradsD=numpy.ma.masked_array(numpy.array([
#(deprecated):        0.5*(thisPhase[1:,1:]-thisPhase[1:,:-1])
#(deprecated):         +0.5*(thisPhase[:-1,1:]-thisPhase[:-1,:-1]),
#(deprecated):        0.5*(thisPhase[1:,1:]+thisPhase[1:,:-1])
#(deprecated):         -0.5*(thisPhase[:-1,1:]+thisPhase[:-1,:-1]) ]
#(deprecated):        ), [pupilMask==0,pupilMask==0] )
#(deprecated):
#(deprecated):      return calcdGradsD, displayGradsD
#(deprecated):   
#(deprecated):   def checkType2():
#(deprecated):      # checking code, Type 2
#(deprecated):      gO=gradientOperatorType2(pupilMask)
#(deprecated):     
#(deprecated):         # \/ constrict the phase screen by averaging (makes it
#(deprecated):         #  consistent with Type 3
#(deprecated):      thisPhase=\
#(deprecated):         onePhase[1:,1:]+onePhase[:-1,1:]+onePhase[:-1,:-1]+onePhase[1:,:-1]
#(deprecated):      onePhaseV=thisPhase.ravel()[gO.illuminatedCornersIdx] # phases
#(deprecated):
#(deprecated):      gM=gO.returnOp()
#(deprecated):      gradV=numpy.dot(gM,onePhaseV)
#(deprecated):      
#(deprecated):      displayGradsD=[ numpy.zeros([(nfft+1)*nfft],numpy.float64),
#(deprecated):                      numpy.zeros([(nfft+1)*nfft],numpy.float64) ]
#(deprecated):      displayGradsD[0][numpy.flatnonzero(gO.maskXIdx.ravel()>-1)]=\
#(deprecated):         gradV[:gO.numberXSubaps]
#(deprecated):      displayGradsD[1][numpy.flatnonzero(gO.maskYIdx.ravel()>-1)]=\
#(deprecated):         gradV[gO.numberXSubaps:]
#(deprecated):      displayGradsD[0].resize([nfft+1,nfft])
#(deprecated):      displayGradsD[1].resize([nfft,nfft+1])
#(deprecated):      displayGradsD[0]=\
#(deprecated):         numpy.ma.masked_array( displayGradsD[0], gO.subapXMask==0 )
#(deprecated):      displayGradsD[1]=\
#(deprecated):         numpy.ma.masked_array( displayGradsD[1], gO.subapYMask==0 )
#(deprecated):
#(deprecated):      # theoretical gradients 
#(deprecated):      calcdGradsD=[
#(deprecated):         numpy.ma.masked_array(thisPhase[:,1:]-thisPhase[:,:-1],
#(deprecated):            gO.subapXMask==0 ),\
#(deprecated):         numpy.ma.masked_array(thisPhase[1:]-thisPhase[:-1],
#(deprecated):            gO.subapYMask==0 )
#(deprecated):            ]
#(deprecated):      return calcdGradsD, displayGradsD
#(deprecated):
#(deprecated):
#(deprecated):
#(deprecated):   def checkType3C():
#(deprecated):      # checking code, Type 3, centred
#(deprecated):      gO=gradientOperatorType3Centred(pupilMask)
#(deprecated):     
#(deprecated):      thisPhase=onePhase # the edge values are never used
#(deprecated):      onePhaseV=thisPhase.ravel()[gO.allCentresIdx]
#(deprecated):      
#(deprecated):      gM=gO.returnOp()
#(deprecated):      gradV=numpy.dot(gM,onePhaseV)
#(deprecated):
#(deprecated):      displayGradsD=numpy.zeros([2*nfft**2],numpy.float64)
#(deprecated):      thisIdx=( (gO.illuminatedCentresIdx//(nfft+2)-1)*nfft
#(deprecated):               +(gO.illuminatedCentresIdx%(nfft+2)-1) )
#(deprecated):      displayGradsD[thisIdx.tolist()+(thisIdx+nfft**2).tolist()]=gradV
#(deprecated):      displayGradsD.resize([2,nfft,nfft])
#(deprecated):      displayGradsD=\
#(deprecated):         numpy.ma.masked_array( displayGradsD, [pupilMask==0,pupilMask==0] )
#(deprecated):      
#(deprecated):      # recast phase and gradients onto square matrices
#(deprecated):      calcdGradsD=numpy.ma.masked_array( numpy.array([
#(deprecated):         0.5*(thisPhase[1:-1,2:]-thisPhase[1:-1,:-2]),
#(deprecated):         0.5*(thisPhase[2:,1:-1]-thisPhase[:-2,1:-1])]),
#(deprecated):            [pupilMask==0,pupilMask==0] )
#(deprecated):
#(deprecated):      return calcdGradsD, displayGradsD
#(deprecated):
#(deprecated):   def checkType2F():
#(deprecated):      # checking code, Type 2 (Fried)
#(deprecated):      gO=gradientOperatorType2Fried(pupilMask)
#(deprecated):     
#(deprecated):         # \/ constrict the phase screen by averaging (makes it
#(deprecated):         #  consistent with Type 3
#(deprecated):      thisPhase=\
#(deprecated):         onePhase[1:,1:]+onePhase[:-1,1:]+onePhase[:-1,:-1]+onePhase[1:,:-1]
#(deprecated):      onePhaseV=thisPhase.ravel()[gO.illuminatedCornersIdx] # phases
#(deprecated):
#(deprecated):      gM=gO.returnOp()
#(deprecated):      gradV=numpy.dot(gM,onePhaseV)
#(deprecated):      
#(deprecated):      displayGradsD=numpy.zeros([2*nfft**2],numpy.float64)
#(deprecated):      displayGradsD[
#(deprecated):       numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1).tolist()+
#(deprecated):       (nfft**2+numpy.flatnonzero(gO.subapMaskSequential.ravel()>-1)).tolist()]\
#(deprecated):            =gradV
#(deprecated):      displayGradsD.resize([2,nfft,nfft])
#(deprecated):      displayGradsD=\
#(deprecated):         numpy.ma.masked_array( displayGradsD, [pupilMask==0,pupilMask==0] )
#(deprecated):
#(deprecated):      # theoretical gradients 
#(deprecated):      calcdGradsD=numpy.ma.masked_array(
#(deprecated):         numpy.array([ thisPhase[1:,1:]-thisPhase[:-1,:-1],
#(deprecated):                       thisPhase[:-1,1:]-thisPhase[1:,:-1] ]
#(deprecated):        ), [pupilMask==0,pupilMask==0] )
#(deprecated):
#(deprecated):      return calcdGradsD, displayGradsD
#(deprecated):
#(deprecated):   def display(title, calcdGradsD, displayGradsD):
#(deprecated):      pg.figure()
#(deprecated):      pg.subplot(221)
#(deprecated):      pg.imshow( displayGradsD[0], interpolation='nearest' )
#(deprecated):      pg.title(title+": gradOp, x")
#(deprecated):      pg.subplot(222)
#(deprecated):      pg.imshow( calcdGradsD[0], interpolation='nearest' )
#(deprecated):      pg.title("direct, x")
#(deprecated):
#(deprecated):      pg.subplot(223)
#(deprecated):      pg.imshow( displayGradsD[1], interpolation='nearest' )
#(deprecated):      pg.xlabel("gradOp, y")
#(deprecated):      pg.subplot(224)
#(deprecated):      pg.imshow( calcdGradsD[1], interpolation='nearest' )
#(deprecated):      pg.title("direct, y")
#(deprecated):
#(deprecated):   for dat in (checkType1,'type 1'),(checkType2F,'type 2, Fried'),\
#(deprecated):            (checkType2,'type 2'),(checkType3C,'type 3, centred'):
#(deprecated):      theory,opVals=dat[0]()
#(deprecated):      display(dat[1], theory,opVals)
#(deprecated):   import sys
#(deprecated):   print("Waiting for button press..."), ; sys.stdout.flush()
#(deprecated):   pg.waitforbuttonpress()
#(deprecated):   print("(done)")
