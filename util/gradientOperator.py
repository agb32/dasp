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
      if subapMask is not None: self.newSubaperturesGiven(subapMask)
      if pupilMask is not None: self.newPupilGiven(pupilMask)

   def newSubaperturesGiven(self, subapMask):
      self.numberSubaps=int(subapMask.sum()) # =n**2 for all illuminated
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


def genericLaplacianCalcOp_NumpyArray(operator):
   self=operator
   self.op=numpy.zeros( [self.numberPhases]*2, numpy.float32 )
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
   self.op=scipy.sparse.csr_matrix( (data, (row,column)), dtype=numpy.float32 )

class laplacianOperatorType1(gradientOperatorType1):
   '''Define an operator that computes the Laplacian'''
   def calcOp_NumpyArray(self):
      genericLaplacianCalcOp_NumpyArray(self)
   def calcOp_scipyCSR(self):
      genericLaplacianCalcOp_scipyCSR(self)

