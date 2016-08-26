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
#A module for pre-conditioned conjugate gradient algorithms.
import numpy
import numpy.fft as fft
import types
import util.FITS,util.phaseCovariance
try:
    import gist
except:
    pass
import util.createPokeMx
import science.infScrn
import util.spmatrix
import scipy.sparse
#import util.dot as quick
"""Aims to solve Ax=b.
For a tomographic AO system.

For a large SPD (symmetric positive definite) system.  PCG is a
standard conjugate gradient algorithm applied to a symmetrized version
of the transformed system:

C-1Ax=C-1b

where the preconditioning matrix C is also SPD.

For this to be efficient, matrix vector products C-1r must be
inexpsnsive to compute and C-1A myust have a desirable eigenvalue
distribution in the sense that eigenvalues are clustered and/or have a
relatively small range.  The PCG algorithm can then be applied.

See Fourier domain preconditioned conjugate gradient algorithm for
atmospheric tomography, Q Yang, C Vogel, B Ellerbroek for further
details (2006)

When using as part of the simulation, you are probably best to call
the createPcg method initially, use pcg.newPokemx to add a new poke
matrix, and use pcg.solve each simulation iteration.  Actually, no,
use the xinterp_recon object.

Need to make sure that all centroids etc are passed (ie don't mask any
out).  This also means that the pokemx should contain all centroids
and all actuators.  This is designed to work with Fried geometry
DM/SHS, and needs the mirror modes to be pokes, not eg zernikes etc
(ie zonal).  Everything needs to be passed because of the sparse
matrix approximations etc that are made, and which assume everything
is nice and regular.


The estimated phase (actuator values) will be placed in pcg.x

"""



class pcg:
    def __init__(self, noiseCov, phaseCov,centroids,dmList=None,ngsList=None,lgsList=None,minIter=10,maxIter=100,convergenceValue=5e-4,options=None,convergenceWarning=1,fakePokeMxVal=0.31,telDiam=42.):
        """
        Given the poke matrix etc, compute the A matrix and b vector
        which are used for the calculations.  The phaseCovMatrix would
        be the Noll matrix if a zernike DM etc was used (but this
        wouldn't work, as this method requires a modal
        representation...) anyway, if gives you an idea.
        
        Here, typically, A will be GT Cn-1 G + Cp-1 where G is the
        poke matrix (shape ncents, nmodes), Cn is the
        noise covariance matrix and Cp is the phase covariance matrix
        (Noll matrix).  b is GT Cn-1 s where s is the centroid
        measurements.

        A has shape nmodes, nmodes and b has shape nmodes.

        Inputs are:
         - the poke matrix with shape (ncents, nmodes) or None, in which case, a sparse representation will be generated here...
         - noiseCov, a diagonal matrix or vector  of shape ncents.
           Or, a vector with 1 element equal to the mean subaperture noise
           variance.
           Or a float, equal to mean subap noise variance.
         - phaseCovMatrix, of shape (nmodes, nmodes).
         - centroids, the centroid values from this simulation iteration.
         - actPosition (obsolete), the actuator position, an array of coordinates,
           such that [i][0] gives the x coord of the ith actuator and
           [i][1] gives the y coord of the ith actuator.  The coordinates
           should be scaled correctly so that this works(!!!), ie should
           be in the correct units.  Sorry thats not very helpful!
         - innerScale (obsolete) the minimum distance to assume between actuators
           (ie to stop a pole when comparing an actuator with itself).
           This can be none if an inner scale is not equal to "inf".
         - nact - the number of actuators in 1d.  Can use none if using same
           number in x and y (gets calculated automatically)
         - l0, the outer scale. Can be "inf" for kolmogorov turbulence.
         - sparsityWidth, the width of the sparse diagonal matrix that
           should be used.  If this value is a float, then this is a
           fraction of the total width, while if an int, this is the
           number of values.
         - turbStrength, the turbulence strength for each layer
           assumed in reconstruction.  If you only have 1 DM,
           turbStrength will only be 1 unit long.
         - minIter, the minimum number of iterations for the PCG algorithm
         - maxIter, the max number of iterations for the PCG algorithm
         - convergenceValue, the max difference between results from two
           consecutive PCG iterations allowed before the problem is considered
           to be converged.  
         - options, a dictionary containing the options for optimisation, eg diagonalising matricess at various points etc.
         - dmList - list of virtual DMs, util.createPokeMx.DM objects.
         - gsList - list of LGS and NGS.  Not used except to create the sparse theoretical poke matrix and compute ncent
           
        """
        self.options=options#a list of options for what is diagonalised etc - ie tuning the algorithm.
        if self.options==None:
            self.options={}
        defaultOptions={"diagonaliseinvChatReordered":1,#should we diagonalise C-1, or use the whole MX?
                        "diagonaliseA":0.5,#should we diagonalise A or use the whole MX.
                        "removeHiFreqZ":0,#whether to remove the high freq component of z each iteration.
                        "removeHiFreqP":0,#whether to remove the high freq component of p each iteration.
                        "removeHiFreqX":0,#whether to remove the high freq component of x each iteration.
                        "useSparsePokeMx":1#whether to assume pokemx is sparse.
                        }
        for key in defaultOptions.keys():
            if not self.options.has_key(key):
                self.options[key]=defaultOptions[key]
        self.defaultOptions=defaultOptions
        #nmodes=pokemx.shape[1]
        self.fakePokeMxVal=fakePokeMxVal
##         if type(turbStrength)==type([]):
##             self.nLayers=len(turbStrength)
##         elif type(turbStrength) in [type(1),type(1.)]:
##             self.nLayers=1
##         else:
##             self.nLayers=turbStrength.shape[0]
        self.nLayers=len(dmList)
        if len(dmList)!=self.nLayers:
            raise Exception("Wrong number of DMs or layers for tomofdpcg: %d %d"%(len(dmList),self.nLayers))
        self.convergenceWarning=convergenceWarning
        self.noiseCov=noiseCov
        self.telDiam=telDiam
        self.dmList=dmList
        self.ngsList=ngsList
        self.lgsList=lgsList
        nacts=0
        ncent=0
        for gs in ngsList+lgsList:
            ncent+=2*(gs.nsubx**2)
        for dm in dmList:
            nacts+=dm.nact**2
        self.nacts=nacts
        self.ncent=ncent
        
        #Assume noise covariance matrix is diagonal (it should be)... so
        #inverting it is just inverting the diagonal elements.
        if type(noiseCov)==type(0.):#just a single value...
            if noiseCov==0.:#not tested?
                self.invNoiseCovVector=0.
            else:
                self.invNoiseCovVector=1./noiseCov
        else:
            noiseCov=numpy.array(noiseCov)
            if len(noiseCov.shape)>1:#get into vector format.
                noiseCov=numpy.diagonal(noiseCov)
            invNoiseCovVector=numpy.where(noiseCov==0,0,1./noiseCov)
            self.invNoiseCovVector=invNoiseCovVector
        #self.sigmam2=numpy.average(invNoiseCovVector)
        #self.sigma=self.sigmam2**-0.5
        phaseCov=numpy.array(phaseCov)
        self.phaseCov=phaseCov
        if len(self.phaseCov.shape)==1:#already diagonalised...
            self.invPhaseCov=1./self.phaseCov
        else:
            self.invPhaseCov=numpy.linalg.inv(self.phaseCov)
        #self.actPosition=None#no longer used...
        print "updating pcg pokemx"
        self.newPokemxSparse()

        #Depending on what you use for the invPhaseCov matrix (ie the BCCB approximation), you can possibly diagonalise A which then makes the A matrix multiplication order N rather than N^2.  The BCCB approximation just assumes that the inverse phase covariance matrix is diagonal, which though not perfect, may be reasonable.  Actually, typically, the invPhaseCov added to GCG is fairly small, so doesn't have much effect - so usually, if GCG is approx diagonal, A can be diagonalised.
        if type(centroids)==types.NoneType:
            centroids=numpy.zeros((self.pokeNoiseSparse.mx.shape[1],),numpy.float32)
        else:
            centroids=numpy.array(centroids)
        self.centroids=centroids
        #self.b=quick.dot(self.pokeNoise,centroids)
        self.b=numpy.zeros((self.nacts,),self.centroids.dtype)
        print "computing pcg centroids"
        self.newCentroids(self.centroids)#compute self.b

        #Here, assume poke matrix has shape ncents,nmodes.
        #self.b=b#eg the pokemxT times noise covariance-1 times centroid measurements (G^T C^-1 s), shape nmodes.  Pokemx has shape ncents,nmodes.
        #self.A=A#eg the poke matrix transposed times inverse noise covariance matrix times poke matrix plus Noll matrix, final shape (nmodes, nmodes).
        self.x=numpy.zeros((self.nacts,),self.b.dtype)#initial guess, shape nmodes.
        self.estPhase=[]
        layerStart=0
        for i in xrange(self.nLayers):
            layerEnd=self.dmList[i].nact**2
            self.estPhase.append(self.x[layerStart:layerStart+layerEnd].view())
            self.estPhase[-1].shape=(self.dmList[i].nact,self.dmList[i].nact)
            layerStart+=layerEnd
        self.minIter=minIter
        self.maxIter=maxIter
        self.tolerance=numpy.zeros((maxIter,),numpy.float32)
        self.alphaHist=numpy.zeros((maxIter,),numpy.float32)
        self.betaHist=numpy.zeros((maxIter,),numpy.float32)
        self.convergenceValue=convergenceValue
        self.convVal=0.
        print "finished pcg initialisation"
        

    def newPokemxSparse(self):
        """Pass a new poke matrix in... in sparse format... (so that larger problems can fit in memory!)
        This basically does the same as newPokemx for some specific allowed cases (should raise an error
        if not allowed), but only ever uses sparse matricees, to avoid having to create the full matrix
        (enough memory may not be available!).
        """
        for key in self.options.keys():
            if self.defaultOptions[key]!=self.options[key]:
                raise Exception("fdpcg - newPokemxSparse - must use default options...")
        spmx=util.createPokeMx.sparseTomo(ngsList=self.ngsList,lgsList=self.lgsList,dmList=self.dmList,telDiam=self.telDiam)
        spmx.mx.data*=self.fakePokeMxVal
        spmx2=spmx.mx
        spmx2.data=spmx.mx.data.copy()
        spmx2=spmx2.transpose()
        #spmx.create(self.fakePokeMxVal)
        #spmxT=util.createPokeMx.sparsePMX(nact)
        #spmxT.create(self.fakePokeMxVal)
        if type(self.invNoiseCovVector)==type(0.):
            spmx.mx.data*=self.invNoiseCovVector
        else:#size probably ncents... this is the diagonal of a matrix to multiply with spmx.
            incv=scipy.sparse.lil_matrix((self.ncent,self.ncent))
            incv.setdiag(self.invNoiseCovVector)
            incv=incv.tocsc()
            spmx.mx=spmx.dot(incv)#multiply each poke matrix value by the noise.
        self.pokeNoiseSparse=spmx
        self.GCGsparse=spmx.dot(spmx2)
        print "todo - check strip of gcgsparse is ok..."
        #util.spmatrix.strip(self.GCGsparse,0.1)
        if type(self.invPhaseCov)==type(0.) or len(self.invPhaseCov.shape)==1:#diagonal...
            #add invPhaseCov vector/value to the diagonal elements of GCGsparse
            print "invphasecov %s"%str( self.invPhaseCov)
            self.GCGsparse=util.spmatrix.addToDiag(self.GCGsparse,self.invPhaseCov)#add to the diagonal element.
        else:
            raise Exception("TODO: tomofdpcg - newPokemxSparse - adding invPhaseCov array to GCGsparse")
        self.diagonalA=util.spmatrix.getDiagonal(self.GCGsparse)
        #diagA=A.theDiagonal#used if options["diagonaliseA"] is 1.
        self.sparseA=self.GCGsparse#already sparse!
        #SPDAsparse=(A+AT)/2#symmetric anyway, so this probably isn't necessary.
        #Note, the 2D FFT of a diagonal matrix is the same as a 1D FFT of the diagonal...
        #ie the 1d array is same as the first row of the 2d array, which is then circulant backwards.
        self.ChatBig2=numpy.fft.fft(self.diagonalA)#This is circulant... but backwards.  However, this doesn't matter, because in reorder it would make it circulant (except here we only have the first row anyway!).
        self.reorderChat(self.ChatBig2,circulant=1)
        


    def reorderChat(self,chat,circulant=1):
        """Reorder Chat matrix, so that can be used as banded or diagonal.
        You should check that this has worked for your situation.
        """
        
        print "reordering Chat so that it is banded or diagonal - you should check this has worked!"
        if len(chat.shape)==1:#circulant, so reordered anyway...
            self.ChatReordered=chat
        else:
            self.ChatReordered=numpy.zeros(chat.shape,chat.dtype)
            size=chat.shape[0]
            for i in range(size):
                shift=(i*2)%size
                if shift==0:
                    shift=size
                self.ChatReordered[i,:shift]=chat[i,-shift:]
                self.ChatReordered[i,shift:]=chat[i,:-shift]# The real part should now look pretty much diagonal.  But, not sure what to do about the imag part.
        #Note, the values on the diagonal of ChatReordered seem to equal sum(numpy.diagonal(A)).  So, this could be useful for a quick approximation...
        print "fdpcg.reorderChat - Matrix inverse ChatReordered (shape=%s)"%str(self.ChatReordered.shape)
        if circulant:
            invChatReordered=util.phaseCovariance.invertCirculantMatrix(self.ChatReordered)
        else:
            self.invChatReordered=numpy.linalg.inv(self.ChatReordered)

        if self.options["diagonaliseinvChatReordered"]:#take just the diagonal - which is the same everywhere...
            #self.invChatReordered=numpy.diagonal(self.invChatReordered)
            #self.invChatReordered=numpy.zeros(invChatReordered.shape,chat.dtype)
            self.invChatReordered=invChatReordered[0]#just keep a single value!!!
        else:#expand up the invChatReordered matrix.
            self.invChatReordered=util.phaseCovariance.expandToCirculant(invChatReordered)
        print "fdpcg.reorderChat - done inverse"
        #The invChatReordered matrix is only used in computeinvChatr, and at the moment, this just uses the diagonal...
        #Need to apply the same reordering to any vector multiplied by this.
        #so now use invChatReordered multiplied by FFT of r to give FFT of z...
        #Assume Chat is diagonal?  Not a good approx for imag part.  But does this matter?  Also reasonably poor approx for real part, though not so bad.  Possibly could use as banded.
    def computeinvChatr(self,invChat,r,mode="diag"):
        """Compute the result of C-1r
        This uses the inverse reordered chat matrix, which is created from the FFT of the diagonalised A matrix.
        
        """
        if mode=="diag":
            #create array for FFTs...
            #layerSize=r.shape[0]/self.nLayers
            #if layerSize!=self.nact*self.nact:
            #    print "WARNING: fdpcg - layersize calculated wrongly..."
            #fft the L blocks of r...
            fftr=numpy.empty(r.shape,numpy.complex64)
            layerStart=0
            for i in xrange(self.nLayers):
                #rhat=numpy.fft.fft2(numpy.reshape(r[i*layerSize:(i+1)*layerSize],(self.nact,self.nact))).flat
                ## then reorder rhat.
                #rhatreordered=fftr[i*olayerSize:(i+1)*olayerSize]#numpy.zeros(rhat.shape,rhat.dtype)
                #rhatreordered[0]=rhat[0]
                #rhatreordered[1:]=rhat[-1:0:-1]
                layerEnd=self.dmList[i].nact**2
                rhat=numpy.fft.fft2(numpy.reshape(r[layerStart:layerStart+layerEnd],(self.dmList[i].nact,self.dmList[i].nact))).flat
                # then reorder rhat.
                rhatreordered=fftr[layerStart:layerStart+layerEnd]#numpy.zeros(rhat.shape,rhat.dtype)
                rhatreordered[0]=rhat[0]
                rhatreordered[1:]=rhat[-1:0:-1]
                layerStart+=layerEnd
            #now perform the matrix multiplication... (in fourier space).
            if type(invChat)==type(0.) or len(invChat.shape)==1:#just the diagonal...
                zhatreordered=invChat*fftr#rhatreordered*invChatdiag
            else:
                print "slow mmx multiply"
                zhatreordered=quick.dot(invChat,fftr)#Slow MMX multiply
            # now reorder zhat back to coord sys and FFT back...
            layerStart=0
            for i in xrange(self.nLayers):
                layerEnd=self.dmList[i].nact**2
                z=fftr[layerStart:layerStart+layerEnd]#instead of creating new arr
                zhat=rhat
                zhatr=zhatreordered[layerStart:layerStart+layerEnd]
                zhat[0]=zhatr[0]
                zhat[1:]=zhatr[-1:0:-1]
                z[:]=numpy.fft.ifft2(numpy.reshape(zhat,(self.dmList[i].nact,self.dmList[i].nact))).ravel()
                layerStart+=layerEnd
            z=fftr.real#set reference to start of array...
            rtval=z[:layerStart]
            #print rtval.shape,rtval.dtype.char
##             if self.options["removeHiFreqZ"]:
##                 #attempt to remove the high frequency components of z...
##                 for i in range(self.nLayers):
##                     zfft=numpy.fft.fft2(numpy.reshape(rtval[i*layerSize:(i+1)*layerSize],(self.nact,self.nact)))
##                     f=(self.nact-1)/2
##                     t=(self.nact+1)/2+1
##                     zfft[f:t,f:t]=numpy.average(zfft[f:t,f:t].flat)
##                     rtval[i*layerSize:(i+1)*layerSize]=numpy.fft.ifft2(zfft).ravel().real
            return rtval
##         elif mode=="banded":
##             print "todo: pgc - FFT the L blocks of r separately"
##             rhat=numpy.fft.fft(r)#FFT the L blocks of r.
##             zhat=self.solve_banded(self.Chat,rhat,self.sparsitywidth,iscompact=1)#solve Cz=r (sparse, linear system).
##             z=numpy.fft.ifft(zhat)#inverse FFT the L blocks of zhat.
##             return z.real
    def newCentroids(self,centroids):
        """At start of next simulation iteration - have a new set of
        centroids to reconstruct the actuator values from."""
        #centroids=numpy.array(centroids)
        self.centroids[:]=centroids.astype(self.centroids.dtype.char)
        if self.options["useSparsePokeMx"]:
            #print self.pokeNoiseSparse.mx.shape,centroids.shape,type(numpy.array(centroids))
            self.b[:]=self.pokeNoiseSparse.dot(centroids)
##         else:
##             print "warning: pcg.newCentroids() - not sparse matrix"
##             self.b[:]=quick.dot(self.pokeNoise,centroids)

    def resetX(self):
        """reset the initial guess... can do this when the loop is
        opened or closed.
        Problem: if set to zero, end up dividing by zero in algorithm."""
        self.x[:]=0.
                    
    def initialise(self,usePrevious=1):
        """Initialise the PCG algorithm ready for iterating"""
        #self.x[:]=0.#initial guess?  Or maybe leave as was from last
        # simulation iteration...
        if usePrevious==0:
            self.x[:]=0.
        if self.options["diagonaliseA"]==1:#A has been diagonalised.
            self.r=self.b-self.diagonalA*self.x
        elif self.options["diagonaliseA"]==0.5:#sparse
            self.r=self.b-self.sparseA.matvec(self.x)
            #self.r=self.b.copy()
            #for key in self.sparseA.keys():
            #    if key==0:
            #        self.r-=self.sparseA[key]*self.x
            #    elif key<0:
            #        self.r[-key:]-=self.sparseA[key]*self.x[:key]
            #    else:
            #        self.r[:-key]-=self.sparseA[key]*self.x[key:]
        else:#full implementation
            print "big mmx multiply"
            self.r=self.b-quick.dot(self.A,self.x)
        #self.z=quick.dot(self.Cm1,self.r)
        #Actually, do the previous line in FFT space (next 3 lines).
        #print "todo: pgc - FFT the L blocks of r separately"
        #self.rhat=numpy.fft.fft(self.r)#FFT the L blocks of r.
        #self.zhat=self.solve_banded(self.Chat,self.rhat,self.sparsitywidth,iscompact=1)#solve Cz=r (sparse, linear system).
        #self.z=ifft(self.zhat)#inverse FFT the L blocks of zhat.
        self.z=self.computeinvChatr(self.invChatReordered,self.r,mode="diag").real
        self.p=self.z
        self.convergenceIter=0
    def compute(self):
        """Perform iterations until convergence (note - what is
        convergence, and how do we test for it?  Here, we just do 5
        iterations."""
        converged=0
        self.tolerance[:]=0.
        self.alphaHist[:]=0.
        self.betaHist[:]=0.
        while self.convergenceIter<self.minIter or (converged==0 and self.convergenceIter<self.maxIter):
            self.nextIter()
            self.convVal=max(numpy.fabs(self.xprev-self.x))
            self.tolerance[self.convergenceIter]=self.convVal
            self.alphaHist[self.convergenceIter]=self.alpha
            self.betaHist[self.convergenceIter]=self.beta
            #print max(numpy.fabs(self.xprev-self.x))
            if self.convVal<self.convergenceValue:
                converged=1
            else:
                converged=0
            print self.convVal
            self.convergenceIter+=1
        if converged==0 and self.convergenceWarning==1:
            print "Warning, PCG algorithm not yet converged.  Continuing anyway. (convVal=%g, requested %g)"%(self.convVal,self.convergenceValue)
        #elif converged==1:
        #    print "PCG converged after %d iters to %g"%(self.convergenceIter,self.convVal)
    def nextIter(self):
        """Perform a PCG iteration"""
        self.p=self.p.real
        self.xprev=self.x.copy()
        if self.options["diagonaliseA"]==1:#A has been diagonalised.
            self.q=self.diagonalA*self.p
        elif self.options["diagonaliseA"]==0.5:#sparse...
            self.q=self.sparseA.matvec(self.p)
            #self.q=numpy.zeros(self.p.shape,self.p.dtype)
            #for key in self.sparseA.keys():
            #    if key==0:
            #        self.q+=self.sparseA[key]*self.p
            #    elif key<0:
            #        self.q[-key:]+=self.sparseA[key]*self.p[:key]
            #    else:
            #        self.q[:-key]+=self.sparseA[key]*self.p[key:]
            
        else:#full implementation
            print "big mmx multiply for q=Ap"
            self.q=quick.dot(self.A,self.p)#dominant cost
        qp=quick.dot(self.q,self.p)#complex
        self.rz=quick.dot(self.r,self.z)
        if qp!=0:
            #print type(self.r),type(self.z),type(qp)#agbhome
            #print self.r.shape,self.z.shape,self.r.dtype.char,self.z.dtype.char
            #tmp=quick.dot(self.r,self.z)
            #print type(tmp)
            self.alpha=(self.rz/qp).real#complex
        else:
            self.alpha=0.#if q or p are zero, implies z must be, so r must be.
        self.x+=(self.alpha*self.p).real#complex, but seems that x.imag tends to zero???  Can we ignore it right from the start?
##         if self.options["removeHiFreqX"]:
##             lsize=self.nact*self.nact
##             for i in range(self.nLayers):
##                 fftx=numpy.fft.fft2(numpy.reshape(self.x[i*lsize:(i+1)*lsize],(self.nact,self.nact)))
##             f=(self.nact-1)/2
##             t=(self.nact+1)/2+1
##             fftx[f:t,f:t]=numpy.average(fftx[f:t,f:t].flat)
##             self.x[i*lsize:(i+1)*lsize]=numpy.fft.ifft2(fftx).ravel().real
            
        
        self.r-=self.alpha*self.q
        #self.z=quick.dot(self.Cm1,self.r)#dominant cost
        #Prev line actually computed in FFT space (next line...).
        self.z=self.computeinvChatr(self.invChatReordered,self.r,mode="diag").real
        if self.rz!=0:
            self.beta=(quick.dot(self.r,self.z)/self.rz).real
        else:
            self.beta=0.#if rz is zero, so is new r or z...
        #print "alpha,beta",self.alpha,self.beta
        self.p=(self.z+self.beta*self.p).real
##         if self.options["removeHiFreqP"]:
##             lsize=self.nact*self.nact
##             for i in range(self.nLayers):
##                 fftp=numpy.fft.fft2(numpy.reshape(self.p[i*lsize:(i+1)*lsize],(self.nact,self.nact)))
##             f=(self.nact-1)/2
##             t=(self.nact+1)/2+1
##             fftp[f:t,f:t]=numpy.average(fftp[f:t,f:t].flat)
##             self.p[i*lsize:(i+1)*lsize]=numpy.fft.ifft2(fftp).ravel().real
    def finalise(self):
        #self.estPhase=numpy.reshape(self.x*self.gainfactor,(self.nact,self.nact))
        #av=numpy.average(self.x)
        #print "Average pcg.x: %g - removing it."%av
        #self.x-=av
        #estphase gets updated automatically, so no need to do anything here (provided self.x is the same array).
        pass
    def solve(self,centroids,usePrevious=1):
        """Get an estimate for the phase...
        This doesn't reset the x vector - assumes that best guess for new
        phase is the previous phase (from last simulation iteration).
        
        """
        self.newCentroids(centroids)
        self.initialise(usePrevious=usePrevious)
        self.compute()
        self.finalise()

        




    


def createPcg(pokemx,noisecov=0.25,nact=9,telDiam=4.2,r0=0.2,l0=30.):
    """Probably, this can be used from the simulation...
    pokemx is the poke matrix.
    noisecov is the noise covariance value (can also be a 1d or 2d (diagonal) array).
    nact is number of actuators (assume fried geometry) in 1 direction.
    """

    bccb=numpy.array(util.phaseCovariance.bccb(nact=nact,telDiam=telDiam,r0=r0,l0=l0))
    phasecovVector=1./numpy.diagonal(numpy.linalg.inv(bccb))
    Pcg=pcg(pokemx,noisecov,phasecovVector,None)
    return Pcg
    
def phase2acts(phase):
    """Expand phase array by 1 in each direction onto actuators..."""
    acts=numpy.zeros((phase.shape[0]+1,phase.shape[1]+1),"d")
    for i in range(acts.shape[0]):
        for j in range(acts.shape[1]):
            n=0
            if i>0:
                if j>0:
                    acts[i,j]+=phase[i-1,j-1]
                    n+=1
                if j<phase.shape[1]:
                    acts[i,j]+=phase[i-1,j]
                    n+=1
            if i<phase.shape[0]:
                if j>0:
                    acts[i,j]+=phase[i,j-1]
                    n+=1
                if j<phase.shape[1]:
                    acts[i,j]+=phase[i,j]
                    n+=1
            if n!=0:
                acts[i,j]/=n
    return acts

def getFocusCents(zernlist={3:1.},nact=9,readnoise=10.,usePoisson=1,sig=500./64,plot=0,fullpupil=0):
    """Compute centroids that would be obtained for a focus zernike term.
    zernlist can also be a phasescreen...
    """
    import util.zernikeMod,util.centroid,util.tel
    nphs=(nact-1)*8
    pupil=util.tel.Pupil(nphs,nphs/2,0)
    if fullpupil:
        pupil.fn[:]=1
    avphase=numpy.zeros((nact-1,nact-1),numpy.float64)
    if type(zernlist)==numpy.ndarray:
        phase=zernlist
    else:
        zern=util.zernikeMod.Zernike(pupil,max(zernlist.keys())+1)
        phase=numpy.zeros((nphs,nphs),"d")
        for key in zernlist.keys():
            phase+=zernlist[key]*zern.zern[key]#focus zernike
    for i in range(nact-1):
        for j in range(nact-1):
            avphase[i,j]=numpy.average(numpy.array(phase[i*8:i*8+8,j*8:j*8+8].flat))
    c=util.centroid.centroid(nact-1,pupil.fn,readnoise=readnoise,addPoisson=usePoisson,readbg=0.,sig=sig)
    c.printmax=1
    cents=c.calc(phase)
    if plot:
        print "Max counts: %g, max, min per subap %g %g, plotting SHS image..."%(max(c.tile.flat),max(c.photPerSubap.flat),min(c.photPerSubap.flat))
        gist.fma();gist.pli(c.tile)
    #cents=util.centroid.calc(phase,8,pupil.fn)
    return c,phase,avphase

def randPhase(nact=9,zmax=55):
    """Create a random phase using random zernikes.
    """
    import util.zernikeMod
    nsubx=nact-1
    zerns=util.zernikeMod.Zernike(nsubx*8,zmax+1).zern
    weights=numpy.random.random(zmax)*2-1
    phase=numpy.zeros((nsubx*8,nsubx*8),numpy.float64)
    for i in range(3,zmax+1):
        phase+=zerns[i]*weights[i-1]
        #gist.fma();gist.pli(phase)
        #raw_input()
    return phase
def chi2(arr1,arr2,scale=0):
    """compute chi squared"""
    arr1=arr1.copy().ravel()
    arr2=arr2.copy().ravel()
    if scale:
        arr1-=min(arr1)
        arr2-=min(arr2)
        arr1/=max(arr1)
        arr2/=max(arr2)
    chi=numpy.sum((arr1-arr2)**2)
    return float(chi)

def doit(zernlist=None,nact=9,cents=None,avphase=None,readnoise=10.,usePoisson=1,sig=1000.,fullpupil=0,monteNoiseCovariance=0,phaseCov="bccb",diagonaliseinvChatReordered=1,diagonaliseA=0.5,useChatFromDiagonalisedA=1,oversampleFFT=0,fft2d=1,oversampleAndBinFFTofr=0,removeHiFreqZ=0,removeHiFreqP=0,removeHiFreqX=0,convToPxl=1):
    """Create a zernike mode and get the centroids for this.  Put these centroids into the PCG algorithm, and recreate the phase, and then compare with the original input phase
    defaultOptions are:
    diagonaliseinvChatReordered:1#should we diagonlise C-1, or use whole MX?
    diagonaliseA:0#should we diagonalise A or use the whole MX.
    useChatFromDiagonalisedA:1#should we compute invChat using a diagonalised chat or not?
    zernlist can be a dict of zernikes eg {3:1} would be focus, amplitude 1...
    or it can be the phase, or it can be none (in which case a phasescreen is
    created).

    This now seems to be working fairly well, iterates to an exact solution after about 10 iters, and 4-5 should be sufficient for a good solution.

    """
    options={"diagonaliseinvChatReordered":diagonaliseinvChatReordered,
             "diagonaliseA":diagonaliseA,
             "useChatFromDiagonalisedA":useChatFromDiagonalisedA,
             "oversampleFFT":oversampleFFT,
             "2dfft":fft2d,
             "oversampleAndBinFFTofr":oversampleAndBinFFTofr,
             "removeHiFreqZ":removeHiFreqZ,
             "removeHiFreqP":removeHiFreqP,
             "removeHiFreqX":removeHiFreqX
             }
    gist.window(0)
    gist.palette("gray.gp")
    pupfn=util.tel.Pupil(nact-1,(nact-1)/2.,0).fn
    actfn=util.tel.Pupil(nact,nact/2.,0).fn
    noiseCov=None
    if type(cents)==type(None):
        if type(zernlist) in [type(1),type(1.)]:
            zernlist={zernlist:1.}
        elif type(zernlist)==type([]):
            tmp={}
            for i in zernlist:
                tmp[i]=1.
            zernlist=tmp
        elif type(zernlist)==type(None):
            zernlist=science.infScrn.makeInitialScreen(dpix=(nact-1)*8,Dtel=4.2,L0=30.,scrnXPxls=None,scrnYPxls=None,seed=None,tstep=0.05,globR0=0.2,strLayer=1.,natype=numpy.float64,windDirection=0.,vWind=10.)[:-1,1:].copy()
        c,phase,avphase=getFocusCents(zernlist,nact=nact,readnoise=readnoise,usePoisson=usePoisson,sig=sig,plot=1,fullpupil=fullpupil)
        print phase.shape
        cents=c.cent
        if monteNoiseCovariance:
            noiseCov=c.computeNoiseCovariance(20,convertToRad=1)
    if convToPxl:#convert from radians to pixel values for the noises...
        radToPxlFactor=1./c.convFactor**2
    else:
        radToPxlFactor=1.
    Pcg=testfull(fromdisk=0,nact=nact,pupfn=pupfn,actfn=actfn,fullpmx=1,noiseCov=noiseCov,phaseCov=phaseCov,options=options,radToPxlFactor=radToPxlFactor)#agbhome
    Pcg.newCentroids(cents)
    Pcg.initialise()

    #now inverse A, multiply with b, to get the MVM reconstructed phase...
    invA=numpy.linalg.inv(Pcg.A)
    gist.fma();gist.pli(invA)
    raw_input("Displaying inverse of A... press return")
    gist.fma();gist.pli(phase)
    raw_input("The phase... press a key")
    recphase=quick.dot(invA,Pcg.b)
    print "Reconstructed phase min/max:",min(recphase.flat),max(recphase.flat)
    recphase.shape=(9,9)
    gist.window(4);gist.palette("gray.gp");gist.fma();gist.pli(recphase)
    gist.palette("gray.gp")
    gist.window(0)
    chires=numpy.zeros((100,),"d")
    #also compute what a traditional MVM with pokemx would give...
    #invpokemx=numpy.linalg.pinv(Pcg.pokemx)#same as generalised_inverse
    #pmphase=quick.dot(invpokemx,cents)
    print "Press return for next iteration or key+return to quit"
    #gist.fma()
    #gist.pli(numpy.reshape(pmphase,(nact,nact)))
    if type(avphase)!=type(None):
        print "press a key"
        raw_input()
        gist.fma()
        gist.pli(avphase)
    actphase=phase2acts(avphase)
    actphase-=numpy.average(actphase.flat)#remove piston
    smallpupfn=util.tel.Pupil(nact,nact/2.-2,0).fn
    raw_input("press return (min, max actphase is %g %g"%(min(actphase.flat),max(actphase.flat)))
    gist.fma();gist.pli(actphase)
    gist.window(1)
    gist.palette("gray.gp")
    gist.window(2)
    gist.palette("gray.gp")
    niter=0
    while len(raw_input())==0:
        niter+=1
        gist.window(1);gist.fma()
        img=numpy.reshape(Pcg.x.real,(nact,nact))
        gist.pli(img)
        gist.window(2);gist.fma()
        img=smallpupfn*img
        #gist.pli(numpy.where(img==0,min(img.flat),img))
        gist.pli(numpy.reshape(Pcg.p,(nact,nact)))
        gist.window(3);gist.fma()
        fftimg=numpy.fft.fft2(numpy.reshape(Pcg.x,(nact,nact)))
        gist.pli(fftimg.real)
        

        Pcg.nextIter()
        chirespos=niter-1
        if chirespos>99:
            chirespos=99
        chires[chirespos]=chi2(Pcg.x,recphase,scale=0)#changed from actphase
        Pcg.alphaHist[chirespos]=Pcg.alpha
        Pcg.betaHist[chirespos]=Pcg.beta
        Pcg.tolerance[chirespos]=max(numpy.fabs(Pcg.xprev-Pcg.x))
        print niter,Pcg.tolerance[chirespos],chires[chirespos],min(Pcg.x.real),max(Pcg.x.real),min(Pcg.p),max(Pcg.p),Pcg.alphaHist[chirespos],Pcg.betaHist[chirespos]
        print "Press return for next iteration or key+return to quit"

    gist.fma();gist.plg(chires[:chirespos])
    gist.window(2);gist.fma();gist.plg(Pcg.tolerance[:niter])
    gist.window(0)
    return Pcg,actphase,phase,recphase,c
#actuators=util.FITS.Read("/home/ali/cvsstuff/aosim/example/testinterp/actuators.fits")[1]
#simcents=util.FITS.Read("/home/ali/cvsstuff/aosim/example/testinterp/centroidsMixed.fits")[1]
#simcents2=Numeric.zeros((128,),simcents.typecode())
#simcents2[:64]=Numeric.reshape(simcents[:,:,0],(64,))
#simcents2[64:]=Numeric.reshape(simcents[:,:,1],(64,))
#Pcg=pcg.doit(9,simcents2,Numeric.reshape(actuators,(9,9))) 
    
