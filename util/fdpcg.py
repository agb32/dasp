#A module for pre-conditioned conjugate gradient algorithms.
import numpy
import numpy.fft as fft
import types
import util.FITS,util.phaseCovariance
import util.dot as quick
try:
    import gist
except:
    pass
#import util.createPokeMx
try:
    import science.infScrn
except:
    pass
"""Aims to solve Ax=b.

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

This code currently isn't fully tested, and is a bit scrappy.  When
using as part of the simulation, you are probably best to call the
createPcg method initially, use pcg.newPokemx to add a new poke
matrix, and use pcg.solve each simulation iteration.
Actually, no, use the xinterp_recon object.

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
    def __init__(self,pokemx, noiseCov, phaseCov,centroids,nact=None,sparsityWidth=0.25,turbStrength=numpy.ones((1,),numpy.float32),minIter=10,maxIter=100,convergenceValue=5e-4,options=None,gainfactor=1.,convergenceWarning=1,fakePokeMxVal=0.31):
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
         - gainfactor - a factor used to multiply the phase at the end of the
           solve, to get it closer P-V to the true value.
        """
        self.options=options#a list of options for what is diagonalised etc - ie tuning the algorithm.
        if self.options==None:
            self.options={}
        defaultOptions={"diagonaliseinvChatReordered":1,#should we diagonalise C-1, or use the whole MX?
                        "diagonaliseA":0.5,#should we diagonalise A or use the whole MX.
                        "useChatFromDiagonalisedA":1,#should we compute invChat using a diagonalised chat or not?
                        "oversampleFFT":0,#over sample the fourier step
                        "2dfft":1,#do a 2D FFT of x to get z instead of 1d.
                        "oversampleAndBinFFTofr":0,#whether to oversample fft in some intermediate points.  This can be 0, 1, 2 or 3 (see code for details)
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
        if type(turbStrength)==type([]):
            self.nLayers=len(turbStrength)
        elif type(turbStrength) in [type(1),type(1.)]:
            self.nLayers=1
        else:
            self.nLayers=turbStrength.shape[0]
        self.convergenceWarning=convergenceWarning
        self.noiseCov=noiseCov
        if nact==None:
            self.nact=int(numpy.sqrt(pokemx.shape[1]))
            if self.nact*self.nact!=pokemx.shape[1]:
                raise Exception("Number of actuators not supplied, and pokemx appears to give a non-square shape")
        else:
            self.nact=nact
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
        if type(pokemx)==type(None):
            self.newPokemxSparse()
        else:
            self.newPokemx(pokemx)

        #Depending on what you use for the invPhaseCov matrix (ie the BCCB approximation), you can possibly diagonalise A which then makes the A matrix multiplication order N rather than N^2.  The BCCB approximation just assumes that the inverse phase covariance matrix is diagonal, which though not perfect, may be reasonable.  Actually, typically, the invPhaseCov added to GCG is fairly small, so doesn't have much effect - so usually, if GCG is approx diagonal, A can be diagonalised.
        if type(centroids)==types.NoneType:
            if type(self.pokemx)==type(None):#will be using sparse representation...
                centroids=numpy.zeros(((self.nact-1)**2*2,),numpy.float32)
            else:
                centroids=numpy.zeros((pokemx.shape[0],),numpy.float32)
        else:
            centroids=numpy.array(centroids)
        self.centroids=centroids
        #self.b=quick.dot(self.pokeNoise,centroids)
        self.b=numpy.zeros((self.nact**2),self.centroids.dtype)
        print "computing pcg centroids"
        self.newCentroids(self.centroids)#compute self.b

        #Here, assume poke matrix has shape ncents,nmodes.
        #self.b=b#eg the pokemxT times noise covariance-1 times centroid measurements (G^T C^-1 s), shape nmodes.  Pokemx has shape ncents,nmodes.
        #self.A=A#eg the poke matrix transposed times inverse noise covariance matrix times poke matrix plus Noll matrix, final shape (nmodes, nmodes).
        self.x=numpy.zeros((self.nact*self.nact,),self.b.dtype)#initial guess, shape nmodes.
        self.estPhase=self.x.view()
        self.estPhase.shape=(self.nact,self.nact)
##         self.l0=l0
##         self.kappa=computeTurbSpectrum(nmodes/nLayers,self.actPosition,l0=l0,innerScale=innerScale)
##         if len(self.invPhaseCov.shape)==1:#diagonal
##             invPhaseCov=numpy.zeros((self.invPhaseCov.shape[0],self.invPhaseCov.shape[0]),self.invPhaseCov.dtype)
##             for i in range(self.invPhaseCov.shape[0]):
##                 invPhaseCov[i,i]=self.invPhaseCov[i]
##         else:
##             invPhaseCov=self.invPhaseCov
##         self.kappa=numpy.fft.inverse_fft2(invPhaseCov)
        #self.ChatBig=computeChatFromPokeMatrix2(pokemx,self.sigma,turbStrength,self.kappa)
        #The question then is whether this ChatBig can be used... ie is it correct?  Then, if so, how does one reorder it, if at all?
#        if type(sparsityWidth)==types.FloatType:
#            sparsityWidth=int(self.ChatBig.shape[0]*sparsityWidth)
#        self.sparsityWidth=sparsityWidth
#        self.Chat=compactToBanded(self.ChatBig,sparsityWidth/2,sparsityWidth/2)
        self.minIter=minIter
        self.maxIter=maxIter
        self.tolerance=numpy.zeros((maxIter,),numpy.float32)
        self.alphaHist=numpy.zeros((maxIter,),numpy.float32)
        self.betaHist=numpy.zeros((maxIter,),numpy.float32)
        self.convergenceValue=convergenceValue
        self.gainfactor=gainfactor
        self.convVal=0.
        print "finished pcg initialisation"
    def newPokemx(self,pokemx):
        """Pass a new poke matrix in..."""
        pokemx=numpy.array(pokemx)
        self.pokemx=pokemx
        pokemxT=numpy.transpose(pokemx)
        #pn=quick.dot(pokemxT,invNoiseCovVector
        #multiply pokemxT with invNoiseCovVector...
        pn=numpy.zeros(pokemxT.shape,pokemx.dtype.char)
        if type(self.invNoiseCovVector)==type(0.):
            pn[:,:,]=pokemxT*self.invNoiseCovVector
        else:
            for i in range(pn.shape[1]):
                pn[:,i]=pokemxT[:,i]*self.invNoiseCovVector[i]
        self.pokeNoise=pn
        self.pokeNoiseSparse=util.createPokeMx.sparsePMX(self.nact,pn)
        if self.options["useSparsePokeMx"]:
            self.GCG=self.pokeNoiseSparse.dot(pokemx)
        else:
            print "warning - pcg: not using sparse version of GCG creation"
            self.GCG=quick.dot(pn,pokemx)#with sparse approx, speed this up.
        if len(self.invPhaseCov.shape)==1:#diagonal...
            self.A=self.GCG.copy()
            #print self.A.shape,self.invPhaseCov.shape
            #for i in range(self.A.shape[0]):
            #    self.A[i,i]+=self.invPhaseCov[i]
            self.A.flat[::self.A.shape[0]+1]+=self.invPhaseCov#add to diagonal.
        else:
            self.A=self.GCG+self.invPhaseCov
        self.diagonalA=numpy.diagonal(self.A)#this is used if options["diagonaliseA"] is set.
        #here, we treat the sparse matrix as the diagonal, one column each side, and the 3 columns around each semi-diagonal... (look at self.A and you'll see what I mean).
        self.sparseA={0:numpy.diagonal(self.A)}
        self.sparseA[1]=numpy.diagonal(self.A[:,1:,])
        self.sparseA[-1]=numpy.diagonal(self.A[1:,])
        self.sparseA[self.nact]=numpy.diagonal(self.A[:,self.nact:,])
        self.sparseA[-self.nact]=numpy.diagonal(self.A[self.nact:,])
        self.sparseA[self.nact-1]=numpy.diagonal(self.A[:,self.nact-1:,])
        self.sparseA[-(self.nact-1)]=numpy.diagonal(self.A[self.nact-1:,])
        self.sparseA[self.nact+1]=numpy.diagonal(self.A[:,self.nact+1:,])
        self.sparseA[-(self.nact+1)]=numpy.diagonal(self.A[self.nact+1:,])

        
        self.SPDA=(self.A+numpy.transpose(self.A))/2#symmetric positive definite version of A (at least, symmetric, and PD in at least some situations...) - invact, virtually identical to A...
        if self.options["oversampleFFT"]:
            if self.options["2dfft"]:
                SPDA=numpy.zeros(map(lambda x:4*x,self.SPDA.shape),self.SPDA.dtype)
                SPDA[:self.SPDA.shape[0],:self.SPDA.shape[1]]=self.SPDA
            else:
                SPDA=numpy.zeros(map(lambda x:2*x,self.SPDA.shape),self.SPDA.dtype)
                SPDA[:self.SPDA.shape[0],:self.SPDA.shape[1]]=self.SPDA
            self.SPDA=SPDA
        print "fdpcg: Creating ChatBig"
        self.ChatBig=numpy.fft.fft2(self.SPDA)
        self.Adiag2d=numpy.zeros(self.A.shape,self.A.dtype)#A diagonal matrix
        #copy only the diagonal values from A to Adiag2d
        self.Adiag2d.flat[::self.A.shape[0]+1]=self.A.flat[::self.A.shape[0]+1]
        if self.options["oversampleFFT"]:
            if self.options["2dfft"]:
                Adiag2d=numpy.zeros(map(lambda x:4*x,self.A.shape),self.A.dtype)
                Adiag2d[:self.Adiag2d.shape[0],:self.Adiag2d.shape[1]]=self.Adiag2d
            else:
                Adiag2d=numpy.zeros(map(lambda x:2*x,self.A.shape),self.A.dtype)
                Adiag2d[:self.Adiag2d.shape[0],:self.Adiag2d.shape[1]]=self.Adiag2d
            self.Adiag2d=Adiag2d
        print "fdpcg: Creating ChatBig2"
        self.ChatBig2=numpy.fft.fft2(self.Adiag2d)
        if self.options["useChatFromDiagonalisedA"]:
            self.reorderChat(self.ChatBig2,circulant=1)
        else:
            self.reorderChat(self.ChatBig,circulant=0)


    def newPokemxSparse(self):
        """Pass a new poke matrix in... in sparse format... (so that larger problems can fit in memory!)
        This basically does the same as newPokemx for some specific allowed cases (should raise an error
        if not allowed), but only ever uses sparse matricees, to avoid having to create the full matrix
        (enough memory may not be available!).
        """
        self.pokemx=None
        nact=self.nact
        for key in self.options.keys():
            if self.defaultOptions[key]!=self.options[key]:
                raise Exception("fdpcg - newPokemxSparse - must use default options...")
        spmx=util.createPokeMx.sparsePMX(nact)
        spmx.create(self.fakePokeMxVal)
        spmxT=util.createPokeMx.sparsePMX(nact)
        spmxT.create(self.fakePokeMxVal)
        if type(self.invNoiseCovVector)==type(0.):
            spmx.storage*=self.invNoiseCovVector
        else:#size probably ncents...
            nsubx=int(numpy.sqrt(self.invNoiseCovVector.shape[0]/2))
            spmx.storage[:4]*=numpy.reshape(self.invNoiseCovVector[:nsubx*nsubx],(nsubx,nsubx))
            spmx.storage[4:,]*=numpy.reshape(self.invNoiseCovVector[nsubx*nsubx:,],(nsubx,nsubx))
            #raise Exception("TODO: fdpcg - newPokemxSparse - invNoiseCovVector")
        self.pokeNoiseSparse=spmx
        self.GCGsparse=self.pokeNoiseSparse.dotWithSparsePmx(spmxT)#returns a dictionary of the diagonal elements that are required...
        if type(self.invPhaseCov)==type(0.) or len(self.invPhaseCov.shape)==1:#diagonal...
            self.GCGsparse[0]+=self.invPhaseCov#add to the diagonal element.
            #A=GCGsparse.copy()
            #A.theDiagonal+=self.invPhaseCov#add to diagonal.
        else:
            #A=GCGSparse+self.invPhaseCov
            raise Exception("TODO: fdpcg - newPokexxSparse - adding invPhaseCov array to GCGsparse")
        self.diagonalA=self.GCGsparse[0]
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
                self.ChatReordered[i,:shift]=chat[i,-shift:,]
                self.ChatReordered[i,shift:,]=chat[i,:-shift]# The real part should now look pretty much diagonal.  But, not sure what to do about the imag part.
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
            layerSize=r.shape[0]/self.nLayers
            if layerSize!=self.nact*self.nact:
                print "WARNING: fdpcg - layersize calculated wrongly..."
            #fft the L blocks of r...
            #print "todo: pgc - FFT the L blocks of r separately"
            if self.options["oversampleFFT"]:
                if self.options["2dfft"]:
                    tmpr=numpy.zeros((self.nact*2,self.nact*2,),r.dtype)
                    fftr=numpy.empty((r.shape[0]*4,),numpy.complex64)
                    olayerSize=layerSize*4#oversampled layer size
                else:
                    tmpr=numpy.zeros((layerSize*2,),r.dtype)
                    fftr=numpy.empty((r.shape[0]*2,),numpy.complex64)
                    olayerSize=layerSize*2#the oversampled layer size...
            else:
                fftr=numpy.empty(r.shape,numpy.complex64)
                olayerSize=layerSize
            for i in range(self.nLayers):
                if self.options["oversampleFFT"]:
                    if self.options["2dfft"]:
                        tmpr[:self.nact,:self.nact]=numpy.reshape(r[i*layerSize:(i+1)*layerSize],(self.nact,self.nact))
                        rhat=numpy.fft.fft2(tmpr).flat
                    else:
                        tmpr[:layerSize]=r[i*layerSize:(i+1)*layerSize]
                        rhat=numpy.fft.fft(tmpr)
                else:
                    if self.options["2dfft"]:
                        if self.options["oversampleAndBinFFTofr"]&1:
                            tmpfft=numpy.zeros((self.nact*2,self.nact*2),r.dtype)
                            tmpfft[:self.nact,:self.nact]=numpy.reshape(r[i*layerSize:(i+1)*layerSize],(self.nact,self.nact))
                            tmpfft=numpy.fft.fft2(tmpfft)
                            rhatbinned=numpy.zeros((self.nact,self.nact),tmpfft.dtype)
                            for j in range(self.nact*2):
                                for k in range(self.nact*2):
                                    rhatbinned[j/2,k/2]+=tmpfft[j,k]
                            rhat=rhatbinned.ravel()
                        else:
                            rhat=numpy.fft.fft2(numpy.reshape(r[i*layerSize:(i+1)*layerSize],(self.nact,self.nact))).flat
                    else:
                        rhat=numpy.fft.fft(r[i*layerSize:(i+1)*layerSize])
                # then reorder rhat.
                rhatreordered=fftr[i*olayerSize:(i+1)*olayerSize]#numpy.zeros(rhat.shape,rhat.dtype)
                rhatreordered[0]=rhat[0]
                rhatreordered[1:,]=rhat[-1:0:-1]
            #now perform the matrix multiplication... (in fourier space).
            if type(invChat)==type(0.) or len(invChat.shape)==1:#just the diagonal...
                zhatreordered=invChat*fftr#rhatreordered*invChatdiag
            else:
                print "slow mmx multiply"
                zhatreordered=quick.dot(invChat,fftr)#Slow MMX multiply
            # now reorder zhat back to coord sys and FFT back...
            for i in range(self.nLayers):
                z=fftr[i*olayerSize:(i+1)*olayerSize]#instead of creating new arr
                zhat=rhat
                zhatr=zhatreordered[i*olayerSize:(i+1)*olayerSize]
                zhat[0]=zhatr[0]
                zhat[1:,]=zhatr[-1:0:-1]
                if self.options["oversampleFFT"]:
                    if self.options["2dfft"]:
                        z[:,]=numpy.fft.ifft2(numpy.reshape(zhat,(self.nact*2,self.nact*2))).ravel()
                    else:
                        z[:,]=numpy.fft.ifft(zhat)
                else:
                    if self.options["2dfft"]:
                        if self.options["oversampleAndBinFFTofr"]&2:
                            tmpfft=numpy.zeros((self.nact*2,self.nact*2),z.dtype)
                            tmpfft[:self.nact,:self.nact]=numpy.reshape(zhat,(self.nact,self.nact))
                            tmpfft=numpy.fft.ifft2(tmpfft)
                            zhatbinned=numpy.zeros((self.nact,self.nact),tmpfft.dtype)
                            for j in range(self.nact*2):
                                for k in range(self.nact*2):
                                    zhatbinned[j/2,k/2]=tmpfft[j,k]
                            z[:,]=zhatbinned.ravel()
                        else:
                            z[:,]=numpy.fft.ifft2(numpy.reshape(zhat,(self.nact,self.nact))).ravel()
                    else:
                        z[:,]=numpy.fft.ifft(zhat)
            z=fftr.real#set reference to start of array...
            if self.options["oversampleFFT"]:
                if self.options["2dfft"]:
                    #gist.fma();gist.pli(numpy.reshape(z,(2*self.nact,2*self.nact)))
                    #raw_input("displaying z")
                    z=numpy.array(numpy.reshape(z,(2*self.nact,2*self.nact))[:self.nact,:self.nact]).ravel()
                else:
                    #gist.fma();gist.pli(numpy.reshape(z[:layerSize],(self.nact,self.nact)))
                    #raw_input("displaying z")
                    #gist.fma();gist.pli(numpy.reshape(z[layerSize:,],(self.nact,self.nact)))
                    #raw_input("displaying z")
                    # now bin z...
                    for i in range(olayerSize):
                        if i*2==0:
                            z[i/2]=z[i]
                        else:
                            z[i/2]+=z[i]
                    #gist.fma();gist.pli(numpy.reshape(z[:layerSize],(self.nact,self.nact)))
                    #raw_input("displaying z")
            else:
                if self.options["2dfft"]:
                        
                    #gist.fma();gist.pli(numpy.reshape(z[:layerSize],(self.nact,self.nact)))
                    #raw_input("displaying z")
                    pass
            rtval=z[:layerSize]
            #print rtval.shape,rtval.dtype.char
            if self.options["removeHiFreqZ"]:
                #attempt to remove the high frequency components of z...
                for i in range(self.nLayers):
                    zfft=numpy.fft.fft2(numpy.reshape(rtval[i*layerSize:(i+1)*layerSize],(self.nact,self.nact)))
                    f=(self.nact-1)/2
                    t=(self.nact+1)/2+1
                    zfft[f:t,f:t]=numpy.average(zfft[f:t,f:t].flat)
                    rtval[i*layerSize:(i+1)*layerSize]=numpy.fft.ifft2(zfft).ravel().real
            return rtval
        elif mode=="banded":
            print "todo: pgc - FFT the L blocks of r separately"
            rhat=numpy.fft.fft(r)#FFT the L blocks of r.
            zhat=self.solve_banded(self.Chat,rhat,self.sparsitywidth,iscompact=1)#solve Cz=r (sparse, linear system).
            z=numpy.fft.ifft(zhat)#inverse FFT the L blocks of zhat.
            return z.real
    def newCentroids(self,centroids):
        """At start of next simulation iteration - have a new set of
        centroids to reconstruct the actuator values from."""
        #centroids=numpy.array(centroids)
        self.centroids[:,]=centroids.astype(self.centroids.dtype.char)
        if self.options["useSparsePokeMx"]:
            self.b[:,]=self.pokeNoiseSparse.dot(centroids)
        else:
            print "warning: pcg.newCentroids() - not sparse matrix"
            self.b[:,]=quick.dot(self.pokeNoise,centroids)

    def resetX(self):
        """reset the initial guess... can do this when the loop is
        opened or closed.
        Problem: if set to zero, end up dividing by zero in algorithm."""
        self.x[:,]=0.
                    
    def initialise(self,usePrevious=1):
        """Initialise the PCG algorithm ready for iterating"""
        #self.x[:,]=0.#initial guess?  Or maybe leave as was from last
        # simulation iteration...
        if usePrevious==0:
            self.x[:,]=0.
        if self.options["diagonaliseA"]==1:#A has been diagonalised.
            self.r=self.b-self.diagonalA*self.x
        elif self.options["diagonaliseA"]==0.5:#sparse
            self.r=self.b.copy()
            for key in self.sparseA.keys():
                if key==0:
                    self.r-=self.sparseA[key]*self.x
                elif key<0:
                    self.r[-key:,]-=self.sparseA[key]*self.x[:key]
                else:
                    self.r[:-key]-=self.sparseA[key]*self.x[key:,]
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
        self.tolerance[:,]=0.
        self.alphaHist[:,]=0.
        self.betaHist[:,]=0.
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
            self.convergenceIter+=1
        if converged==0 and self.convergenceWarning==1:
            print "Warning, PCG algorithm not yet converged.  Continuing anyway. (convVal=%g, requested %g)"%(self.convVal,self.convergenceValue)
        #elif converged==1:
        #    print "PCG converged after %d iters to %g"%(self.convergenceIter,self.convVal)
    def nextIter(self):
        """Perform a PCG iteration"""
        self.p=self.p.real
        self.xprev=self.x
        if self.options["diagonaliseA"]==1:#A has been diagonalised.
            self.q=self.diagonalA*self.p
        elif self.options["diagonaliseA"]==0.5:#sparse...
            self.q=numpy.zeros(self.p.shape,self.p.dtype)
            for key in self.sparseA.keys():
                if key==0:
                    self.q+=self.sparseA[key]*self.p
                elif key<0:
                    self.q[-key:,]+=self.sparseA[key]*self.p[:key]
                else:
                    self.q[:-key]+=self.sparseA[key]*self.p[key:,]
            
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
        self.x=(self.x+self.alpha*self.p).real#complex, but seems that x.imag tends to zero???  Can we ignore it right from the start?
        if self.options["removeHiFreqX"]:
            lsize=self.nact*self.nact
            for i in range(self.nLayers):
                fftx=numpy.fft.fft2(numpy.reshape(self.x[i*lsize:(i+1)*lsize],(self.nact,self.nact)))
            f=(self.nact-1)/2
            t=(self.nact+1)/2+1
            fftx[f:t,f:t]=numpy.average(fftx[f:t,f:t].flat)
            self.x[i*lsize:(i+1)*lsize]=numpy.fft.ifft2(fftx).ravel().real
            
        
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
        if self.options["removeHiFreqP"]:
            lsize=self.nact*self.nact
            for i in range(self.nLayers):
                fftp=numpy.fft.fft2(numpy.reshape(self.p[i*lsize:(i+1)*lsize],(self.nact,self.nact)))
            f=(self.nact-1)/2
            t=(self.nact+1)/2+1
            fftp[f:t,f:t]=numpy.average(fftp[f:t,f:t].flat)
            self.p[i*lsize:(i+1)*lsize]=numpy.fft.ifft2(fftp).ravel().real
    def finalise(self):
        self.estPhase=numpy.reshape(self.x*self.gainfactor,(self.nact,self.nact))
        #av=numpy.average(self.x)
        #print "Average pcg.x: %g - removing it."%av
        #self.x-=av
            
    def solve(self,centroids,usePrevious=1):
        """Get an estimate for the phase...
        This doesn't reset the x vector - assumes that best guess for new
        phase is the previous phase (from last simulation iteration).
        
        """
        self.newCentroids(centroids)
        self.initialise(usePrevious=usePrevious)
        self.compute()
        self.finalise()

        
def computeChat():
    """Compute the C^ matrix (C-hat).

    Shape of Gam will be (nxny, nxny), assumed to operate on the
    phase that has been flattened to a vector.

    """
    dim=nLayers*nx*ny
    Chat=numpy.zeros((dim,dim),numpy.float64)
    GamxHat=xxx
    GamyHat=xxx
    Mhat=xxx#??? FFT of array with flattened subflag on diagonal?

    GamxHatT=numpy.conjugate(numpy.transpose(GamxHat))
    GamyHatT=numpy.conjugate(numpy.transpose(GamyHat))
    MhatT=numpy.conjugate(numpy.transpose(Mhat))

    for Ly in range(nLayers):
        for Lx in range(nLayers):
            for j in range(Ngs):
                P[j,Ly]=xxx
                P[j,Lx]=xxx
                PHat=[xxx,xxx]
                PHatT=numpy.conjugate(numpy.transpose(PHat))
                xbit=quick.dot(quick.dot(quick.dot(GamxHatT,MhatT),Mhat),GamxHat)
                ybit=quick.dot(quick.dot(quick.dot(GamyHatT,MhatT),Mhat),GamyHat)

                Chat[xxx,xxx]+=quick.dot(quick.dot(PHatT,xbit+ybit),PHat)
            Chat[xxx,xxx]*=SHnoise**2
            if Ly==Lx:
                Chat[xxx,xxx]+=1./turbStrength(Lx)**2 * kappa
    return Chat

def reorderChat(Chat):
    """Reorder Chat so that it becomes block diagonal.  This makes
    the solution of Cz=r faster.
    """
    print "TODO"

def computeShiftOperator(dimy,dimx):
    """Compute the shift operator, and its fourier representation.
    Argument dimx,y gives the linear dimension of the phase screen.
    It is assumed that this matrix will operate on a vector, ie the
    phase screen has been unwrapped to a vector size dimx times dimy.
    """
    dim=dimx*dimy
    Sx=numpy.zeros((dim,dim),numpy.int8)
    for i in range(dim):
        Sx[i,(i+dimx)%dim]=1
    Sy=numpy.zeros((dim,dim),numpy.int8)
    for i in range(dim):
        Sy[i,(i+1)%dimx+(i/int(dimx))*dimx]=1
    return Sx,Sy
def computeGradientOperator(dimy,dimx):
    """Compute the gradient operator x and y components.  Here, dimy
    and dimx should be the same.
    I think this is only applicable for Fried geometry.  
    """
    if dimy!=dimx:
        print "pcg: Warning - dimy not equal to dimx, may lead to problems"
    Sx,Sy=computeShiftOperator(dimy,dimx)
    I=numpy.identity(dimy*dimx)
    Gamx=0.5*quick.dot(Sx-I,Sy+I)
    Gamy=0.5*quick.dot(Sy-I,Sx+I)
    return Gamx,Gamy

def computeFFTmatrix(size):
    """Compute the matrix representor for an FFT.  Pre and post
    multiply a 2d array by arr and arrinv to get its FFT."""
    arr=numpy.zeros((size,size),numpy.complex64)
    ss=1/numpy.sqrt(size)
    for y in range(size):
        for x in range(size):
            arr[y,x]=ss*(numpy.cos(2*numpy.pi*y*x/size)-1j*numpy.sin(2*numpy.pi*y*x/size))
    arrinv=numpy.linalg.pinv(arr)#same as generalised_inverse
    return arr,arrinv
def computeMask(subflag):
    """Compute the mask matrix... the input is array of 1 or 0 depending on whether subap is used."""
    subflag=numpy.array(subflag)
    subaps=subflag.flatten()
    dim=subaps.shape[0]
    M=numpy.zeros((dim*2,dim*2),numpy.float64)
    for i in range(dim):
        if subaps[i]:
            M[i,i]=1
            M[i+dim,i+dim]=1
    Mhat=fft.fft2(M)
    return M,Mhat
    
def computeChatFromPokeMatrix(pokemx,nLayers,Ngs,sigma,turbStrength,kappa):
    """This is an AGB attempt to compute the Chat matrix using the
    poke matrix.  This relies on the poke matrix, which I think is an
    okay thing to do, though not quite sure where the propagator comes
    in to it.
    Here, assume pokemx.shape=(nlayers, nGS, ncents, nmodes)
    (note, nmodes==nactuators)
    Why do we have a poke matrix for each layer?  And, why do we have one
    for each guide star?  Is this possible, and what does it mean?  How would
    one create the poke matricees?
    Note, this doesn't work with a zernike reconstructor, needs to be zonal.
    Note, actually, have 1 poke matrix per DM (ie one per layer), each of which
    contains ncents * nGS values... so, some reshaping needed...
    

    nLayers is the number of atmos layers (DMs) to be used in
    tomography, Ngs is the number of guide stars, sigma is the
    standard deviation of the subap noise, turbStrength is the
    turbulance strength for each layer, and kappa is an array
    representing the magnitude of the wave number (appropriately
    scaled) and raised to appropriate power for the power spectrum in
    question (Kol or von Kar).
    
    Once the Chat has been created, can make use of sparsity and get
    rid of some of it (use compactToBanded method).  Chat can probably
    be treated as a banded matrix (you should verify this),
    particularly after an fftshift has been applied.  Note, you may
    need to do some sparse matrix reordering first to get it into a
    banded form.  This will particularly be the case if you have more
    than 1 layer (todo).
    """
    nmodes=pokemx.shape[3]
    fftpokemx=numpy.zeros(pokemx.shape,numpy.complex64)
    Chat=numpy.zeros((nLayers*nmodes,nLayers*nmodes),numpy.complex64)
    for L in range(nLayers):
        for j in range(Ngs):
            fftpokemx[L,j]=numpy.fft.fft2(pokemx[L,j])
    sigmam2=sigma**-2
    for Ly in range(nLayers):
        for Lx in range(nLayers):
            Cpart=Chat[Ly*nmodes:(Ly+1)*nmodes,Lx*nmodes:(Lx+1)*nmodes]
            for j in range(Ngs):
                tmp=sigmam2*quick.dot(numpy.conjugate(numpy.transpose(fftpokemx[Ly,j])),fftpokemx[Lx,j])
                Cpart+=tmp
            if Ly==Lx:#turb covariance stuff...
                Cpart+=1./turbStrength(Lx)**2 * kappa
    #Now that Chat has been created, need to make use of sparsity, and get rid
    #of some of it.  This can probably be treated as a banded matrix, particularly if apply an fftshift. ie compactToBanded(fftshift(Chat),xxx,xxx)
    
    return Chat
def computeChatFromPokeMatrix2(pokemx,sigma,turbStrength,kappa):
    """Similar to computeChatFromPokeMatrix, except the pokemx here
    just has shape ncents, nmodes, removing the need for nlayers and
    nGS dimensions, this being swallowed into the pokemx.

    This effectively relies on the result that we are trying to
    compute Chat from GT C-1 G + Cphi-1.  And notes that GT G is the
    same as F G^* G^ F-1, so what we compute here is effectively G^*
    G^.

    After computing this, one can then make use of sparsity and get
    rid of some of it (use compactToBanded method).  Chat can probably
    be treated as banded (you should verify this), though may need
    some sparse matrix reordering first.

    Note, the pokemx is likely to be block diagonal, because poking
    one mirror does not affect the phase on other mirrors.

    We assume here that the noise on each subap is identical, ie equal
    to sigma, so that the noise covariance matrix is sigma^2 I with I
    being the identity matrix.

    Not really quite sure how this would plan out if using more than 1
    layer - maybe need to FFT the blocks of the pokemx, not the entire
    pokemx (todo).
    """
    turbStr=1./turbStrength**2
    sigmam2=1./(sigma*sigma)
    fftpokemx=numpy.fft.fft2(pokemx)
    fftpokemxCT=numpy.conjugate(numpy.transpose(fftpokemx))
    Chat=sigmam2*quick.dot(fftpokemxCT,fftpokemx)
    modesPerLayer=float(Chat.shape[0])/turbStrength.shape[0]
    if numpy.ceil(modesPerLayer)!=numpy.floor(modesPerLayer):
        print "pcg: warning - modes per layer not an integer."
    modesPerLayer=int(modesPerLayer)
    for i in range(turbStrength.shape[0]):#for each turbulent layer...
        p=i*modesPerLayer
        q=p+modesPerLayer
        Chat[p:q,p:q]+=turbStr[i]*kappa
    return Chat

def solve_banded(A,b,sparsitywidth,iscompact=0):
    """Solves Ax=b for x (A matrix, b vector).  Assumes that A is
    sparse diagonal matrix with a width equal to sparsitywidth, which
    typically will be odd.  Note, A should be stored in banded matrix
    format, ie stored in diagonal ordered form.  This probably means
    that sparsitywidth isn't needed?  Have to test this...
    """
    l=u=sparsitywidth/2#number of elements below and to right of diagonal.
    if iscompact==0:
        Ac=compactToBanded(A,l,u)
    else:
        Ac=A
    x=scipy.linalg.solve_banded((l,u),Ac,b)
    return x

def compactToBanded(a,l,u):
    """Compact a banded matrix into an array
    which can hold the compacted version.  Here, a is the banded
    matrix, l and u refer to the number of elements below and above
    (lower, upper) that are to be used.  The return from this can then be used in scipy.linalg.solve_banded."""
    if a.shape[0]!=a.shape[1]:
        print "pcg: warning - compactToBanded not square matrix."
    banded=numpy.zeros((l+u+1,a.shape[1]),a.dtype.char)
    for i in range(a.shape[1]):
        istart=i-l
        iend=i+u+1
        ostart=0
        oend=banded.shape[0]
        if iend>a.shape[1]:
            oend-=iend-a.shape[1]
            iend-a.shape[1]
        if istart<0:
            ostart=-istart
            istart=0
        banded[ostart:oend,i]=a[i,istart:iend]
    return banded
def computeTurbSpectrum(nmode,actPosition,l0="inf",innerScale=0.1):
    """Compute the grid function, kappa, which represents the
    magnitude of the 2d vector wave number.  This is interpreted as a
    pointwise multiplication operator.

    This will actually depend on DM geometry, scale of DM to the pupil
    etc.

    Basically, compute the distance of each actuator from each other
    one, and form a matrix of this distance.  Then invert it to get
    the wavenumbers...

    actPosition an array of coordinates, such that [i][0] gives the x
    coord of the ith actuator and [i][1] gives the y coord of the ith
    actuator.  The coordinates should be scaled correctly so that this
    works(!!!), ie should be in the correct units.  Sorry thats not
    very helpful, but you can probably work it out from your geometry!

    innerScale is the value of inner scale for turbulence, in same
    units as the actuator coordinates.  If the actuator spacing is
    then less than this value, it is forced to this value, to avoid a
    pole.

    Can specify l0 as a floating point number for a von Karman
    spectrum (in which case, innerScale is ignored).  If using a von
    Karman spectrum, we are computing l0^2/(1+l0^2*r^2), if
    kolmogorov, just computing 1/r^2.  And then raising to power of
    11/6 at the end.
    """
    waveno=numpy.zeros((nmode,nmode),numpy.float32)
    innerScale2=innerScale*innerScale
    if l0=="inf":#kolmogorov spectrum
        for i in range(nmode):
            for j in range(nmode):
                dist2=(actPosition[i][0]-actPosition[j][0])**2+(actPosition[i][1]-actPosition[j][1])**2
                if dist2<innerScale2:
                    dist2=innerScale2
                waveno[i,j]=1./dist2
    else:#von Karman spectrum
        L02=l0*l0
        for i in range(nmode):
            for j in range(nmode):
                dist2=(actPosition[i][0]-actPosition[j][0])**2+(actPosition[i][1]-actPosition[j][1])**2
                if dist2<innerScale2:
                    dist2=innerScale2
                waveno[i,j]=(1+L02/dist2)/L02
    #finally raise to a power...
    waveno=waveno**(11./6)
    #maybe scale by some power of 2pi?
    return waveno
    
def computePropagationMatrix(npup,ncent):
    """Prepares a propagation matrix.  Not exactly sure what this is, though
    suspect that it should act on a phase VECTOR (not matrix), and produce
    another matrix which can then be centroided (or have the LAMBDA matrix
    applied to it).
    """
    mx=numpy.zeros((ncent,npup*npup),numpy.float32)
    nsubx=numpy.sqrt(ncent/2)
    n=npup/nsubx
    #need to combine n elements into 1 for the new matrix, so that the Lambda
    #matrix will work.
    
def setToBanded(mx,width):
    """Set to zero all elements of a matrix except those within width of the
    diagonal"""
    size=min(mx.shape)
    newmx=numpy.zeros((size,size),mx.dtype.char)
    for i in range(size):
        newmx[i,i]=mx[i,i]
        for j in range(width):
            if i+j<size:
                newmx[i,i+j]=mx[i,i+j]
                newmx[i+j,i]=mx[i+j,i]
    return newmx

def test():
    """Create a test PCG object.
    Notes:
    One thing that could be tried is to use the phaseCovariance matrix from
    the bccb approximation."""
    
    pokemx=numpy.array(util.FITS.Read("pokemx.fits")[1]).transpose()
    #phasecov=numpy.array(util.FITS.Read("phaseCovariance.fits")[1])#this one is from a simulation.
    phasecov=numpy.array(util.FITS.Read("phaseCovBCCBTruncated.fits")[1])#this one is created using the BCCB approximation (util.phaseCovariance.bccb), inversed, truncated to 68x68, diagonalised, and reinversed.
    noisecov=numpy.array(util.FITS.Read("noiseCovariance.fits")[1])
    #these matricess were for a 4.2m telescope with 9x9 actuator mirror, with the corner 3 actuators, and central one not used.
    actPosition=numpy.zeros((68,2),numpy.float32)
    print pokemx.shape,phasecov.shape,noisecov.shape
    dmpupil=numpy.array(
        [[0,0,1,1,1,1,1,0,0],
         [0,1,1,1,1,1,1,1,0],
         [1,1,1,1,1,1,1,1,1],
         [1,1,1,1,1,1,1,1,1],
         [1,1,1,1,0,1,1,1,1],
         [1,1,1,1,1,1,1,1,1],
         [1,1,1,1,1,1,1,1,1],
         [0,1,1,1,1,1,1,1,0],
         [0,0,1,1,1,1,1,0,0]])
    pos=0
    ipos=0
    scale=4.2/8
    while pos<68:
        x=ipos%9
        y=ipos/9
        if dmpupil[y,x]==1:
            actPosition[pos,0]=x*scale
            actPosition[pos,1]=y*scale
            pos+=1
        ipos+=1
    Pcg=pcg(pokemx, noisecov, phasecov,None,sparsityWidth=0.25,turbStrength=numpy.ones((1,),numpy.float32))
    return Pcg

def testfull(fromdisk=1,nact=9,pupfn=None,actfn=None,fullpmx=0,noiseCov=None,phaseCov=None,options={},radToPxlFactor=1 ):
    """Create test PCG object
    Notes:
    Uses all actuators (dm_minarea=0.0).
    This is to see whether the scheme needs everything in...
    
    """
    if fromdisk:
        # pokemx=numpy.array(util.FITS.Read("pokemxfull.fits")[1]).transpose()#a poke matrix using all 81 actuators (but not full number of centroids).
        pokemx=numpy.array(util.FITS.Read("pokemxfullfull.fits")[1]).transpose()#a poke matrix using all 81 actuators and 128 centroids.
        # noisecov=numpy.zeros((pokemx.shape[0],),numpy.float64)
        if type(noiseCov)==type(None):
            noiseCov=float(numpy.average(numpy.diagonal(numpy.array(util.FITS.Read("noiseCovariance.fits")[1]))))
        bccb=numpy.array(util.FITS.Read("phaseCovBCCB.fits")[1])#this one is created using the BCCB approximation (util.phaseCovariance.bccb), maintaining info about all 81 actuators
        # Now, inverse the bccb matrix, take the diagonal, and re-inverse it, giving the approximate phase covariance matrix...
        phasecov=1./numpy.diagonal(numpy.linalg.inv(bccb))
    else:
        pokemx=util.createPokeMx.straightPokes(nact,pupfn=pupfn,actmap=actfn,returnfull=fullpmx).transpose()
        if type(noiseCov)==type(None):
            noiseCov=0.25
        if type(phaseCov)==type(None) or phaseCov=="bccb" or phaseCov=="bccbAll":
            bccb=util.phaseCovariance.bccb(nact=nact)
            #bccb.savespace(1)
            # Now, inverse the bccb matrix, take the diagonal, and re-inverse it, giving the approximate phase covariance matrix...
            if phaseCov=="bccbAll":
                phasecov=bccb
            else:
                phasecov=1./numpy.diagonal(numpy.linalg.inv(bccb))
        elif phaseCov=="bttb":
            bttb=util.phaseCovariance.bttb(nact=nact)
            phasecov=bttb
            #phasecov.savespace(1)
        else:
            phasecov=phaseCov
    phasecov*=radToPxlFactor#to get into pixels...
    noiseCov*=radToPxlFactor#to get into pixels (from radians)
    Pcg=pcg(pokemx, noiseCov, phasecov,None,sparsityWidth=0.25,turbStrength=numpy.ones((1,),numpy.float32),options=options)
    return Pcg

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
        pupil.fn[:,]=1
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
            zernlist=science.infScrn.makeInitialScreen(dpix=(nact-1)*8,Dtel=4.2,L0=30.,scrnXPxls=None,scrnYPxls=None,seed=None,tstep=0.05,globR0=0.2,strLayer=1.,natype=numpy.float64,windDirection=0.,vWind=10.)[:-1,1:,].copy()
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
#simcents2[64:,]=Numeric.reshape(simcents[:,:,1],(64,))
#Pcg=pcg.doit(9,simcents2,Numeric.reshape(actuators,(9,9))) 
    
print "TODO: try with a pokemx that has unused modes ignored."
