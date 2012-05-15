import base.aobase as aobase
from util.flip import fliparray
import cmod.mkimg,cmod.imgnoise,cmod.utils,cmod.binimg
import numpy as na
#import Numeric as Nmeric
import numpy.fft as fft
import numpy.random as nara
import scipy.linalg as LA
import numpy.linalg as NALA

class mapRecon(aobase.aobase):
    """
    A class which carries out a MAP (minimum a-posteri) reconstruction.
    The reconstruction matrix (which is applied to the centroids, to give the
    new DM shape is given by Francois Assemat PhD thesis as:
    W = G^-1 <b_G p^T> (<pp^T> + C_b)^-1
    where G is the geometry covariance matrix for the basis forming the
    influence functions of the DM, <pp^T> is the covariance of slopes and
    C_b is the covariance of noise.
    The part in parenthesis, <pp^T> + C-b is given as
    <mm^T> in the thesis.

    Further note:
    The (<pp>^T + C_b) matrix can be computed by montecarlo when sim is running, simply by getting slope covariance from eg 1000-10000 cycles with a flat mirror.
    The <b_G p^T> matrix can be computed from the poke matrix when sim is running, ie by poking the mirror with zernikes, and getting the (noiseless) slopes.  This is because <bp^T> = <ff^TD^T> = <ff^T>D^T where f is the mirror mode and D is the poke matrix since p = Df and f = sum (a_i Z_i).

    Note, the standard definition for a MAP reconstructor is:
    ( H^t C_n^-1 H + C_c^-1 )^-1 H^T C_n^-1
    Here, H is a poke matrix or similar.
    C_n is the wfs noise covariance matrix.
    C_c is the phase covariance matrix.
    However, this hasn't been implemented since FA's method was supposedly superior (although I should think significantly harder to get a working implementation).
    
    """
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None):
        if type(parent)!=type({}):
            parent={"cent":parent}
        aobase.aobase.__init__(self,parent,config)
        self.dataValid=1#data is valid even before first iter as assume mirror starts zerod.
        self.nMirrorModes=self.config.getVal("nMirrorModes")
        if forGUISetup==1:
            self.outputData=[(self.nMirrorModes,),na.float64]
        else:
            self.nsubx=self.config.getVal("wfs_nsubx")
            n=self.wfs_n=self.config.getVal("wfs_n")
            self.outputData=na.zeros((self.nMirrorModes,),na.float64)
            self.dmType=self.config.getVal("dmType",default="zernike")
            self.pupil=self.config.getVal("pupil")
            self.subarea=self.pupil.subarea#na.zeros((self.nsubx,self.nsubx),na.float64)
            self.subflag=self.pupil.subflag#na.zeros((self.nsubx,self.nsubx),na.int8)
            self.pupsub=self.pupil.pupsub#na.zeros((self.nsubx,self.nsubx,self.wfs_n,self.wfs_n),na.float64)
            self.dmpupil=self.pupil.dmpupil#na.zeros(self.pupil.fn.shape,na.int32)
            if self.subflag.itemsize==8:
                print "WARNING: untested with 8 byte longs...(wfs)"
            wfs_minarea=self.config.getVal("wfs_minarea")#0.5...                                      # Min unvignetted subap area to use - why here ????
##             for i in range(self.nsubx):        
##                 for j in range(self.nsubx):
##                     self.pupsub[i][j]=self.pupil.fn[i*n:(i+1)*n,j*n:(j+1)*n].astype(na.float64)    # Get pupil fn over subaps
##                     self.subarea[i][j]=na.sum(na.sum(self.pupsub[i][j]))
##                     if(self.subarea[i][j]>(wfs_minarea*n*n)):
##                         self.subflag[i][j]=1# Flag vignetted subaps
##                         self.dmpupil[i*n:(i+1)*n,j*n:(j+1)*n]=self.pupil.fn[i*n:(i+1)*n,j*n:(j+1)*n].astype(na.int32)
            #get the number of atmospheric layers to be used in reconstruction.  For a single WFS case, this is just one (no 3d profiling).
            self.nReconLayers=self.config.getVal("nReconLayers")
            self.wfsn=self.config.getVal("wfs_n")
            self.nsubaps=self.nsubx*self.nsubx
            self.r0=self.config.getVal("r0")
            self.npup=self.config.getVal("npup")
            self.telDiam=self.config.getVal("telDiam")
            if self.dmType=="zernike":
                import zernikeDM
                self.dm=zernikeDM.dm(self.pupil,self.nMirrorModes)#create a DM object
            else:
                print "DM type %s not known - assuming zernike"%self.dmType
                self.dm=zernikeDM.dm(self.pupil,self.nMirrorModes)
            #get the mirror mode area - for a continuous mirror, this will be the whole mirror, while for a segmented mirror, it will be the area of a segment.
            self.mirrorModeArea=self.dm.computeModeArea()
            self.turbProfileType=self.config.getVal("turbProfileType",default="kolmogorov")
            #get the turb profile information.  This is a dictionary, with keys being the order of differentiation - only 1 and 2 are used here.
            #The values are then dictionaries of what is the differentiating variable, ie x or y.
            #The values of these are then arrays of the function.
            self.turbProfile=[]
            if self.turbProfileType=="kolmogorov":#actually, create one of these per layer (or 1 if using single wfs).
                for i in range(self.nReconLayers):#actually, do something here to apply a different r0 for each layer...
                    self.turbProfile.append(self.kolmogorovProfile(self.r0/self.telDiam*self.npup,self.npup*4))
            else:
                self.turbProfile=self.turbProfileType#assume type stores the actual profile.
            self.wfsSearchPath=self.config.getVal("wfsSearchPath")#this will be something like [["wfscent","globals"],...] depending on the wfs module used.  It is used to create a dummy wfs to compute noise statistics, and so should allow the values for a given wfs to be obtained.  When using more than one wfs, they should be appended in the list, eg [["wfscent_1","wfscent","globals"],["wfscent_2",...]]



            # compute the G^-1 part.
            print "computing dm inv geom cov matrix"#shape=(nmode-1,nmode-1) - piston removed
            self.invGeomCovarianceMatrix=self.dm.computeInvGeomCovarianceMatrix()
            #compute the noise covariance matrix - shape=(ncent,)
            print "computing noise covariances"
            self.noiseCovarianceMatrix=self.computeNoiseCovariance()
            #compute the slope covariance matrix - shape=(ncent,ncent)
            print "computing slope covariance matrix"
            self.slopeCovarianceMatrix=self.computeSlopeCovarianceMatrix()
            #add noise to slope
            print "adding noise and slope covariances"
            for i in range(self.nsubaps*2):
                self.slopeCovarianceMatrix[i,i]+=self.noiseCovarianceMatrix[i]
                if self.noiseCovarianceMatrix[i]==na.NaN:
                    print "mapRecon - noiseCovarianceMatrix has NaN at %d"%i
            #compute the invrse of the noise and slope covariance matrix
            #print na.array(self.slopeCovarianceMatrix)
            #print na.sum(na.sum(self.slopeCovarianceMatrix-na.transpose(self.slopeCovarianceMatrix)))
            print "inverting covariance matrix"#shape=(ncent,ncent)
            try:
                self.invNoiseSlopeMatrix=LA.cho_solve(LA.cho_factor(self.slopeCovarianceMatrix),na.identity(self.slopeCovarianceMatrix.shape[0]))
            except:
                print "noiseSlopeCovarianceMatrix may not be positive definite - trying generatlised inverse"
                self.invNoiseSlopeMatrix=NALA.generalized_inverse(self.slopeCovarianceMatrix)
            #compute the slope/dm covariance matrix - shape=(nmode-1,ncent)
            print "computing slope/dm cov mat"
            self.slopeDMCovMatrix=self.computeSlopeDMCovarianceFunction()
            #multiply together to get reconstructor
            print "multiplying to create reconstructor"
            #print self.slopeDMCovMatrix.shape,self.invNoiseSlopeMatrix.shape
            #print na.sum(na.sum(self.slopeDMCovMatrix)),na.sum(na.sum(self.invNoiseSlopeMatrix))
            tmpProd=na.dot(self.slopeDMCovMatrix,self.invNoiseSlopeMatrix)#shape=(nmode-1,ncent)
            #print self.invGeomCovarianceMatrix.shape,tmpProd.shape
            #print na.sum(na.sum(self.invGeomCovarianceMatrix)),na.sum(na.sum(tmpProd))
            self.reconstructor=na.dot(self.invGeomCovarianceMatrix,tmpProd)#shape=(nmode-1,ncent)
            self.inputData=na.zeros((self.nsubaps*2,),na.float64)
            self.centx=self.inputData[:self.nsubaps]
            self.centy=self.inputData[self.nsubaps:]

        
    def computeNoiseCovariance(self,detectorNo=0):
        """Compute the noise covariance matrix - this is a diagonal matrix since noise is uncorrelated.  This is the C_b part.  The result is in pixels I think."""
        arr=na.zeros((2*self.nsubaps,),na.float64)
        searchOrder=self.config.searchOrder
        self.config.setSearchOrder(self.wfsSearchPath[detectorNo])
        wfs=simpleWFS(self.config)
        self.config.setSearchOrder(searchOrder)
        wfsn=self.wfsn
        for i in range(self.nsubaps):
            subx=i%self.nsubx
            suby=i/self.nsubx
            pupsub=self.pupil.fn[suby*wfsn:suby*wfsn+wfsn,subx*wfsn:subx*wfsn+wfsn].astype(na.float64)
            x,y=wfs.createStats(20,pupsub)
            arr[i]=x
            arr[i+self.nsubaps]=y
        return arr
    
    def computeSlopeCovarianceMatrix(self):
        """Compute the slope covariance matrix.  This is the <pp^T> part.  The result is in pixels."""
        n=self.nsubaps
        wfsn=self.wfsn
        nsubx=self.nsubx
        matAll=na.zeros((2*n,2*n),na.float64)
        matrixpiX=matAll[:n,:n]
        matrixpiY=matAll[n:,n:]
        matrixpiXY=matAll[:n,n:]
        #create arrays for input to convolution of the 2 subaps.  Note, these are oversized for aliasing etc.
        mask0=na.zeros((self.nsubx*self.wfsn*2*4,self.nsubx*self.wfsn*2*4),na.float32)
        mask1=na.zeros((self.nsubx*self.wfsn*2*4,self.nsubx*self.wfsn*2*4),na.float32)
        #print self.pupil.fn
        D=self.wfsn*self.npup/self.telDiam#subap diam
        for k0 in range(n):#iterate over all subaps, comparing with all others
            subx0=k0%nsubx
            suby0=k0/nsubx
            # mask the pupil function for the current subaperture.  We change only the central portion of mask, making it equal to the subap shape.
            s=mask0.shape
            #put the current subap shapes into the centre of an array
            mask0[s[0]/2-wfsn/2:s[0]/2+wfsn/2,s[1]/2-wfsn/2:s[1]/2+wfsn/2]=self.dmpupil[suby0*wfsn:suby0*wfsn+wfsn,subx0*wfsn:subx0*wfsn+wfsn].astype(na.float32)
            #compute the total number of pxls.
            S0=na.sum(na.sum(self.dmpupil[suby0*wfsn:suby0*wfsn+wfsn,subx0*wfsn:subx0*wfsn+wfsn]))
            for k1 in range(n):
                #get coords for subap x and y position
                subx1=k1%nsubx
                suby1=k1/nsubx
                # mask the pupil function for the current subaperture.  We change only the central portion of mask, making it equal to the subap shape.
                s=mask1.shape
                mask1[s[0]/2-wfsn/2:s[0]/2+wfsn/2,s[1]/2-wfsn/2:s[1]/2+wfsn/2]=self.dmpupil[suby1*wfsn:suby1*wfsn+wfsn,subx1*wfsn:subx1*wfsn+wfsn].astype(na.float32)
                #compute the total number of pxls.
                S1=na.sum(na.sum(self.dmpupil[suby1*wfsn:suby1*wfsn+wfsn,subx1*wfsn:subx1*wfsn+wfsn]))
                #get the intercorrelation of the subaperture shapes (in pxls).
                I=fft.fftshift(fft.ifft2(fft.fft2(mask0/S0)*na.conjugate(fft.fft2(mask1/S1)))).real
                
                
                for l in range(self.nReconLayers):#note - for a sngle dm system, use the global r0, and only 1 layer...
                    #compute the coords we're interested in...
                    posy=I.shape[0]/2+(suby1-suby0)*self.wfsn
                    posx=I.shape[1]/2+(subx1-subx0)*self.wfsn
                    #convolve the 2nd derivative of turbulence profile with subap function.  Then, take only the interesting part of this.
                    tmp=fft.fftshift(self.convolve(self.turbProfile[l]["XX"],I))
                    tmp=tmp[posy,posx].real#imag part should be zero...
                    matrixpiX[k0,k1]+=tmp
                    tmp=fft.fftshift(self.convolve(self.turbProfile[l]["YY"],I))
                    tmp=tmp[posy,posx].real
                    matrixpiY[k0,k1]+=tmp
                    tmp=fft.fftshift(self.convolve(self.turbProfile[l]["XY"],I))
                    tmp=tmp[posy,posx].real
                    matrixpiXY[k0,k1]+=tmp
                if S0>0 and S1>0:
                    matrixpiX[k0,k1]*=0.5*D*D/S0/S1
                    matrixpiY[k0,k1]*=0.5*D*D/S0/S1
                    matrixpiXY[k0,k1]*=0.5*D*D/S0/S1
                else:
                    matrixpiX[k0,k1]=0.
                    matrixpiY[k0,k1]=0.
                    matrixpiXY[k0,k1]=0.
                    
        #copy into last quadrant...
        matAll[n:,:n]=matrixpiXY
        
        return matAll
        

    def agbComputeSlopeDMCovarianceFunction(self):
        """Compute the slope/mirror covariance function.  This is the <b_G p^T> part
        This tries to compute cross correlation of turb profile with subap mask, and then cross correlate the mode influence function with the result.  This is what AGB worked it out to be, but doesn't seem to work.
        """
        n=self.nsubaps
        nsubx=self.nsubx
        wfsn=self.wfsn
        subapmask=na.zeros((self.nsubx*wfsn*2*4,self.nsubx*wfsn*2*4),na.float32) 
        modemask=na.zeros((self.nsubx*wfsn*2*4,self.nsubx*wfsn*2*4),na.float64)
        D=self.wfsn*self.npup/self.telDiam#subap diam / m
        matAll=na.zeros((self.nMirrorModes,2*n),na.float64)
        matrixX=matAll[:,:n]
        matrixY=matAll[:,n:]
        #iterate over subaps
        for k in range(self.nsubaps):
            subx0=k%nsubx
            suby0=k/nsubx
            #create the mask
            s=subapmask.shape
            S=na.sum(na.sum(self.dmpupil[suby0*wfsn:suby0*wfsn+wfsn,subx0*wfsn:subx0*wfsn+wfsn]))
            subapmask[s[0]/2-wfsn/2:s[0]/2+wfsn/2,s[1]/2-wfsn/2:s[1]/2+wfsn/2]=self.dmpupil[suby0*wfsn:suby0*wfsn+wfsn,subx0*wfsn:subx0*wfsn+wfsn].astype(na.float32)
            #convolve subap with phase...
            #Qx=fft.fftshift(self.convolve(self.turbProfile[0]["X"],subapmask/S))
            #FFTQx=fft.fft2d(Qx)
            #Qy=fft.fftshift(self.convolve(self.turbProfile[0]["Y"],subapmask/S))
            #FFTQy=fft.fft2d(Qy)
            #cross correlate phase with subap
            #Qx=fft.fftshift(fft.fft2d(na.conjugate(fft.inverse_fft2d(self.turbProfile[0]["X"]))*fft.inverse_fft2d(subapmask/S)))
            #Qy=fft.fftshift(fft.fft2d(na.conjugate(fft.inverse_fft2d(self.turbProfile[0]["Y"]))*fft.inverse_fft2d(subapmask/S)))
            #Qx=fft.fftshift(self.crossCorrelate(self.turbProfile[0]["X"],subapmask/S))
            #Qy=fft.fftshift(self.crossCorrelate(self.turbProfile[0]["Y"],subapmask/S))
            Qx=self.crossCorrelate(self.turbProfile[0]["X"],subapmask/S)
            Qy=self.crossCorrelate(self.turbProfile[0]["Y"],subapmask/S)
            iFFTQx=fft.ifft(Qx)
            iFFTQy=fft.ifft(Qy)
            # store some values for debugging
            self.subapmask=subapmask
            self.Qx=Qx
            self.iFFTQx=iFFTQx
            for m in range(1,self.nMirrorModes):
                s=modemask.shape
                # create mirror mode...
                modemask[s[0]/2-self.nsubx*wfsn/2:s[0]/2+self.nsubx*wfsn/2,s[1]/2-self.nsubx*wfsn/2:s[1]/2+self.nsubx*wfsn/2]=self.dm.getMode(m)/self.dm.normFact[m]
                posy=modemask.shape[0]/2+(suby0-self.nsubx/2)*self.wfsn
                posx=modemask.shape[1]/2+(subx0-self.nsubx/2)*self.wfsn
                #compute the cross correlation of the mirror mode conjugated, with the previously created Q function.

                #tmp=fft.fftshift(fft.fft2d(na.conjugate(fft.inverse_fft2d(modemask))*iFFTQx))
                tmp=fft.fft2(na.conjugate(fft.ifft2(modemask))*iFFTQx)
                # store values for debugging
                self.tmp=tmp
                self.posy=posy
                self.posx=posx
                self.modemask=modemask
                tmp=tmp[posy,posx].real
                matrixX[m,k]+=tmp
                #tmp=fft.fftshift(fft.fft2d(na.conjugate(fft.inverse_fft2d(modemask)))*iFFTQy)
                tmp=fft.fft2(na.conjugate(fft.ifft2(modemask))*iFFTQy)
                tmp=tmp[posy,posx].real
                matrixY[m,k]+=tmp
                if S==0:
                    matrixX[m,k]=0.
                    matrixY[m,k]=0.
                else:
                    matrixX[m,k]*=-0.5*D/self.mirrorModeArea/S
                    matrixY[m,k]*=-0.5*D/self.mirrorModeArea/S
        return matAll
                
    def crossCorrelate(self,a,b):
        """Computes the cross correlation of f* x g where f* is the conjugate of a function f and x represents the cross correlation"""
        if a.shape!=b.shape:
            print "WARNING in mapRecon.crossCorrelate - different shapes for convolution: %s %s"%(str(a.shape),str(b.shape))
        return fft.fft2(na.conjugate(fft.ifft2(a))*fft.ifft2(b))

    def computeSlopeDMCovarianceFunction(self,mode="working"):
        """Compute the slope/mirror covariance function.  This is the <b_G p^T> part, taken from FA thesis.  Seems not quite correct...
        The parameter mode can be "correct" (which has equations as those from FA thesis), "working" (which gives best agreement with monte-carlo).
        """
        n=self.nsubaps
        nsubx=self.nsubx
        wfsn=self.wfsn
        subapmask=na.zeros((self.nsubx*wfsn*2*4,self.nsubx*wfsn*2*4),na.float32) 
        modemask=na.zeros((self.nsubx*wfsn*2*4,self.nsubx*wfsn*2*4),na.float64)
        D=self.wfsn*self.npup/self.telDiam#subap diam / m
        matAll=na.zeros((self.nMirrorModes,2*n),na.float64)
        matrixY=matAll[:,:n]
        matrixX=matAll[:,n:]
        #iterate over modes, comparing each with every other.
        for m in range(1, self.nMirrorModes):#node mode 0 is piston...
            s=modemask.shape
            #create mirror mode...
            modemask[s[0]/2-self.nsubx*wfsn/2:s[0]/2+self.nsubx*wfsn/2,s[1]/2-self.nsubx*wfsn/2:s[1]/2+self.nsubx*wfsn/2]=self.dm.getMode(m)/self.dm.normFact[m]
            if m in [3]:#defocus - not sure why this is needed.  Bit of a bodge really!  
                modemask*=-1.
            self.modemask=modemask#save for debugging
            # does the mode mask need dividing by anything?
            if mode=="working":
                conjFFTmirrormode=na.conjugate(fft.fft2(modemask))#removed shift
            elif mode=="correct":
                conjFFTmirrormode=na.conjugate(fft.ifft2(modemask))#removed shift
            self.conjFFTmirrormode=conjFFTmirrormode
            for k in range(n):
                subx0=k%nsubx
                suby0=k/nsubx
                #create the mask
                s=subapmask.shape
                subapmask[s[0]/2-wfsn/2:s[0]/2+wfsn/2,s[1]/2-wfsn/2:s[1]/2+wfsn/2]=self.dmpupil[suby0*wfsn:suby0*wfsn+wfsn,subx0*wfsn:subx0*wfsn+wfsn].astype(na.float32)
                
                S=na.sum(na.sum(self.dmpupil[suby0*wfsn:suby0*wfsn+wfsn,subx0*wfsn:subx0*wfsn+wfsn]))
                self.subapmask=subapmask/S
                if mode=="working":
                    FFTsubapmask=fft.fft2(subapmask/float(S))
                elif mode=="correct":
                    FFTsubapmask=fft.ifft2(subapmask/float(S))
                self.FFTsubapmask=FFTsubapmask
                prod=na.conjugate(conjFFTmirrormode)*na.conjugate(FFTsubapmask)#this one is wrong, but seems to give better results!
                #prod=conjFFTmirrormode*FFTsubapmask
                self.prod=prod
                if mode=="working":
                    I=fft.fftshift(fft.ifft2(prod))
                elif mode=="correct":
                    I=fft.fftshift(fft.fft2(prod))
                self.I=I
                #print "I:",na.sum(na.sum(I))
                #iterate over all layers... (or just one if only using single wfs)
                for l in range(self.nReconLayers):
                    posy=I.shape[0]/2+(suby0-self.nsubx/2)*self.wfsn
                    posx=I.shape[1]/2+(subx0-self.nsubx/2)*self.wfsn
                    tmp=fft.fftshift(self.convolve(self.turbProfile[l]["X"],I))
                    self.tmp=tmp
                    self.posy=posy
                    self.posx=posx
                    tmp=tmp[posy,posx].real
                    matrixX[m,k]+=tmp
                    tmp=fft.fftshift(self.convolve(self.turbProfile[l]["Y"],I))
                    self.tmp2=tmp
                    tmp=tmp[posy,posx].real
                    matrixY[m,k]+=tmp
                if S==0:
                    matrixX[m,k]=0.
                    matrixY[m,k]=0.
                else:
                    matrixX[m,k]*=0.5*D/self.mirrorModeArea/S
                    matrixY[m,k]*=0.5*D/self.mirrorModeArea/S
        return matAll[1:]#remove piston


    def kolmogorovProfile(self,r0,n):
        """Compute the functions for derivatives of the structure function for
        a kolmog atmosphere.  Input r0 is the r0 for this layer (in pixels) (or global), and n is the number of pixels in the telescope aperture (oversampling is advisable)."""
        #Create individual functions - these are from Assemat thesis.
        dims=(n*2,n*2)
        
        def dDdx(x,y,r0=r0,sx=n,sy=n):
            #if x==0 and y==0:
            #    return 0.
            tmp=6.88*r0**(-5/3.)*5*(x-sx)/(3*((x-sx)*(x-sx)+(y-sy)*(y-sy))**(1/6.))
            tmp[sx,sy]=0.
            return tmp
        def dDdy(x,y,r0=r0,sx=n,sy=n):
            #if x==0 and y==0:
            #    return 0.
            tmp=6.88*r0**(-5/3.)*5*(y-sy)/(3*((x-sx)*(x-sx)+(y-sy)*(y-sy))**(1/6.))
            tmp[sx,sy]=0.
            return tmp
        def ddDdxdx(x,y,r0=r0,sx=n,sy=n):
            #if x==0 and y==0:
            #    x=1e-6
            tmp=6.88*r0**(-5/3.)*5*(2*(x-sx)*(x-sx)+3*(y-sy)*(y-sy))/(9*((x-sx)*(x-sx)+(y-sy)*(y-sy))**(7/6.))
            tmp[sx,sy]=6.88*r0**(-5/3.)*5*(2*1e-12)/(9*(1e-12)**(7/6.))
            tmp[sx,sy]=(tmp[sx-1,sy]+tmp[sx+1,sy]+tmp[sx,sy-1]+tmp[sx,sy+1])/4
            return tmp
        def ddDdxdy(x,y,r0=r0,sx=n,sy=n):
            #if x==0 and y==0:
            #    return 0.
            tmp=6.88*r0**(-5/3.)*-5*(x-sx)*(y-sy)/(9*((x-sx)*(x-sx)+(y-sy)*(y-sy))**(7/6.))
            tmp[sx,sy]=0.
            return tmp
        def ddDdydy(x,y,r0=r0,sx=n,sy=n):
            #if x==0 and y==0:
            #    x=1e-6
            tmp=6.88*r0**(-5/3.)*5*(2*(x-sx)*(x-sx)+3*(y-sy)*(y-sy))/(9*((x-sx)*(x-sx)+(y-sy)*(y-sy))**(7/6.))
            tmp[sx,sy]=6.88*r0**(-5/3.)*5*(2*1e-12)/(9*(1e-12)**(7/6.))
            tmp[sx,sy]=(tmp[sx-1,sy]+tmp[sx+1,sy]+tmp[sx,sy-1]+tmp[sx,sy+1])/4
            return tmp
        #now create the structure - dictionary of dictionary of arrays.  First key is the order of differentiation, 2nd is the variable used for differentiation.
        d={"X":na.fromfunction(dDdx,dims),"Y":na.fromfunction(dDdy,dims),"XX":na.fromfunction(ddDdxdx,dims),"YY":na.fromfunction(ddDdydy,dims),"XY":na.fromfunction(ddDdxdy,dims)}
        #print d
        return d
    def convolve(self,a,b):
        if a.shape!=b.shape:
            print "WARNING in mapRecon.convolve - different shapes for convolution: %s %s"%(str(a.shape),str(b.shape))
        return fft.ifft2(fft.fft2(a)*fft.fft2(b))



    def newParent(self,parent,idstr=None):
        if type(parent)!=type({}):
            parent={"cent":parent}
        self.parent=parent

    def generateNext(self,msg=None):
        """Data coming in from parents can be of 2 types:
        - centroid data (lgs)
        - tilt sensitivity data (ngs)
        Determine which is which from the name of the parent object
        (the dictionary key).  If cent*, is centroid data, if
        tiltsens*, is tilt sensitivity data.
        """
        if self.generate==1:
            if self.newDataWaiting:
                nin=0
                self.useMonteRecon=0
                for key in self.parent.keys():
                    if key=="reconMatrix":
                        if self.parent[key].dataValid==1:
                            self.useMonteRecon=1
                    else:
                        if self.parent[key].dataValid==1:
                            nin+=1
                        else:
                            print "zernike Recon: Waiting for data from wfs, but not valid"
                if nin>0:
                    if nin==len(self.parent.keys())-self.parent.has_key("reconMatrix"):
                        self.dataValid=1
                    else:
                        print "zernike Recon: Warning - got some data but not all, setting dataValid=0"
                        self.dataValid=0
                else:
                    self.dataValid=0
            if self.dataValid:
                self.calc()# Calculate DM updates
            else:
                self.dataValid=0

    def calc(self):
        """Carry out the reconstruction at each simulation iteration"""
        #self.centx[:]=0.
        #self.centy[:]=0.
        self.inputData[:]=0.
        for key in self.parent.keys():  #average the centroids together...
            if key!="reconMatrix":
                self.inputData+=na.array(self.parent[key].outputData.flat)
            #self.centx += self.parent[key].outputData[0]#tempcentx
            #self.centy += self.parent[key].outputData[1]#tempcenty
        self.inputData/=float(len(self.parent.keys())-self.parent.has_key("reconMatrix"))
        #self.centx /= float(len(self.parent.keys()))  #centproc_list))
        #self.centy /= float(len(self.parent.keys()))  #centproc_list))
        self.doReconstruction()
    def doReconstruction(self):
        """Actually do the matrix multiplication..."""
        #xxx don't know what sign to use here.
        if self.useMonteRecon:
            self.outputData[1:]=-na.array(na.dot(self.parent["reconMatrix"].outputData,self.inputData))
        else:
            self.outputData[1:]=-na.array(na.dot(self.reconstructor,self.inputData))
        
class simpleWFS:
    def __init__(self,config):
        """A class to implement a simple centroider - purely for getting noise stats...
        Note, config here should have a getVal method.  If using the standard simulation config, make sure you set the correct searchOrder before creating this object."""
        self.config=config
        self.wfs_n=self.config.getVal("wfs_n")
        self.wfs_nfft=self.config.getVal("wfs_nfft")
        self.wfs_nimg=self.config.getVal("wfs_nimg")
        self.wfs_ncen=self.config.getVal("wfs_ncen")
        self.xtiltfn=((na.fromfunction(self.tilt,(self.wfs_n,self.wfs_n))-float(self.wfs_n)/2.+0.5)/float(self.wfs_n)).astype(na.float64)# subap tilt fn
        self.ytiltfn=na.transpose(self.xtiltfn)
        self.subimg=na.zeros((self.wfs_nfft,self.wfs_nfft),na.float64)      # SH sub-images (high LL)
        self.fftwPlan=cmod.mkimg.setup(self.subimg) # Setup imaging FFTs for subaps
        self.pupsub=na.zeros((self.wfs_n,self.wfs_n),na.float64)
        self.bimg=na.zeros((self.wfs_nimg,self.wfs_nimg),na.float64)
        self.skybrightness=self.config.getVal("wfs_skybrightness")
        self.wfs_read_mean=self.config.getVal("wfs_read_mean")
        self.wfs_read_sigma=self.config.getVal("wfs_read_sigma")
        self.wfs_bandwidth=self.config.getVal("wfs_bandwidth")
        self.wfs_thruput=self.config.getVal("wfs_thruput")
        self.wfs_phot_rate_factor=self.config.getVal("wfs_phot_rate_factor")
        self.telDiam=self.config.getVal("telDiam")
        self.wfs_nsubx=self.config.getVal("wfs_nsubx")
        self.wfs_int=self.config.getVal("wfs_int")
        self.wfs_mag=self.config.getVal("wfs_mag")
        self.wfs_floor=self.config.getVal("wfs_floor")
        self.cenmask=na.zeros((self.wfs_nimg,self.wfs_nimg),na.float32)             # Define centroiding mask
        self.cenmask[self.wfs_nimg/2-self.wfs_ncen/2:self.wfs_nimg/2+self.wfs_ncen/2,self.wfs_nimg/2-self.wfs_ncen/2:self.wfs_nimg/2+self.wfs_ncen/2]=1.
        self.sig=self.wfs_sig()
        
    def tilt(self,x,y):
        return y

    def wfs_sig(self):
        """wfs signal per exposure"""
        bandwidth=self.wfs_bandwidth                                  # Optical bandwidth (Angstrom)
        thruput=self.wfs_thruput                                      # Optical thruput incl. DQE
        rate_factor=self.wfs_phot_rate_factor                         # Photons/cm^2/s/Angstrom
        pupil_area=(self.telDiam/float(self.wfs_nsubx)*100.)**2.      # Pupil area/cm^2
        sig=rate_factor*thruput*bandwidth*self.wfs_int*pupil_area/(2.5**self.wfs_mag) # Detected image photons
        #print 'WFS  photons/subap/integration: ',sig
        return sig
        
    def createStats(self,nimg,pupsub):
        """nimg is the number of times to do it to get the stats, and pupsub is float64 array size wfs_n x wfs_n defining the pupil"""
        bimg=self.bimg
        tmp=0.5*float(self.wfs_n)/float(self.wfs_nfft)*2.*na.pi
        phssub=(-tmp*self.xtiltfn-tmp*self.ytiltfn).astype(na.float64)
        sumx=sumx2=0.
        sumy=sumy2=0.
        for i in range(nimg):
            cmod.mkimg.mkimg(self.fftwPlan,phssub,self.subimg,pupsub.astype("d"))
            img=fliparray(self.subimg)
            cmod.binimg.binimg(img,bimg)
            totsig=na.sum(na.sum(bimg))
            if totsig>0:
                bimg*=self.sig/totsig
            bimg+=self.skybrightness
            if totsig>0 or self.skybrightness>0:
                cmod.imgnoise.shot(bimg,bimg)
            bimg+=nara.normal(self.wfs_read_mean,self.wfs_read_sigma,bimg.shape)
            cimg=na.where(bimg<self.wfs_floor,0,self.bimg-self.wfs_floor)*self.cenmask
            xcent,ycent=cmod.utils.centroid(cimg.astype("f"))
            #print xcent,ycent
            sumx+=xcent
            sumy+=ycent
            sumx2+=xcent*xcent
            sumy2+=ycent*ycent
        sumx/=nimg
        sumy/=nimg
        varx=sumx2/nimg-sumx*sumx
        stdx=na.sqrt(varx)
        vary=sumy2/nimg-sumy*sumy
        stdy=na.sqrt(vary)
        return varx,vary
        
