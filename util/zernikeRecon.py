"""Implement a theoretical zernike reconstructor (no poke matrix needed).
"""
import os,numpy,numpy.linalg
import util.tel
from cmod.zernike import zern
import cmod.binimg
import util.FITS
#import util.dot as quick

class zernikeRecon:
    def __init__(self,npup,nsubx,nzrad=None,nz=None,recon_jmin=1,pupil=None,zern_tilts=1,computeReconmx=1):
        """recon_jmin - minimum zernike to use for reconstruction.  1 if want
        to include tip/tilt and 3 otherwise.
        
        """
        self.pupil=pupil
        if type(self.pupil)==type(None):
            self.pupil=util.tel.Pupil(npup,npup/2,0,nsubx)
        self.npup=npup
        self.nsubx=nsubx
        n=npup/nsubx
        self.ndata=self.pupil.ndata
        self.subarea=self.pupil.subarea
        self.subflag=self.pupil.subflag
        self.data=numpy.zeros(self.ndata,numpy.float64)
        self.dmphs=numpy.zeros((npup,npup),numpy.float64)
        self.xtiltfn=(numpy.fromfunction(lambda x,y:y,(n,n))-float(n)/2.+0.5)/float(n)
        self.ytiltfn = numpy.transpose(self.xtiltfn)          # Subap tilt fns in x & y
        self.outputData=numpy.zeros((3,nsubx,nsubx),numpy.float64)
        self.pist=self.outputData[0]  # numpy.zeros((nsubx,nsubx),numpy.float64)
        self.xtilt=self.outputData[1] # numpy.zeros((nsubx,nsubx),numpy.float64)
        self.ytilt=self.outputData[2] # numpy.zeros((nsubx,nsubx),numpy.float64)

        self.recon_jmin=recon_jmin
        if nzrad==None:
            nzrad=int(nsubx*0.75)
        self.nzrad=nzrad
        nzmax=self.nzmax=(self.nzrad+1)*(self.nzrad+2)/2#self.config.getVal("nzmax")#max number of zernikes during initialisation.
        if nz==None:
            nz=self.nzmax
        self.nz=nz
        self.zern_tilts=zern_tilts#whether to use tilts during fit of zernikes.
        self.coeff=numpy.zeros((nzmax),numpy.float64)
        self.SRCDIR=os.getcwd()+"/"
        #now look for a shz file...
        files=os.listdir(self.SRCDIR)
        files.sort()
        filefound=[]
        filestart="shz%d_%d_"%(npup,nsubx)
        for file in files:
            if file[:len(filestart)]==filestart and file[-5:]==".fits":
                try:
                    tmpnzmax=int(file[len(filestart):-5])
                    if tmpnzmax>=nzmax:
                        filefound.append(tmpnzmax)
                except:
                    pass
        filefound.sort()

        if len(filefound)==0:#not os.path.exists(self.SRCDIR+"shz%d.fits"%nsubx):
            # (Re)calc and save subap averaged Zern gradients
            print "Computing zernike gradients"
            self.subzpist,self.subzxgrad,self.subzygrad=self.zgrads(npup,nsubx,nzmax,numpy.where(self.subarea==0,0.,1./self.subarea))
        else:# Read zernike pist and grads from file
            subgrad=util.FITS.Read(self.SRCDIR+"shz%d_%d_%d.fits"%(npup,nsubx,filefound[0]))[1]
            self.subzpist=subgrad[0]#numpy.zeros((nzmax,nsubx,nsubx),numpy.float64)
            self.subzxgrad=subgrad[1]#numpy.zeros((nzmax,nsubx,nsubx),numpy.float64)
            self.subzygrad=subgrad[2]#numpy.zeros((nzmax,nsubx,nsubx),numpy.float64)
        k=-2
        subzgrad=numpy.zeros((self.ndata,nzmax),numpy.float64)
        for i in range(nsubx): # Zernike SVD fit setup
            for j in range(nsubx):
                if(self.subflag[i][j]==1):
                    k+=2
                    for m in range(nzmax):
                        subzgrad[k,m]=   self.subzxgrad[m,i,j]
                        subzgrad[k+1,m]= self.subzygrad[m,i,j]
        print "Performing zernike setup (SVD/matrix inverse)"
        #status=cmod.zfit.setup(subzgrad) # C-call to do SVD setup (does svd decomposition into UWV, and sets small elements of W to zero.
        self.subzgrad=subzgrad
        if computeReconmx:
            self.zernikeReconMx,self.eigenVals=self.calcZernikeRecon(subzgrad,1e-6)
        else:
            self.zernikeReconMx=None
            self.eigenVals=None
        # save some multiplication later on in a loop...
        self.subzpist*=self.nsubx*numpy.pi
        self.subzxgrad*=2.*numpy.pi
        self.subzygrad*=2.*numpy.pi
    def calcZernikeRecon(self,pokemx,corr_thresh=None):
        """Make control matrix from poke matrix
        corr_thresh is the fractional value of the threshold,
        ie the threshold will be placed at corr_thresh*max value."""
        fail=0
        try:
            svec=numpy.linalg.svd(pokemx)
            #svec=LinearAlgebra.singular_value_decomposition(pokemx)
            u,a,vt=svec
        except:
            print "svd failed to converge..."
            fail=1
        if fail==0:
            n_removed=0
            thresh=max(a)*corr_thresh
            for i in range(len(a)):
                if a[i] < thresh:
                    a[i]=0.
                    n_removed += 1
            print 'Removed',n_removed,'modes from control matrix'
            #print 'Eigenvalues:',a
            ut=numpy.transpose(u)
            v=numpy.transpose(vt)
            id=numpy.identity(len(a))
            ai=numpy.multiply(a,id)
            ai=numpy.where(ai != 0, 1/ai, 0)
            #print v.shape,ai.shape,ut.shape
            reconmx = quick.dot(v, quick.dot(ai, ut[:ai.shape[0]]))
        else:#but this uses svd anyway, so wahts the point?
            reconmx=numpy.linalg.pinv(pokemx)
            a=None
        return reconmx,a

    def zgrads(self,npup,nsubx,nzmax,invarea):
        """Get mean zernike pist and grads over subaps + write to file
        invarea is 1/subarea, except ==0 where subarea==0.
        """
        n=npup/nsubx
        j=0                                                               # Make the Zernike fns 
        pp=[0]
        qq=[0]
        tt=[1]
        trig=0
        for p in range(1,self.nzrad+1):
            for q in range(p+1):
                if(numpy.fmod(p-q,2)==0):
                    if(q>0):
                        pp.append(p)
                        qq.append(q)
                        trig=not(trig)
                        tt.append(trig)
                        pp.append(p)
                        qq.append(q)
                        trig=not(trig)
                        tt.append(trig)
                    else:
                        pp.append(p)
                        qq.append(q)
                        tt.append(1)
                        trig=not(trig)
        gamx=None
        gamy=None
        #look for suitable gamma file...
        files=os.listdir(self.SRCDIR)
        filex=[]
        filey=[]
        for file in files:
            if file[:4]=="gamx" and file[-5:]==".fits":
                try:
                    nzradtmp=int(file[4:-5])
                    if nzradtmp>=nzrad:
                        filex.append(nzradtmp)
                except:
                    pass
            if file[:4]=="gamy" and file[-5:]==".fits":
                try:
                    nzradtmp=int(file[4:-5])
                    if nzradtmp>=nzrad:
                        filey.append(nzradtmp)
                except:
                    pass
        if len(filex)==0 or len(filey)==0:
            print "util.zernikeRecon - making gammas"
            gamx,gamy=util.zernikeMod.makegammas(self.nzrad,returntype="numpy")
        else:
            filex.sort()
            filey.sort()
            gamx=util.FITS.Read(self.SRCDIR+"gamx%d.fits"%filex[0])[1]
            gamy=util.FITS.REad(self.SRCDIR+"gamy%d.fits"%filey[0])[1]
        subgrad=numpy.zeros((3,nzmax,nsubx,nsubx),numpy.float64)
        subzpist=subgrad[0]#numpy.zeros((nzmax,nsubx,nsubx),numpy.float64)
        subzxgrad=subgrad[1]#numpy.zeros((nzmax,nsubx,nsubx),numpy.float64)
        subzygrad=subgrad[2]#numpy.zeros((nzmax,nsubx,nsubx),numpy.float64)
        ztmp=numpy.zeros((npup,npup),numpy.float64)
        tmp=numpy.zeros((nsubx,nsubx),numpy.float64)
        for m in range(nzmax):
            print "Creating zern %d/%d"%(m,nzmax-1)
            zern(ztmp,pp[m],qq[m],tt[m])
            #print "Binning..."
            cmod.binimg.binimg(ztmp,tmp)#bin the zernike...
            subzpist[m]=tmp*invarea#divide by the area...
            for i in range(nzmax):
                #cmod.binimg.binimg(gamx[i,m]*ztmp,tmp)
                subzxgrad[i]+=tmp*gamx[i,m]
                #cmod.binimg.binimg(gamy[i,m]*ztmp,tmp)
                subzygrad[i]+=tmp*gamy[i,m]
        for m in range(nzmax):
            subzxgrad[m]*=invarea
            subzygrad[m]*=invarea
        util.FITS.Write(subgrad,self.SRCDIR+'shz%d_%d_%d.fits'%(npup,nsubx,nzmax))
        return subzpist,subzxgrad,subzygrad

    def calc(self):
        """Reconstruct piston errors for current tilt data"""
        nsubx=self.nsubx
        ndata=self.ndata  
        subflag=self.subflag
        #first reorganise the centroids...
        idata=-1
        for i in range(nsubx):
            for j in range(nsubx):
                if(subflag[i][j]==1):
                    idata+=1
                    self.data[idata]=self.centx[i,j]
                    idata+=1
                    self.data[idata]=self.centy[i,j]
        #and then fit them.
        #cmod.zfit.fit(self.data,self.coeff)# C-call for SVD Zernike fit
        self.coeff[:,]=quick.dot(self.zernikeReconMx,self.data)
        self.pist*=0.
        self.xtilt*=0.
        self.ytilt*=0.
        for m in range(self.recon_jmin,min(self.nz+1,self.nzmax)):
            self.pist+=self.subzpist[m]*self.coeff[m]
        self.pist*=self.subflag#apply the subap mask.
        for i in range(nsubx):                                            # Pistons + tilt updates for DM:
            for j in range(nsubx):
                if(subflag[i][j]==1):
                    for m in range(self.recon_jmin,min(self.nz+1,self.nzmax)):
                        #self.pist[i][j]+=self.subzpist[m,i,j]*self.coeff[m]#*(float(self.wfs_nsubx)*numpy.pi)-now done at start...
                        if self.zern_tilts:    # Piston + tilts, or
                            self.xtilt[i][j]+=self.subzxgrad[m,i,j]*self.coeff[m]#*2.*numpy.pi-now done at start
                            self.ytilt[i][j]+=self.subzygrad[m,i,j]*self.coeff[m]#*2.*numpy.pi-now done at start
                    if not self.zern_tilts:  # Piston only, Tilts direct
                        self.xtilt[i,j]= self.centx[i,j]*2.*numpy.pi
                        self.ytilt[i,j]= self.centy[i,j]*2.*numpy.pi


    def dofit(self):                                           # Add new pistons and tilts to mirror
        """Add new pistons and tilts to mirror - ie fit the pistons with the
        tilts to the true mirror size."""
        nsubx=    self.nsubx
        n=    self.npup/self.nsubx
        phssub=numpy.zeros((n,n),numpy.float64)
        for i in range(nsubx):
            for j in range(nsubx):
                if(self.subflag[i,j]==1):
                    phssub=self.xtilt[i,j]*self.xtiltfn # Get segment phase fn
                    phssub+=self.ytilt[i,j]*self.ytiltfn  # for updated piston
                    phssub+=self.pist[i,j]                # and tilt
                    self.dmphs[i*n:(i+1)*n,j*n:(j+1)*n]=phssub # Tessalate onto overall DM phs map

    def run(self,centroidObject):
        """centroidObject is a util.centroid.centroid object, for which the
        run method has already been called..."""
        self.centx=centroidObject.centx
        self.centy=centroidObject.centy
        self.calc()
        self.dofit()
        return self
