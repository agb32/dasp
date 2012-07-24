#A module to compute phase covariance matricees (ie Noll matrix etc).
#Uses Kolmogorov and von Kalman specta.

#Also can be used to analyse turbulence, ie gets the zernike variance etc.

#The matrix computed here is defined as
#<a_i b_j> =
#\int dr_a dr_b <phi_a(Rr_a) phi_b(Rr_b)> Z_i(r_a) Z_j(r_b) W(r_a) W(r_b)
#where W is the aperture wieghting function, Z is the mirror modes
#(zernikes, or poked actuators etc) and phi is the phase.

#We can use Wilson and Jenkins eq B3 to compute the phi part...
#<phi phi> = 0.022883r_0^(-5/3)L_0^(8/3) x 0.841309 - D(R(r_b - r_a))
#Note that sigma2 is obtained by integrating the von Kalman spectrum (with
#a numerical integral for cos^(5/3)... giving the 0.841309 factor).
#sigma2 is the variance of the phase.
#After discussion with FA, turns out that 2sigma2 is equal to
#0.17253*(L0/r0)**(5./3), ie just the first part of the phase structure func.
#This is because of the definition of phase structure func at infinity...
import numpy as na
import numpy.fft
import scipy
import scipy.special#agbhome
#import cmod.besselkv#agbhome
import util.tel
import util.FITS
import util.zernikeMod
import cmod.phaseCov#agbhome
import types
import sys
#import FFT
#use the make() method...
def computeCov(modes,scale,r0,l0):
    """compute the covariance matrix.  Modes is a 3d array of [nmodes, x,y]
    and includes the (DM) pupil mask.
    Scale is the distance in m between 2 pixels on the modes matrix.
    """
    print "Use cmod.phaseCov.covariance() instead of this..."
    sigma2=0.022883*r0**(-5./3)*l0**(8./3)*0.841309
    xmax=modes.shape[2]
    ymax=modes.shape[1]
    nmodes=modes.shape[0]
    covMat=na.zeros((nmodes,nmodes),na.float64)
    for g in range(nmodes):
        for h in range(nmodes):
            for i in range(ymax):
                for j in range(xmax):
                    for k in range(ymax):
                        for l in range(xmax):
                            tmp=modes[g,i,k]*modes[h,k,l]
                            sep=scale*na.fabs((na.sqrt(i*i+j*j)-na.sqrt(k*k+l*l)))
                            tmp2=sigma2-vonKal(sep,l0,r0)
                            covMat[g,h]+=tmp2*tmp

    return covMat
                            

def computeCov2(modes,scale,r0,l0):
    """compute the covariance matrix.  Modes is a 3d array of [nmodes, x,y]
    and includes the (DM) pupil mask.
    Scale is the distance in m between 2 pixels on the modes matrix.
    """
    print "Use cmod.phaseCov.covariance() instead of this.  Need to create a phase array (vonKalArr or kolArr), the mirror modes array, and the output first."
    #sigma2=0.022883*r0**(-5./3)*l0**(8./3)*0.841309
    sigma2=0.5*0.17253*(l0/r0)**(5./3)
    xmax=modes.shape[2]
    ymax=modes.shape[1]
    nmodes=modes.shape[0]
    #first compute a von Kalman matrix...
    vk=na.zeros((ymax,xmax),na.float64)
    for y in range(ymax):
        for x in range(xmax):
            vk[y,x]=vonKal(scale*na.sqrt(x*x+y*y),l0,r0)


    covMat=na.zeros((nmodes,nmodes),na.float64)
    for g in range(nmodes):
        for h in range(nmodes):
            print g,h
            for i in range(ymax):
                for j in range(xmax):
                    partint=0.
                    for k in range(ymax):
                        for l in range(xmax):
                            #tmp=modes[h,k,l]
                            #sep=scale*na.fabs((na.sqrt(i*i+k*k)-na.sqrt(j*j+l*l)))
                            tmp2=sigma2-vk[int(abs(k-i)),int(na.fabs(l-j))]
                            partint+=tmp2*modes[h,k,l]
                    covMat[g,h]+=partint*modes[g,i,j]
    return covMat


def computeCov3(modes,nact,dmflag,r0,l0,telDiam):
    """
    Python version of quick algorith, which, if works, will be implemented in c
    """

    tmp=numpy.zeros((nact,nact),numpy.float32)
    mpos=0
    nmodes=dmflag.sum()
    out=numpy.zeros((nmodes,nmodes),numpy.float32)
    npup=modes.shape[1]
    pc=vonKalArr(r0,l0,npup,telDiam)
    import sys
    for g in range(nact):
        for h in range(nact):
            print "%d,%d    \r"%(g,h),
            sys.stdout.flush()
            if dmflag[g,h]:
                for i in range(npup):
                    for j in range(npup):
                        partint=0.
                        for k in range(npup):
                            for l in range(npup):
                                partint+=pc[numpy.abs(k-i),numpy.abs(l-j)]*modes[mpos,k,l]
                        tmp[g,h]+=partint*modes[mpos,i,j]
                mpos+=1
    mpos1=0
    mpos2=0
    for g1 in range(nact):
        for h1 in range(nact):
            if dmflag[g1,h1]:
                mpos2=0
                for g2 in range(nact):
                    for h2 in range(nact):
                        if dmflag[g2,h2]:
                            out[mpos1,mpos2]=tmp[numpy.abs(g1-g2),numpy.abs(h1-h2)]
                            mpos2+=1
                mpos1+=1
    return out

def computeCov4(mode,xcoord,ycoord,dmflag,modes,modeCoords,vig,npup,nact,r0,l0,telDiam,dovignetted=1):
    """
    Python version of quick algorith, which, if works, will be implemented in c
    First, takes 1 mode (assumes all identical), and computes covariance with itself in all positions.  Then this is used to create a BCCB or BTTB or something matrix.  Then any actuator influence functions which are vignetted in anyway are recalculated.
    Finally, unwanted parts are stripped out.
    mode is the local mode.  All modes must be the same when unvignetted for this to work (eg use pspline or something)
    x/ycoords are in pixels, the shift between modes.
    """
    tmp=numpy.zeros((nact,nact),numpy.float32)
    pc=vonKalArr(r0,l0,npup+mode.shape[0],telDiam*(npup+mode.shape[0])/float(npup))
    ny=mode.shape[0]
    nx=mode.shape[1]
    tmpmode=mode.copy()
    if ny!=nx:
        raise Exception("Non square mode")
    #print ycoord,xcoord
    for g in range(nact):
        for h in range(nact):
            yshift=ycoord[g]%1
            xshift=xcoord[h]%1
            print "%d,%d   shift %g %g \n"%(g,h,yshift,xshift),
            sys.stdout.flush()
            if yshift!=0 or xshift!=0:
                #non-integer pxl shift so interpolate the mode.
                #print "\nShifting by %g %g"%(yshift,xshift)
                mode2=tmpmode
                yin=numpy.arange(ny).astype(numpy.float32)
                xin=numpy.arange(nx).astype(numpy.float32)
                yout=yin-yshift
                xout=xin-xshift
                cmod.interp.gslCubSplineInterp_UB(mode,yin,xin,yout,xout,mode2,4)
            else:
                mode2=mode
            for i in range(ny):
                    for j in range(nx):
                        partint=0.
                        for k in range(ny):
                            for l in range(nx):
                                partint+=pc[int(numpy.abs(k-i+int(ycoord[g]))),int(numpy.abs(l-j+int(xcoord[h])))]*mode2[k,l]
                        tmp[g,h]+=partint*mode[i,j]
            #tmp[h,g]=tmp[g,h]

    #Now ravel this and copy it...
    #tmp=tmp.ravel()
    nacts=nact*nact
    phasecov=numpy.zeros((nact,nact,nact,nact),numpy.float32)
    #phasecov[0]=tmp
    for i in range(nact):
        for j in range(nact):
            for k in range(nact):
                for l in range(nact):
                    phasecov[i,j,k,l]=tmp[numpy.abs(i-k),numpy.abs(j-l)]

    phasecov.shape=(nacts,nacts)
    #and now strip out modes that are unused...
    dmflag=dmflag.ravel()
    nmodes=dmflag.sum()
    phs=numpy.zeros((nmodes,nmodes),numpy.float32)
    pos=0
    for i in range(nacts):
        if dmflag[i]:
            if i!=pos:
                phasecov[pos]=phasecov[i]
                phasecov[:,pos]=phasecov[:,i]
            pos+=1
    phasecov=phasecov[:nmodes,:nmodes]
    #And now redo the modes that are vignetted...
    if dovignetted:
        makeWithLocalModes(npup,modes,modeCoords,r0=r0,l0=l0,telDiam=telDiam,typ="vk",nthreads=8,lam=500.,vig=vig,output=phasecov)
    return phasecov

def computeCov5(mode,xcoord,ycoord,dmflag,modes,modeCoords,vig,npup,nact,r0,l0,telDiam,dovignetted=1,nthreads=1):
    """As computeCov4, but using a c module"""
    tmp=numpy.zeros((nact,nact),numpy.float32)
    pc=vonKalArr(r0,l0,npup+mode.shape[0],telDiam*(npup+mode.shape[0])/float(npup))
    if not mode.flags.contiguous:
        mode=mode.copy()
    ny=mode.shape[0]
    nx=mode.shape[1]
    if ny!=nx:
        raise Exception("Non square mode")
    xyCoords=numpy.zeros((2,nact),numpy.float32)
    xyCoords[0]=ycoord
    xyCoords[1]=xcoord
    phasecov=numpy.zeros((nact,nact,nact,nact),numpy.float32)
    cmod.phaseCov.covarianceQuick(pc,mode,xyCoords,tmp,phasecov,nthreads)
    #Now ravel this and copy it...
    #tmp=tmp.ravel()
    nacts=nact*nact
    #phasecov[0]=tmp
    #for i in range(nact):
    #    for j in range(nact):
    #        for k in range(nact):
    #            for l in range(nact):
    #                phasecov[i,j,k,l]=tmp[numpy.abs(i-k),numpy.abs(j-l)]

    phasecov.shape=(nacts,nacts)
    #and now strip out modes that are unused...
    dmflag=dmflag.ravel()
    nmodes=dmflag.sum()
    phs=numpy.zeros((nmodes,nmodes),numpy.float32)
    pos=0
    for i in range(nacts):
        if dmflag[i]:
            if i!=pos:
                phasecov[pos]=phasecov[i]
                phasecov[:,pos]=phasecov[:,i]
            pos+=1
    phasecov=phasecov[:nmodes,:nmodes]
    #And now redo the modes that are vignetted...
    if dovignetted:
        makeWithLocalModes(npup,modes,modeCoords,r0=r0,l0=l0,telDiam=telDiam,typ="vk",nthreads=nthreads,lam=500.,vig=vig,output=phasecov)
    return phasecov



def vonKal(sep,l0,r0):
    """compute von Kalman spectrum at separation sep (in m).
    From Jenkins 1998, MNRAS 294, the spatial structure function."""
    bessel=scipy.special.kv#Bessel function of the second kind.
    #bessel=cmod.besselkv.besselkv#agbhome
    gamma=1.1287870299081262#scipy.special.gamma(5./6)
    if sep==0:
        return 0.
    return 0.17253*(l0/r0)**(5./3)*(1-2*(na.pi*sep/l0)**(5./6)/gamma*float(bessel(5./6,2*na.pi*sep/l0)))


def vonKalArr(r0=0.2,l0=30.,npup=32,telDiam=4.,cx=0,cy=0):
    """scale is telDiam/npup.
    This computes the von kalman phase.
    """
    scale=float(telDiam)/npup
    #sigma2=0.022883*r0**(-5./3)*l0**(8./3)*0.841309
    sigma2=0.5*0.17253*(l0/r0)**(5./3)
    
    vk=na.zeros((npup,npup),na.float32)
    for y in range(npup):
        for x in range(npup):
            r=float(scale*na.sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy)))
            vk[y,x]=sigma2-vonKal(r,l0,r0)/2
    return vk

def kolArr(r0=0.2,npup=32,telDiam=4.,pupilfn=None):
    """This computes the kolmogorov phase, minus the piston (see RWW paper).
    This takes a while because of the double integral (20secs for 32x32 pupil function), while the vonKalman version is very fast.  So probably best not to use this!!!  If you can be bothered, create a c version - should be less than a second then."""
    scale=float(telDiam)/npup
    vk=na.zeros((npup,npup),na.float32)
    pfn=pupilfn/na.sum(na.sum(pupilfn))

    for y in range(npup):
        for x in range(npup):
            r=na.sqrt(x*x+y*y)*scale
            vk[y,x]=-0.5*6.88*(r/r0)**(5./3)
            tmp=0.
            for i in range(npup):
                for j in range(npup):
                    a=i-npup/2
                    b=j-npup/2
                    r=scale*na.sqrt(a*a+b*b)
                    rs=(r/r0)**(5./3)
                    tmp+=pfn[i,j]*rs
                    a=i-y-npup/2
                    b=j-y-npup/2
                    r=scale*na.sqrt(a*a+b*b)
                    rs=(r/r0)**(5./3)
                    tmp+=pfn[i,j]*rs
            vk[y,x]=vk[y,x]+0.5*6.88*tmp
    print "doing double integral"
    tmp=0.
    for i in range(npup):
        for j in range(npup):
            for ii in range(npup):
                for jj in range(npup):
                    r=scale*na.sqrt((i-ii)*(i-ii)+(j-jj)*(j-jj))
                    rs=(r/r0)**(5./3)
                    tmp+=pfn[i,j]*pfn[ii,jj]*rs
    vk-=na.array(0.5*6.88*tmp,na.float32)
    return vk

def version():
    """Here, return a version string, which should be altered if any of these functions here or in the c module are changed - this means that the saved ones can no longer be used...
    For version 2, changed sigma2 definition."""
    return "2"

def make(typ="vk",npup=32,nmode=16,pupil=None,modes=None,r0=0.2,telDiam=4.,l0=30.,vk=None,nthreads=8,scale=0,lam=500.):
    """Compute covariance matrix, either for vonKalman spectrum or kolmogorov.
    Modes is defined on the pupil, such that orthogonal modes (note - your basis may not be orthogonal) will sum to one when multiplied together.
    """
    if type(pupil)==types.NoneType:
        pupil=util.tel.Pupil(npup,npup/2,npup/8)#use pupil.fn
    if type(modes)==types.NoneType:
        zerns=util.zernikeMod.Zernike(pupil,nmode)#use zerns.zern
        modes=na.array(zerns.zern/na.sqrt(na.sum(na.sum(pupil.fn)))).astype(na.float32)
    if type(pupil)!=na.ndarray:
        pupil=pupil.fn
    if type(vk)==types.NoneType:
        if typ=="vk":
            vk=vonKalArr(r0,l0,npup,telDiam)
        else:
            vk=kolArr(r0,npup,telDiam,pupilfn=pupil)
    output=na.zeros((nmode,nmode),na.float32)
    cmod.phaseCov.covariance(vk,modes,output,nthreads)
    if scale:
        output/=pupil.sum()
        #s=int(na.sum(na.sum(pupil.fn)))
        output*=4.4862630963/5.111101472#/s#these factors here are to place in agreement with the noll matrix when zernikes are used. (after the noll matrix has been scaled by (D/ro)**(5/3))  No idea where these factors come from...
    output*=(500./lam)**2#wavelength scaling
    return pupil,modes,vk,output

def makeWithLocalModes(npup,modes,modeCoords,r0=0.12,l0=30.,telDiam=4.2,typ="vk",nthreads=8,lam=500.,vig=None,output=None):
    """Same as make(), but uses mirror modes that are local to each actuator, and not the whole mirror... these are created using makeLocalMirrorModes in util.dm.
    modeCoords is also returned by makeLocalMirrorModes().
    if vig!=None, must be a 1d int32 array of flags for subaps to do.
    if vig!=None, you will also want to specify output... which should be partly filled already...
    """
    if typ=="vk":
        vk=vonKalArr(r0,l0,npup,telDiam)
    else:
        vk=kolArr(r0,npup,telDiam,pupilfn=pupil)
    nmode=modes.shape[0]
    if output==None:
        output=na.zeros((nmode,nmode),na.float32)
    cmod.phaseCov.covarianceLocal(vk,modes,modeCoords,vig,output,nthreads)
    output*=(500./lam)**2
    return output

def makeWithLocalModesFFT(npup,modes,modeCoords,r0=0.12,l0=30.,telDiam=4.2,typ=
"vk",nthreads=8,lam=500.,vig=None,output=None):
    """Attempts to be a faster version using FFTs for convolution.
    """
    nmodes=modes.shape[0]
    ny=modes.shape[1]
    nx=modes.shape[2]
    if ny!=nx:
        raise Exception("util.phaseCovariance: ny==nx not valid")
    if typ=="vk":
        vk=vonKalArr(r0,l0,npup*2,telDiam*2,cx=npup,cy=npup)
        #vk=vonKalArr(r0,l0,npup+ny,telDiam*(npup+ny)/float(npup),cx=nx,cy=ny)
    else:
        vk=kolArr(r0,npup,telDiam,pupilfn=pupil)
        raise Exception("not yet implemented")

    if output==None:
        output=numpy.zeros((nmodes,nmodes),numpy.float32)
    m2=numpy.zeros((ny*2,nx*2),numpy.float32)
    print "doing loop"
    for g in xrange(nmodes):
        gx=modeCoords[g,0]
        gy=modeCoords[g,1]
        for h in xrange(nmodes):
            hx=modeCoords[h,0]
            hy=modeCoords[h,1]
            #xmin=numpy.abs(hx-gx)#modeCoords[h,0]
            #ymin=numpy.abs(hy-gy)#modeCoords[h,1]
            xmin=gx-hx+npup-nx
            ymin=gy-hy+npup-ny
            xmax=xmin+nx*2
            ymax=ymin+ny*2
            vkpart=vk[ymin:ymax,xmin:xmax]
            m2[ny/2:ny*3./2,nx/2:nx*3./2]=modes[h]#zeropad the mode.
            tmp=numpy.fft.ifft2(numpy.fft.fft2(vkpart)*numpy.fft.fft2(numpy.fft.fftshift(m2))).real[ny/2:ny*3/2.,nx/2:nx*3/2.]
            output[g,h]=(tmp*modes[g]).sum()
    output*=(500./lam)**2
    return output

def modeFFTWorker(gstart,gend,modes,fftModes):
    print "Starting modeFFT for %d->%d"%(gstart,gend)
    nmodes=modes.shape[0]
    ny=modes.shape[1]
    nx=modes.shape[2]
    m2=numpy.zeros((ny*2,nx*2),"f")
    for g in xrange(gstart,gend):
        m2[ny/2:ny*3./2,nx/2:nx*3./2]=modes[g]#zeropad the mode.
        fftModes[g]=numpy.fft.fft2(m2)

def localFFTWorker(gstart,gend,modes,modeCoords,fftModes,npup,vk,output):
    import sys
    nmodes=modes.shape[0]
    ny=modes.shape[1]
    nx=modes.shape[2]
    #m2=numpy.zeros((ny*2,nx*2),"f")
    print "Starting for %d->%d (%g elements)"%(gstart,gend,(2*nmodes-gstart-gend)/2.*(gend-gstart))
    for g in xrange(gstart,gend):
        gx=modeCoords[g,0]
        gy=modeCoords[g,1]
        if gstart==0:
            print "%d / %d (%d)   \r"%(g,gend,nmodes),
            sys.stdout.flush()
        for h in xrange(g,nmodes):
            #if gstart==0:
            #    print "%d %d / %d %d   \r"%(g,h,gend,nmodes),
            #    sys.stdout.flush()
            hx=modeCoords[h,0]
            hy=modeCoords[h,1]
            #xmin=numpy.abs(hx-gx)#modeCoords[h,0]
            #ymin=numpy.abs(hy-gy)#modeCoords[h,1]
            xmin=gx-hx+npup-nx
            ymin=gy-hy+npup-ny
            xmax=xmin+nx*2
            ymax=ymin+ny*2
            vkpart=vk[ymin:ymax,xmin:xmax]
            #m2[ny/2:ny*3./2,nx/2:nx*3./2]=modes[h]#zeropad the mode.
            #fftModes=numpy.fft.fft2(m2)
            tmp=numpy.fft.ifft2(numpy.fft.fft2(numpy.fft.fftshift(vkpart))*fftModes[h]).real[ny/2:ny*3/2.,nx/2:nx*3/2.]
            output[g,h]=(tmp*modes[g]).sum()
            output[h,g]=output[g,h]#symmetric.
def makeWithLocalModesFFTThreaded(npup,modes,modeCoords,r0=0.12,l0=30.,telDiam=4.2,typ="vk",nthreads=8,lam=500.,vig=None,output=None):
    """Attempts to be a faster version using FFTs for convolution.
    """
    import os
    import processing
    nmodes=modes.shape[0]
    ny=modes.shape[1]
    nx=modes.shape[2]
    if ny!=nx:
        raise Exception("util.phaseCovariance: ny==nx not valid")
    if typ=="vk":
        vk=vonKalArr(r0,l0,npup*2,telDiam*2,cx=npup,cy=npup)
        #vk=vonKalArr(r0,l0,npup+ny,telDiam*(npup+ny)/float(npup),cx=nx,cy=ny)
    else:
        vk=kolArr(r0,npup,telDiam,pupilfn=pupil)
        raise Exception("not yet implemented")

    if output==None:
        output=numpy.memmap("/dev/shm/phasecovoutput",dtype="f",mode="w+",shape=(nmodes,nmodes))
        os.unlink("/dev/shm/phasecovoutput")
    fftModes=numpy.memmap("/dev/shm/fftModes",dtype="F",mode="w+",shape=(nmodes,ny*2,nx*2))
    gend=0
    plist=[]
    for i in range(nthreads):
        gstart=gend
        gend=gstart+(nmodes-gend)/(nthreads-i)
        if i==nthreads-1:
            gend=nmodes
        plist.append(processing.Process(target=modeFFTWorker,args=(gstart,gend,modes,fftModes)))
        plist[-1].start()
    for p in plist:
        p.join()
    gend=0
    plist=[]
    NperThread2=(nmodes*nmodes)/nthreads
    for i in range(nthreads):
        gstart=gend
        gend=int(nmodes+0.5-numpy.sqrt(nmodes*nmodes-NperThread2+gstart*gstart-2*nmodes*gstart))#gstart+(nmodes-gend)/(nthreads-i)
        if i==nthreads-1:
            gend=nmodes
        plist.append(processing.Process(target=localFFTWorker,args=(gstart,gend,modes,modeCoords,fftModes,npup,vk,output)))
        plist[-1].start()
    for p in plist:
        p.join()
    out=numpy.array(output)*(500./lam)**2
    return out
        



def computeCov3local(modes,modeCoords,npup,nact,dmflag,r0,l0,telDiam):
    """
    Python version of quick algorith, which, if works, will be implemented in c
    """

    tmp=numpy.zeros((nact,nact),numpy.float32)
    mpos=0
    nmodes=dmflag.sum()
    out=numpy.zeros((nmodes,nmodes),numpy.float32)
    pc=vonKalArr(r0,l0,npup,telDiam)
    ny=modes.shape[1]
    nx=modes.shape[2]
    print ny,nx
    import sys
    for g in range(nact):
        goffx=modeCoords[mpos,0]
        goffy=modeCoords[mpos,1]
        for h in range(nact):
            print "%d,%d    \r"%(g,h),
            sys.stdout.flush()
            if dmflag[g,h]:
                hoffx=modeCoords[mpos,0]
                hoffy=modeCoords[mpos,1]
                for i in range(ny):
                    for j in range(nx):
                        partint=0.
                        for k in range(ny):
                            for l in range(nx):
                                try:
                                    t1=pc[numpy.abs(k+hoffy-i-goffy),numpy.abs(l+hoffx-j-goffx)]
                                except:
                                    print k,hoffy,i,goffy,l,hoffx,j,goffx
                                    raise
                                partint+=t1*modes[mpos,k,l]
                        tmp[g,h]+=partint*modes[mpos,i,j]
                mpos+=1
    mpos1=0
    mpos2=0
    for g1 in range(nact):
        for h1 in range(nact):
            if dmflag[g1,h1]:
                mpos2=0
                for g2 in range(nact):
                    for h2 in range(nact):
                        if dmflag[g2,h2]:
                            print "%d %d %d %d %d %d"%(g1,h1,mpos1,g2,h2,mpos2)
                            out[mpos1,mpos2]=tmp[numpy.abs(g1-g2),numpy.abs(h1-h2)]
                            mpos2+=1
                mpos1+=1
    return out

    

def makeNoll(file="/home/ali/cvsstuff/aosim/util/noll.fits",telDiam=4.,r0=0.2,nfft=8,wfsn=8,nsubx=8):
    noll=util.FITS.Read(file)[1]*(telDiam/r0)**(5./3.)
    scale=float(nfft)/wfsn/nsubx/na.pi
    return noll,scale
    

def toeplitz():
    """compute phase covariance as almost a block toeplitz toeplitz block
    matrix"""
    import util.tel
    pupil=util.tel.Pupil(64,128,0)
    import util.fitModes
    actModes=na.zeros((81,9,9),na.float32)
    mirrorModes=na.zeros((81,64,64),na.float32)
    for i in range(81):
        y=i/9
        x=i%9
        actModes[i,y,x]=1.
        mirrorModes[i]=util.fitModes.interp(actModes[i],"bicubic",64,0.1,0.25,9)
    p,m,v,o=make(typ="vk",npup=64,nmode=81,pupil=pupil,modes=mirrorModes,r0=0.2,telDiam=4.2,l0=30.)
    return o
def bttb(nact=9,telDiam=4.2,r0=0.2,l0=30.):
    """
    Block Toeplitz with Toeplitz Blocks.
    compute as a proper block toeplitz matrix.
    This uses a much simpler model for the modes, ie just assuming 1 actuator
    is pushed, and ignoring interpolation etc.
    Approximate phase covariance as a true BTTB matrix.
    The inverse of this matrix is sparse diagonal (ie a few diagonal lines).
    Completely ignores pupil function.
    """
    sigma2=0.5*0.17253*(l0/r0)**(5./3)
    nmode=nact*nact
    actSep=telDiam/(nact-1)#actuator separation
    output=na.zeros((nmode,nmode),na.float32)
    for i in range(nact):#iterate over blocks
        block=output[:nact,i*nact:i*nact+nact]
        for j in range(nact):#iterate within a block
            sep=na.sqrt((j*actSep)**2+(i*actSep)**2)
            block[0,j]=sigma2-vonKal(sep,l0,r0)/2
        for j in range(nact):#copy to rest of block
            block[j,j:,]=block[0,:nact-j]
            block[j:,j]=block[0,:nact-j]
    for i in range(nact):#copy block to rest of blocks
        output[i*nact:i*nact+nact,i*nact:,]=output[:nact,:nact*nact-i*nact]
        output[i*nact:,i*nact:i*nact+nact]=na.transpose(output[:nact,:nact*nact-i*nact])
    return output

def bccb(nact=9,telDiam=4.2,r0=0.2,l0=30.,expand=1):
    """
    Block Circulant with Circulant Blocks.
    Compute phase covariance approximation as block circulant, circulant
    block matrix.
    This looks to be a fairly poor approximation.
    Note, the inverse of this matrix (ie the inverse noise covariance matrix)
    could possibly be treated as diagonal... infact, just the identity matrix
    scaled.
    I think the result is in radians (squared?), so may need to convert this to
    pixels, depending on your application.
    If expand==0, return only the first row...
    """
    sigma2=0.5*0.17253*(l0/r0)**(5./3)
    nmode=nact*nact
    #if nact%2!=1:#seems to work...
        #print "Warning - nact not odd - may cause problems in util.phaseCovariance.bccb (untested)"
    actSep=telDiam/(nact-1)#actuator separation
    if expand:
        output=na.zeros((nmode,nmode),na.float32)
        for i in range(nact/2+1):#iterate over blocks
            block=output[:nact,i*nact:i*nact+nact]
            for j in range(nact/2+1):#iterate within a block
                sep=na.sqrt((j*actSep)**2+(i*actSep)**2)
                block[0,j]=sigma2-vonKal(sep,l0,r0)/2
                if j>0:#make block circulant matrix
                    block[0,nact-j]=block[0,j]
            for j in range(nact):#copy to rest of block
                block[j,j:,]=block[0,:nact-j]
                block[j:,j]=block[0,:nact-j]
            if i>0:#make block circulant
                output[:nact,nact*(nact-i):nact*(nact-i+1)]=block
        for i in range(nact):#copy block to rest of blocks
            output[i*nact:i*nact+nact,i*nact:,]=output[:nact,:nact*nact-i*nact]
            output[i*nact:,i*nact:i*nact+nact]=na.transpose(output[:nact,:nact*nact-i*nact])
    else:#generate first row only.
        output=na.zeros((nmode,),na.float32)
        for i in range(nact/2+1):#iterate over blocks
            block=output[i*nact:(i+1)*nact]
            for j in range(nact/2+1):#iterate within a block
                sep=na.sqrt((j*actSep)**2+(i*actSep)**2)
                block[j]=sigma2-vonKal(sep,l0,r0)/2
                if j>0:#make block circulant
                    block[nact-j]=block[j]
            #copy to rest of blocks
            if i>0:
                output[(nact-i)*nact:(nact-i+1)*nact]=block
    return output
    

def checkCirculant(m):
    """Checks to see whether a matrix is circulant"""
    a=expandToCirculant(m[0])
    s=na.sum((a==m).astype("l").ravel())
    if s==m.shape[0]*m.shape[1]:
        return 1
    else:
        return 0

def expandToCirculant(m):
    """Expand the first row of a circulant matrix (a) to a full circulant matrix."""
    a=na.zeros((m.shape[0],m.shape[0]),m.dtype)
    a[0]=m
    #expand a into a circulant matrix...
    for i in range(1,a.shape[0]):
        a[i,i:,]=a[0,:a.shape[1]-i]
        a[i,:i]=a[0,a.shape[1]-i:,]
    return a
    
def invertCirculantMatrix(c):

    """inverts a circulant matrix (note, not BCCB).  This is
    equivalent to LinearAlgebra.inverse(c)[0], but much faster.  See F. P
    Pijpers (1999) Unbiased image reconstruction as an inverse problem
    Monthly Notices of the Royal Astronomical Society 307 (3),
    659-668. for more info.

    Note, this just returns the first row of the circulant matrix -
    other rows are obtained by wrapping this row round.
    """
    
    #ident=na.identity(c.shape[0],na.int32)
    if len(c.shape)==2:
        cc=c[0]
    else:#just the first row has been passed...
        cc=c
    return numpy.fft.ifft(1/numpy.fft.fft(cc)).real

def invertBCCB(c,nb=None):
    """inverts a BCCB matrix, returning the first row only.
    This can then be copied and wrapped to give other rows.
    nb is the number of blocks...
    This is the same as LinearAlgebra.inverse(c)[0]
    """
    ntot=c.shape[0]
    if nb==None:
        nb=int(na.sqrt(ntot))#number of blocks
    ne=int(ntot/float(nb))#number of elements in a block.
    if ne!=nb:
        print "Warning phaseCovariance.invertBCCB not tested for nblocks!=nelements..."
    if len(c.shape)==2:
        m=c[0]
    else:
        m=c
    m.shape=ne,nb
    #ident=na.zeros((ne,nb))
    #ident[0,0]=1
    rt=na.reshape(numpy.fft.ifft2(1/numpy.fft.fft2(m)).real,(ntot,))
    m.shape=(ntot,)
    return rt

def multBCCB(a,b,nb=None):
    """multiplies a and b, two BCCB matricees.
    nb is the number of blocks.
    This does the same as na.dot(a,b)[0]
    """
    if nb==None:
        nb=int(na.sqrt(a.shape[0]))
    ne=int(a.shape[0]/float(nb))
    if ne!=nb:
        print "Warning phaseCovariance.multBCCB not tested for nblocks!=nelements..."
    aa=a[0]
    aa.shape=ne,nb
    bb=b[0]
    bb.shape=ne,nb
    return na.reshape(numpy.ifft2(numpy.fft.fft2(aa)*numpy.fft.fft2(bb)).real,(a.shape[0],))
    
    
def testPhaseCov(npup,nact,actoffset,pupil,dmObj,dmidstr,mirrorModes=None,nthreads=8,dmInterpType=None,actCoupling=None,actFlattening=None):
    """follow the sequence used in xinterp_recon.
    First, make mirror modes, and then use them to make the phase covariance.
    Can run it using (in aosim.eagle/baseComparison):
import base.readConfig
config=base.readConfig.AOXml("paramsForPhaseCovariance.xml",batchno=0)
v=config.getVal
import util.phaseCovariance
res=util.phaseCovariance.testPhaseCov(v("npup"),v("nAct"),0.,config.getVal("pupil"),config.getVal("dmObj"),"0,")

    """
    import util.dm
    if dmInterpType==None:
        dmInterpType="spline"
        print "Warning - util.phaseCovariance.testPhaseCov assuming spline"
    actmode=na.zeros((nact,nact),numpy.float32)
    mirrorSurface=util.dm.MirrorSurface(dmInterpType,npup,nact,actoffset=actoffset,actCoupling=actCoupling,actFlattening=actFlattening)
    dmflag=dmObj.computeDMPupil(dmidstr,centObscuration=pupil.r2,retPupil=0)[0]
    dmindices=na.nonzero(dmflag.ravel())[0]
    nmodes=dmindices.shape[0]
    mirrorScale=None
    if mirrorModes==None:
        mirrorModes=na.zeros((nmodes,npup,npup),numpy.float32)
        mirrorScale=na.zeros((nmodes,),numpy.float32)

        for i in range(nmodes):
            y=dmindices[i]/(nact)
            x=dmindices[i]%(nact)
            actmode[y,x]=1
            mirrorModes[i]=mirrorSurface.fit(actmode)
            mirrorScale[i]=na.sqrt(na.sum(mirrorModes[i]*mirrorModes[i]))
            mirrorModes[i]/=mirrorScale[i]
            actmode[y,x]=0
    p,m,v,o=make(typ="vk",npup=npup,nmode=nmodes,pupil=pupil,modes=mirrorModes,r0=0.12,telDiam=4.2,l0=30.,nthreads=nthreads)
    return mirrorModes,mirrorScale,o
    
def run(nact=9,npup=None,width=-1,nthreads=8,dovignetted=1):
    """Runs a test... used during development, and for timing."""
    import util.atmos
    import util.dm
    import time
    t1=time.time()
    if npup==None:
        npup=(nact-1)*8
    dm=util.dm.dmInfo("dm",["a"],0,nact,interpType="pspline",reconLam=500.)
    atmosGeom=util.atmos.geom({"L0":util.atmos.layer(0,0.,10.,1.,10)},[util.atmos.source("a",0,0,-1,8)],npup,npup,4.2,0.12,30)
    dmObj=util.dm.dmOverview([dm],atmosGeom)
    dm.computeDMPupil(atmosGeom,retPupil=0)
    r2=0.
    phasecov=dm.computePhaseCovariance(atmosGeom,r2,width=width,nthreads=nthreads,dovignetted=dovignetted)
    t2=time.time()
    print "\nTook %gs"%(t2-t1)
    return phasecov,dmObj,atmosGeom,t2-t1

if __name__=="__main__":
    tlist=[]
    tlist2=[]
    tlist3=[]
    alist=[17,33,49,65,85]#8.15,99,362,907,2176,
    #for i in alist:
    #    pc=run(i,(i-1)*8)
    #    tlist.append(pc[3])
    #print "Actuators calculated:",alist
    #print "Time taken with 8x8 pxls per act:",tlist
    #for i in alist:#300,3818,14806,36194,92310
    #    pc=run(i,(i-1)*20)
    #    tlist2.append(pc[3])
    #print "Actuators calculated:",alist
    #print "Time taken with 8x8 pxls per act:",tlist
    #print "Time taken with 20x20 pxls per act:",tlist2
    print "*** testfftThreaded ***"
    for i in alist:#37,447,2043,6526,18905s
        pc=run(i,(i-1)*20,width=-1)
        tlist3.append(pc[3])
    print "Actuators calculated:",alist
    print "Time taken with 20x20 pxls per act using testfftThreaded:",tlist3#210,3096
    #Time taken with 20x20 pxls per act using localThreaded: [21.783562898635864, 236.8443238735199, 1141.2177050113678, 3330.3107211589813, 9762.2584209442139]
    #Now doing fftModes stuff, takes: 19, 203, 953, 2771, 8030
