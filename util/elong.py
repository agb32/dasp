import numpy
import util.FITS
import util.calcPxlScale
import util.tel


class Gauss:
    def __init__(self,n,peakList=None,xoff=0.,yoff=0.):
        """peakList is tuple of (alt,depth,strength)
        xoff,yoff can be used to shift the spots by this many pixels.
        """
        self.n=n
        self.peakList=peakList
        self.xoff=xoff
        self.yoff=yoff
    def fn(self,i,j):
        #x=(i-self.n/2).astype(numpy.float32)+0.5
        #y=(j-self.n/2).astype(numpy.float32)+0.5
        x=(i-self.n/2).astype(numpy.float)+0.5-self.xoff
        y=(j-self.n/2).astype(numpy.float)+0.5-self.yoff
        tmp=numpy.exp(-1.*(self.c*y*y + 2.*self.b*x*y + self.a*x*x))
        if numpy.any(numpy.isnan(tmp)):
            print "nan %d %d"%(i,j)
        return tmp
        #return self.c*y*y + 2.*self.b*x*y + self.a*x*x


    def make(self,fwhmx,fwhmy,th):
        sig2x=(fwhmx**2)/(8.*numpy.log(2.))
        sig2y=(fwhmy**2)/(8.*numpy.log(2.))
        
        thr=th*numpy.pi/180.
        self.a=(numpy.cos(thr)**2)/(2.*sig2x) + (numpy.sin(thr)**2)/(2.*sig2y)
        self.b=(numpy.sin(2.*thr))/(4.*sig2y) - (numpy.sin(2.*thr))/(4.*sig2x)
        self.c=(numpy.sin(thr)**2)/(2.*sig2x) + (numpy.cos(thr)**2)/(2.*sig2y)

        return numpy.fromfunction(self.fn,(self.n,self.n))

class MultiGauss:
    def __init__(self,n,peakList=[]):
        """Make a multiple gaussian peak, with height h1, depth d1, intensity (relative) i1.
        peakList is an array of tuples of (h,d,i)
        h is relative to the beacon_alt (so can be negative).
        If h==0, means this gaussian is centred.
        h is in metres.
        d is a fraction of beacon_depth.
        i is the fractional intensity.
        """
        self.n=n
        self.peakList=peakList

    def fn(self,i,j):
        #x=(i-self.n/2).astype(numpy.float32)+0.5
        #y=(j-self.n/2).astype(numpy.float32)+0.5
        x=(i-self.n/2).astype(numpy.float)+0.5-self.xoff
        y=(j-self.n/2).astype(numpy.float)+0.5-self.yoff
        tmp=numpy.exp(-1.*(self.c*y*y + 2.*self.b*x*y + self.a*x*x))
        if numpy.any(numpy.isnan(tmp)):
            print "nan %d %d"%(i,j)
        #actually should probably do an integration here... esp for narrow peaks.
        return tmp
        #return self.c*y*y + 2.*self.b*x*y + self.a*x*x

    def make(self,fwhmx,fwhmy,theta,phi,beacon_alt,pxlscale,zenith):
        """fwhm are in pixels, theta in degrees, phi in degrees,zenith in radians"""
        if zenith!=0:
            print "Warning - non-zero zenith not yet worked out..."
        img=numpy.zeros((self.n,self.n),numpy.float64)
        thr=theta*numpy.pi/180.
        phi*=numpy.pi/180.
        for i in range(len(self.peakList)):
            off,d,intensity=self.peakList[i]
            #need to scale off to get correct results...
            #pxloff=numpy.arctan(off*numpy.sin(numpy.arctan(r/h))/h)*3600*180/numpy.pi/pxlscale
            h=(beacon_alt-off)#/numpy.cos(zenith)
            r=beacon_alt*numpy.tan(phi)#/numpy.cos(zenith)
            pxloff=numpy.arctan(off/h*r/h)*3600*180/numpy.pi/pxlscale
            sig2x=(fwhmx**2)/(8.*numpy.log(2.))#*d**2
            sig2y=(fwhmy**2)/(8.*numpy.log(2.))*(d/numpy.cos(zenith))**2
        
            # compute the rotation components... and include the widths.
            self.a=(numpy.cos(thr)**2)/(2.*sig2x) + (numpy.sin(thr)**2)/(2.*sig2y)
            self.b=(numpy.sin(2.*thr))/(4.*sig2y) - (numpy.sin(2.*thr))/(4.*sig2x)
            self.c=(numpy.sin(thr)**2)/(2.*sig2x) + (numpy.cos(thr)**2)/(2.*sig2y)
            #and compute the offset...
            #self.xoff=off*numpy.sin(thr)*numpy.sin(phi)
            #self.yoff=off*numpy.cos(thr)*numpy.sin(phi)
            self.xoff=pxloff*numpy.sin(thr)
            self.yoff=pxloff*numpy.cos(thr)
            #print self.xoff,off*numpy.sin(thr)*numpy.sin(phi)
            img+=intensity*numpy.fromfunction(self.fn,(self.n,self.n))
        return img



def make(spotsize=32,nsubx=110,wfs_n=16,wfs_nfft=None,wfs_nimg=None,clipsize=None,lam=640e-9,telDiam=42.,telSec=None,beacon_alt=90000.,beacon_depth=10000.,zenith=0.,photons=1e6,unelong_width=1.,fname=None,launchDist=0.,launchTheta=0.,xoff=0.,yoff=0.):
    """
    unelong_width in arcsec is the unelongated spot width.
    zenith is in degrees.
    For a scaled model, scale nsubx, telDiam, beacon_alt and beacon_depth.
    launchDist(m) and launchTheta(deg) specify the position of the launch telescope.  If launchDist==0, this is onaxis launch.
    """
    if wfs_nfft==None:
        wfs_nfft=wfs_n*2
    if wfs_nimg==None:
        wfs_nimg=wfs_n
    if clipsize==None:
        clipsize=wfs_nfft
    if telSec==None:
        telSec=telDiam/7.
    binfactor=clipsize/wfs_nimg
    pixscale=util.calcPxlScale.pxlScale(lam,float(telDiam)/nsubx,wfs_n,wfs_nfft,clipsize/wfs_nimg)/binfactor#0.05 #.1   #arcsec
    #We want the pixel scale of the spots as generated, not after use.  Hence, the division by binfactor here.
    print "pixscale: %g arcsec/pxl during fft, or %g arcsec/pxl on CCD."%(pixscale,pixscale*binfactor)
    pup=util.tel.Pupil(nsubx*spotsize,nsubx*spotsize/2,nsubx*spotsize/2*telSec/telDiam).fn
    gauss=Gauss(spotsize,xoff=xoff,yoff=yoff)
    # test=gauss.make(unelong_width/pixscale,10.,45.)
    # print test,numpy.isnan(test),numpy.any(numpy.isnan(test))
    # util.FITS.Write(test,'test.fits')

    # Make the elongated PSFs
    psfs=numpy.zeros((nsubx,nsubx,spotsize,spotsize),numpy.float32)
    lx=launchDist*numpy.cos(launchTheta/180*numpy.pi)
    ly=launchDist*numpy.sin(launchTheta/180*numpy.pi)

    for isub in range(nsubx):
        for jsub in range(nsubx):
            # print isub,jsub
            xsub=(float(isub)-float(nsubx)/2.+0.5)*telDiam/float(nsubx)+lx
            ysub=(float(jsub)-float(nsubx)/2.+0.5)*telDiam/float(nsubx)+ly
            r=numpy.sqrt(xsub*xsub+ysub*ysub)
            phi=numpy.arctan2(xsub,ysub)*180./numpy.pi
            
            h1=(beacon_alt-beacon_depth/2.)/numpy.cos(zenith*numpy.pi/180.)
            h2=(beacon_alt+beacon_depth/2.)/numpy.cos(zenith*numpy.pi/180.)
            a1=numpy.arctan(h1/r)
            a2=numpy.arctan(h2/r)
            da=(a2-a1)*3600*180/numpy.pi
            e=numpy.sqrt(da*da+unelong_width*unelong_width)
            
            psfs[isub,jsub]=gauss.make(unelong_width/pixscale,e/pixscale,phi)
            psfs[isub,jsub]*=photons/psfs[isub,jsub].sum()*pup[isub*spotsize:(isub+1)*spotsize,jsub*spotsize:(jsub+1)*spotsize].sum()/spotsize**2
            
    img=numpy.zeros((nsubx*spotsize,nsubx*spotsize),numpy.float)
    for i in range(nsubx):
        for j in range(nsubx):
            img[i*spotsize:(i+1)*spotsize,j*spotsize:(j+1)*spotsize]=psfs[i,j]
    if fname!=None:
        util.FITS.Write(psfs,fname)
        util.FITS.Write(img,fname,writeMode="a")  #human-viewable version
    return psfs,img

def make2(spotsize=32,nsubx=110,wfs_n=16,wfs_nfft=None,wfs_nimg=None,clipsize=None,lam=640e-9,telDiam=42.,telSec=None,beacon_alt=90000.,beacon_depth=10000.,lineProfile=None,zenith=0.,photons=1e6,unelong_width=1.,fname=None):
    """
    unelong_width in arcsec is the unelongated spot width.
    zenith is in degrees.
    For a scaled model, scale nsubx, telDiam, beacon_alt and beacon_depth.
    lineProfile is an object which describes the LGS line profile.  Typically a gaussian or double peaked gaussian.
    """
    if wfs_nfft==None:
        wfs_nfft=wfs_n*2
    if wfs_nimg==None:
        wfs_nimg=wfs_n
    if clipsize==None:
        clipsize=wfs_nfft
    if telSec==None:
        telSec=telDiam/7.
    if lineProfile==None:
        lineProfile=Gauss(spotsize)
    binfactor=clipsize/wfs_nimg
    pixscale=util.calcPxlScale.pxlScale(lam,float(telDiam)/nsubx,wfs_n,wfs_nfft,clipsize/wfs_nimg)/binfactor#0.05 #.1   #arcsec
    print "pixscale: %g arcsec/pxl during fft, or %g arcsec/pxl on CCD."%(pixscale,pixscale*binfactor)
    pup=util.tel.Pupil(nsubx*spotsize,nsubx*spotsize/2,nsubx*spotsize/2*telSec/telDiam).fn
    # Make the elongated PSFs
    psfs=numpy.zeros((nsubx,nsubx,spotsize,spotsize),numpy.float32)
    for isub in range(nsubx):
        for jsub in range(nsubx):
            # print isub,jsub
            xsub=(float(isub)-float(nsubx)/2.+0.5)*telDiam/float(nsubx)
            ysub=(float(jsub)-float(nsubx)/2.+0.5)*telDiam/float(nsubx)
            r=numpy.sqrt(xsub*xsub+ysub*ysub)
            phi=numpy.arctan2(xsub,ysub)*180./numpy.pi
            
            h1=(beacon_alt-beacon_depth/2.)/numpy.cos(zenith*numpy.pi/180.)
            h2=(beacon_alt+beacon_depth/2.)/numpy.cos(zenith*numpy.pi/180.)
            a1=numpy.arctan(r/h1)
            a2=numpy.arctan(r/h2)
            a=numpy.arctan(r/(beacon_alt/numpy.cos(zenith*numpy.pi/180.)))*180./numpy.pi
            da=(a2-a1)*3600*180/numpy.pi
            e=numpy.sqrt(da*da+unelong_width*unelong_width)#elongated width
            
            psfs[isub,jsub]=lineProfile.make(unelong_width/pixscale,e/pixscale,phi,a,beacon_alt,pixscale,zenith*numpy.pi/180.)
            psfs[isub,jsub]*=photons/psfs[isub,jsub].sum()*pup[isub*spotsize:(isub+1)*spotsize,jsub*spotsize:(jsub+1)*spotsize].sum()/spotsize**2
            
    img=numpy.zeros((nsubx*spotsize,nsubx*spotsize),numpy.float)
    for i in range(nsubx):
        for j in range(nsubx):
            img[i*spotsize:(i+1)*spotsize,j*spotsize:(j+1)*spotsize]=psfs[i,j]
    if fname!=None:
        util.FITS.Write(psfs,fname)
        util.FITS.Write(img,fname,writeMode="a")  #human-viewable version
    return psfs,img


if __name__=="__main__":
    make()

