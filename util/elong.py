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

    def plotProfile(self,n,alt,depth,plotDepth):
        """n is number of points to plot.  plotDepth is the depth (centred on alt) of which to display.
        alt/depth are the mean alt/depth."""
        x=(numpy.arange(n)-n/2.)*plotDepth/n+alt
        res=numpy.zeros((n,),numpy.float32)
        for i in range(len(self.peakList)):
            off,d,intensity=self.peakList[i]
            xx=x-alt-off
            dd=d*depth
            res+=intensity*numpy.exp(-(xx**2)/(2*dd**2))
        return res,x



class LineProfile:
    def __init__(self,n,nsubx,centralHeight,profile,heights,pxlscale,telDiam,launchDist=0.,launchTheta=0.,unelong=1.,oversamplefactor=1):
        """profile is the Na layer strengths, at heights.
        xoff,yoff can be used to shift the spots by this many pixels.
        pxlscale is arcsec per pixel.
        centralHeight is the central height of the profile, in m.
        profile is the profile strength in TBC
        heights are the heights at which profile values are valid.
        n is the number of pixels per sub-aperture required for this LGS correlation image.
        nsubx is the number of sub-apertures.
        launchDist in m
        launchTheta in degrees
        unelong is the unelongated width in arcsec.
        oversamplefactor is a factor by which to oversample during generation, with binning by this factor before return to the user (of size n,n)

        To include a zenith angle, just change the profile heights and centralHeight (I think)
        """
        launchTheta=-launchTheta
        self.n=n
        self.nsubx=nsubx
        self.centralHeight=centralHeight
        self.profile=profile
        self.heights=heights
        self.xoff=0.
        self.yoff=0.
        self.pxlscale=pxlscale
        self.telDiam=telDiam
        self.launchTheta=launchTheta*numpy.pi/180.
        self.launchDist=launchDist
        self.unelong=unelong
        self.oversamplefactor=oversamplefactor
        self.no=self.n*self.oversamplefactor
        self.psf=util.centroid.createAiryDisc2(self.no,unelong/self.pxlscale*oversamplefactor,0,-1)
    def generate(self,clearPadding=0):
        """How this works:
        Integrates along the line profile at required resolution.
        Puts as a single line in a 2D array.
        Convolves with a PSF.
        Rotates to where it should be.

        LGS fires up to centralHeight being directly above the telescope (ie if launchDist!=0, plume isn't vertical).

        If clearPadding is set, and oversamplefactor>1, then it will zero anything around the padding area.

        """
        xlaunch=self.launchDist*numpy.cos(self.launchTheta)
        ylaunch=self.launchDist*numpy.sin(self.launchTheta)
        #How many pixels along the spot? Depends on rotation (launchTheta)
        maxspotlen=int(numpy.ceil(self.no*(1+(numpy.sqrt(2)-1)*numpy.sin(self.launchTheta*2)**2)))
        opxlscale=self.pxlscale/self.oversamplefactor
        tmp=numpy.zeros((maxspotlen,maxspotlen),numpy.float32)
        distFromLaunchToCentralHeight=numpy.sqrt(self.launchDist**2+self.centralHeight**2)
        plumeangle=numpy.arctan2(self.launchDist,self.centralHeight)
        img=numpy.zeros((self.nsubx,self.nsubx,self.no,self.no),numpy.float32)
        #tmpimg=numpy.zeros((self.nsubx,self.nsubx,maxspotlen,maxspotlen),"f")
        import scipy.signal
        import scipy.ndimage
        for y in range(self.nsubx):
            if self.nsubx>=20:
                print "LGS PSF: %d/%d"%(y+1,self.nsubx)
            for x in range(self.nsubx):
                #compute dist/angle from launch.
                yy=(y-self.nsubx/2.+0.5)/self.nsubx*self.telDiam
                xx=(x-self.nsubx/2.+0.5)/self.nsubx*self.telDiam
                r=numpy.sqrt(xx*xx+yy*yy)#distance of subap to centre in m.
                #compute distances from centre of subap:
                distToLaunch=numpy.sqrt((xlaunch-xx)**2+(ylaunch-yy)**2)
                distToCentralHeight=numpy.sqrt(r**2+self.centralHeight**2)
                #Now compute the angle subtended by the beam and from centre of subap to beam.  Cosine rule:
                psi=numpy.arccos((distToLaunch**2-distToCentralHeight**2-distFromLaunchToCentralHeight**2)/(-2*distFromLaunchToCentralHeight*distToCentralHeight))
                #psi=numpy.arccos(((xlaunch-xx)**2+(ylaunch-yy)**2-xx**2-yy**2-self.centralHeight**2-self.launchDist**2-self.centralHeight**2)/(-2*distFromLaunchToCentralHeight*distToCentralHeight))
                #angle of the plume within the subap.
                subapTheta=numpy.arctan2(ylaunch-yy,xlaunch-xx)
                subapR=r
                #centre pixel (maxspotlen/2) is at centralHeight.
                for i in range(maxspotlen):
                    #starting and ending heights for this pixel (with 0 being central height)
                    anglestart=opxlscale*(i-maxspotlen/2.)/3600.*numpy.pi/180.
                    angleend=anglestart+(opxlscale/3600.*numpy.pi/180.)
                    #which corresponds to these distances in m along the beam:
                    if anglestart>0:
                        phi=numpy.pi-psi-anglestart#180 degrees in a triangle
                    else:
                        #anglestart=-anglestart
                        phi=psi+anglestart
                    diststart=distToCentralHeight*numpy.sin(anglestart)/numpy.sin(phi)+self.centralHeight
                    if angleend>0:
                        phi=numpy.pi-psi-angleend
                    else:
                        #angleend=-angleend
                        phi=psi+angleend
                    distend=distToCentralHeight*numpy.sin(angleend)/numpy.sin(phi)+self.centralHeight
                    #if i<maxspotlen/2. and distend>self.centralHeight:
                    #    tmp[maxspotlen/2,maxspotlen-1-i]=0
                    #    continue
                    if diststart>distend:
                        diststart=0#I think this can happen if beyond infinity... so just set to zero - since distend will be small in this case.
                    #if (y==3 or y==4) and x==3:
                    #    print "%d %d %d %g %g %g %g %g %g"%(y,x,i,anglestart,angleend,phi,psi,diststart,distend)
                    #and which correspond to these profile heights:
                    diststart=diststart*numpy.cos(plumeangle)
                    distend=distend*numpy.cos(plumeangle)
                    #Now integrate the profile between these heights to get the flux in this pixel
                    tmp[(maxspotlen)/2,maxspotlen-1-i]=self.integrate(diststart,distend)
                #Now convolve with the psf
                res=scipy.signal.convolve2d(tmp,self.psf,mode="same")
                #tmpimg[y,x]=tmp
                #now rotate
                rot=scipy.ndimage.interpolation.rotate(res,-subapTheta*180/numpy.pi)
                sh=rot.shape
                ys=(sh[0]-self.no)/2
                xs=(sh[1]-self.no)/2
                #and select the rotated part
                img[y,x]=rot[ys:ys+self.no,xs:xs+self.no]
        img=numpy.where(img<0,0,img)
        if clearPadding and self.oversamplefactor>1:
            s=(self.no-self.n)/2
            e=s+self.n
            img[:,:,:s]=0
            img[:,:,e:]=0
            img[:,:,:,:s]=0
            img[:,:,:,e:]=0
        return img#,tmpimg

    def integrate(self,h1,h2):
        """integrate the profile from h1 to h2, using trapeziums.
        """
        i=0
        p=self.profile
        h=self.heights
        sig=0.
        while i<h.size and h[i]<h1:
            i+=1
        if i==self.heights.size:
            return 0.#h1 is above all the profile so no flux.
        if i==0:
            #print "Warning: no profile specified at <=%gm"%h1
            pass
            #assume it to be 0.
        else:#compute profile at h1
            i-=1
            p1=p[i]+(p[i+1]-p[i])/(h[i+1]-h[i])*(h1-h[i])
            #now integrate to h[i+1]
            sig=(p1+p[i+1])/2.*(h[i+1]-h1)
            i+=1
        #now integrate up to h2
        while i<h.size-1 and h[i]<h2:
            sig+=(p[i]+p[i+1])/2.*(h[i+1]-h[i])
            i+=1
        #and now do final bit
        if i==h.size-1 or i==0:
            #print "Warning: no profile specified at >=%gm"%h2
            pass
        else:
            p2=p[i-1]+(p[i]-p[i-1])/(h[i]-h[i-1])*(h2-h[i-1])
            sig+=(p2+p[i-1])/2.*(h2-h[i-1])
        if sig<0:
            print "OOOps - sig %g"%sig
        return sig

    def fnOrig(self,i,j):
        """Called for pixel i,j of a given sub-aperture (for which the relevant self.... are set)."""
        y=(i-self.no/2.)+0.5-self.xoff*self.oversamplefactor
        x=(j-self.no/2.)+0.5-self.yoff*self.oversamplefactor
        #Parallel to the profile, interpolate.
        #Perpendicular to the profile, decay exponentially.
        #Beyond the end of the profile, decay exponentially.
        #1.  How far from the profile axis are we?
        r=numpy.sqrt(x*x+y*y)#in pixels.
        phi=numpy.arctan2(y,x)
        subtendAngle=self.lineangle-numpy.pi/2-phi
        r_perp=r*numpy.cos(subtendAngle)
        r_along=r*numpy.sin(subtendAngle)
        sigma=self.unelong/self.pxlscale*self.oversamplefactor
        sigsq=2*sigma**2
        #decay further away from the line profile.  Need to integrate across the pixel I think.
        if self.subapLaunchDist==0:#subap is under launch location
            scale=numpy.exp(-r**2/sigsq)/(2*self.oversamplefactor)**2
        else:
            scale=numpy.exp(-r_perp**2/sigsq)

        #Now work out what height we are at (using r_along).
        theta=r_along/self.oversamplefactor*self.pxlscale/3600./180.*numpy.pi
        #0 arcsec corresponds to central height.
        #subapLaunchDist is the launchDist from this subap.
        #thetaOffset is elevation angle up to central height for this subap.  i.e. is arctan(centralHeight/subapLaunchDist.
            
        h=self.subapLaunchDist*numpy.tan(theta+self.thetaOffset)
        #Now we know the height, interpolate along the profile.
        if self.subapLaunchDist==0:#subap is under launch location
            strength=self.profile.sum()
            print "DIST:",strength,scale
        elif h<self.heights[0]:
            strength=0
        elif h>self.heights[-1]:
            strength=0
        else:
            strength=numpy.interp(h,self.heights,self.profile)
        return strength*scale

    def generateOrig(self):
        xlaunch=self.launchDist*numpy.cos(self.launchTheta)
        ylaunch=self.launchDist*numpy.sin(self.launchTheta)
        img=numpy.zeros((self.nsubx,self.nsubx,self.no,self.no),numpy.float32)
        for y in range(self.nsubx):
            for x in range(self.nsubx):
                #compute dist/angle from launch.
                yy=(y-self.nsubx/2.+0.5)/self.nsubx*self.telDiam
                xx=(x-self.nsubx/2.+0.5)/self.nsubx*self.telDiam
                r=numpy.sqrt(xx*xx+yy*yy)#distance of subap to centre in m.
                self.subapTheta=numpy.arctan2(yy,xx)
                self.subapR=r
                yyl=yy+ylaunch
                xxl=xx+xlaunch
                rr=numpy.sqrt(xxl*xxl+yyl*yyl)#dist of subap from launch in m.
                self.lineangle=numpy.arctan2(yyl,xxl)
                self.subapLaunchDist=rr
                self.thetaOffset=numpy.arctan2(self.centralHeight,self.subapLaunchDist)
                #self.lineangle=self.launchTheta-numpy.arcsin(r/rr*numpy.sin(self.subapTheta-self.launchTheta))##angle of spot across subap.
                print y,x,self.lineangle*180./numpy.pi,self.launchTheta*180/numpy.pi,r,rr,self.subapTheta*180/numpy.pi
                for i in range(self.no):
                    for j in range(self.no):
                        img[y,x,i,j]=self.fnOrig(i,j)
                #img[y,x]=numpy.fromfunction(self.fn,(self.n,self.n))
        if self.oversamplefactor!=1:
            img.shape=self.nsubx,self.nsubx,self.n,self.oversamplefactor,self.n,self.oversamplefactor
            img=img.sum(5).sum(3)
        return img
    def tile(self,img=None):
        if img==None:
            img=self.generate()
        nsubx=img.shape[0]
        n=img.shape[2]
        out=numpy.zeros((nsubx*n,nsubx*n),numpy.float32)
        for i in range(nsubx):
            for j in range(nsubx):
                out[i*n:i*n+n,j*n:j*n+n]=img[i,j]
        return out

    def plotProfile(self,n,alt,depth,plotDepth):
        """n is number of points to plot.  plotDepth is the depth (centred on alt) of which to display.
        alt/depth are the mean alt/depth."""
        x=(numpy.arange(n)-n/2.)*plotDepth/n+alt
        res=numpy.zeros((n,),numpy.float32)
        for i in range(len(self.peakList)):
            off,d,intensity=self.peakList[i]
            xx=x-alt-off
            dd=d*depth
            res+=intensity*numpy.exp(-(xx**2)/(2*dd**2))
        return res,x


def make(spotsize=32,nsubx=110,wfs_n=16,wfs_nfft=None,wfs_nimg=None,clipsize=None,lam=640e-9,telDiam=42.,telSec=None,beacon_alt=90000.,beacon_depth=10000.,zenith=0.,photons=1e6,unelong_width=1.,fname=None,launchDist=0.,launchTheta=0.,xoff=0.,yoff=0.,pup=None):
    """
    unelong_width in arcsec is the unelongated spot width.
    zenith is in degrees.
    For a scaled model, scale nsubx, telDiam, beacon_alt and beacon_depth.
    launchDist(m) and launchTheta(deg, anti clockwise from 3 oclock(!)) specify the position of the launch telescope.  If launchDist==0, this is onaxis launch.
    Was anticlock from 12 oclock - changed 2013/10/11 to standard position (3).

    """
    launchTheta-=90#move from 12 oclock to 3 oclock
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
    if pup is None:
        print "Generating pupil function for elongated spot intensities"
        pup=util.tel.Pupil(nsubx*spotsize,nsubx*spotsize/2,nsubx*spotsize/2*telSec/telDiam).fn
    else:
        if hasattr(pup,"fn"):
            if pup.fn.shape[0]!=nsubx*spotsize or pup.shape[1]!=nsubx*spotsize:
                print "Rescaling pupil function for LGS elongation pattern"
                pup=pup.rescale(nsubx*spotsize)
            pup=pup.fn
    if pup.shape[0]!=nsubx*spotsize or pup.shape[1]!=nsubx*spotsize:
        raise Exception("Error in elong - pup shape should be equal to the lgs psf shape, %d,%d"%(nsubx*spotsize,nsubx*spotsize))
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
            a1=numpy.arctan2(h1,r)
            a2=numpy.arctan2(h2,r)
            da=(a2-a1)*3600*180/numpy.pi
            e=numpy.sqrt(da*da+unelong_width*unelong_width)
            
            psfs[isub,jsub]=gauss.make(unelong_width/pixscale,e/pixscale,phi)
            psfs[isub,jsub]*=photons/psfs[isub,jsub].sum()*pup[isub*spotsize:(isub+1)*spotsize,jsub*spotsize:(jsub+1)*spotsize].sum()/spotsize**2
            
    img=numpy.zeros((nsubx*spotsize,nsubx*spotsize),numpy.float32)
    for i in range(nsubx):
        for j in range(nsubx):
            img[i*spotsize:(i+1)*spotsize,j*spotsize:(j+1)*spotsize]=psfs[i,j]
    if fname!=None:
        util.FITS.Write(psfs,fname)
        util.FITS.Write(img,fname,writeMode="a")  #human-viewable version
    return psfs,img

def make2(spotsize=32,nsubx=110,wfs_n=16,wfs_nfft=None,wfs_nimg=None,clipsize=None,lam=640e-9,telDiam=42.,telSec=None,beacon_alt=90000.,beacon_depth=10000.,lineProfile=None,zenith=0.,photons=1e6,unelong_width=1.,fname=None,launchDist=0.,launchTheta=0.):
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
            xsub=(float(isub)-float(nsubx)/2.+0.5)*telDiam/float(nsubx)-launchDist*numpy.cos(launchTheta/180.*numpy.pi)
            ysub=(float(jsub)-float(nsubx)/2.+0.5)*telDiam/float(nsubx)-launchDist*numpy.sin(launchTheta/180.*numpy.pi)
            r=numpy.sqrt(xsub*xsub+ysub*ysub)
            phi=numpy.arctan2(xsub,ysub)*180./numpy.pi
            
            h1=(beacon_alt-beacon_depth/2.)/numpy.cos(zenith*numpy.pi/180.)
            h2=(beacon_alt+beacon_depth/2.)/numpy.cos(zenith*numpy.pi/180.)
            a1=numpy.arctan2(r,h1)
            a2=numpy.arctan2(r,h2)
            a=numpy.arctan2(r,(beacon_alt/numpy.cos(zenith*numpy.pi/180.)))*180./numpy.pi
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

