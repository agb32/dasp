"""pxlscale=lambda/d*wfs_n/nfft*binfactor*180/numpy.pi*3600
where binfactor= clipsize/nimg
Note, clipsize can be > fftsize, if using a psf.
FoV is then equal to nfft/binfactor*pixel scale - this is the fov over which you can see spot motion.  This needs to be sufficient to be bigger than size of diffraction spot plus spot motion.
Note, this isn't equal to nimg*pixelscale, which could be much larger...
Diam is the subap diameter.
"""
import numpy
def pxlScale(lam,diam,n,nfft,binfactor):
    if lam>1:
        lam*=1e-9
    return lam/diam*n/nfft*binfactor*180*3600/numpy.pi

def spotMotionRMS(lam,r0,d,r0IsAt500nm=1):
    """calculate expected spot motion for kolmogorov turb for square subap (angle of arrival variance) taken from David St Jacques thesis, at r0 specified at wavelength lambda.  Note, r0 proportional to lambda**(6./5)
    """
    if lam>1:
        lam*=1e-9
    if r0IsAt500nm:
        r0*=(lam/500e-9)**(6./5)
        print "Scaling r0 from 500nm to %gnm gives %g"%(lam*1e9,r0)
    return numpy.sqrt(0.162*lam**2*r0**(-5./3)*d**(-1./3)*(3600*180/numpy.pi)**2)

"""Seeing/seeing and turbulence:

http://www.astrosurf.com/cavadore/optique/turbulence/

r0 defines observed seeing:  seeing = 251*lambda/r0 with lambda in um, r0 in mm.
In arcsecs.


Or, 0.0251*lambda/r0 with lambda in nm, r0 in cm.

eg:
0.251 * 500 / 12.9   = .97 arcsec seeing.

According to Eric, this factor should be 0.2057
"""
