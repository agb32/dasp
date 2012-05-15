#$Id: dist.py,v 1.6 2008/03/05 15:43:09 ali Exp $
"""
Module dist

Equivalent of dist IDL function written in Python

Author : F Assemat
"""

import numpy as na
#import Numeric as na#if swapping between numpy/Numeric, remember to sort out nonzero().

def dist(n,m=None,dy=0,dx=0,natype=na.float64,sqrt=1):
    """
    The DIST function creates a rectangular array in which the value of
    each element is proportional to its frequency. This array may be used
    for a variety of purposes, including frequency-domain filtering and
    making pretty pictures.

    Syntax
        Result=dist(n,m=None,dy=0,dx=0,natype=na.float32):


Arguments
         n  The number of lines in the resulting array.

         m  The number of columns in the resulting array. If M is omitted,
            the resulting array will be n by n

         dy The offset in y of the center (0 by default)
         dx The offset in x of the center (0 by default)
         natype : the Numarray type of the return array (Float64 by default)
         
"""
#dimensions of the output array
    if m is None:
        m=n
        
    #Generation of the different axis a grid of horizontales lines
    axe_x=na.arange(m,dtype=natype)
    #axe_x.savespace(1)
    axe_x-=m/2+dx
    axe_y=na.arange(n,dtype=natype)
    #axe_y.savespace(1)
    axe_y-=n/2+dy

    #creation of the grid of square distances
    grid=axe_x[na.newaxis,:,]**2+axe_y[:,na.newaxis]**2
    if sqrt:
        na.sqrt(grid,grid)

    return grid

def profilRadial(psf,radiusOnly=0,nfft=None,dtype=None):
    """Computes the radial profile of the 2D array psf
    The output array contains two arrays : 
    - the 1st array stores the radius (in pixels)
    - the 2nd array stores the radial profile
    An FA function.
    """	
    ####  We first create a map of the square of the distances to the center of the PSF image
    ####  The map of distances is computed with the dist function (see aosim/utils/dist.py for help)
    ####  We use the square because we have then integer indexes
    if type(psf)!=type(None):
        nfft=psf.shape[0]
        try:
            dtype=psf.typecode()
        except:
            dtype=psf.dtype.char
    r2=dist(nfft);r2*=r2 ## square : we have integer numbers
    nnRad=na.argsort(r2.flat) ##we flatten the grid of distances and sort the values
    r2r=na.take(r2.flat,nnRad) ## we sort the grid of distances
    xRad=na.nonzero(na.diff(r2r))[0] ##we look when the grid of distances change of value
    #### Allocation of the return result
    result=na.zeros((2,xRad.shape[0]),dtype=dtype)
    rRad=na.take(na.sqrt(r2r),xRad) ##radius (in pixels) giving the horizontal axis for radial and encircled energy profiles
    #### We store the radial axis in line 0 of the returned array
    result[0]=rRad.astype(dtype)
    if radiusOnly:
        return result[0]

    tabNbPointsRad=na.diff(xRad) ## number of points per bin
    #### We flatten the PSF and sort the values
    psfSort=na.take(psf.flat,nnRad)
    
    #### We compute the encircled energy of the PSF
    tabEncircledEnergy=na.take(na.cumsum(psfSort),xRad)
    #### We compute the radial profile
    tabEnergyPerPixBin=na.diff(tabEncircledEnergy) ##to have the energy per pixel bin
    profil=tabEnergyPerPixBin/tabNbPointsRad
    
    
    #### We store the radial profile in line 1 of the returned array
    result[1,0]=psfSort[0]
    result[1,1:,]=profil[:,].astype(dtype) ##radial profile of PSF
    
    return result
