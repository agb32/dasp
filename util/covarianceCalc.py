"""This is to be used after many instances of covariance.py have run..."""
import os
import sys
import numpy

import util.FITS



#kill=[0,1, 7,8,  9,17,   63,71,  72,73, 79,80]
# Pupil function
def pupil(npup,r1=-1.,r2=-1.):
    if r1 < 0.:
        r1 = float(npup/2)
    pup=numpy.zeros((npup,npup),numpy.int)
    for i in range(npup):
        for j in range(npup):
            r=numpy.sqrt((float(i-npup/2)+0.5)**2+(float(j-npup/2)+0.5)**2)
            #r=sqrt((float(i-npup/2))**2+(float(j-npup/2))**2)
            if ((r<=r1)&(r>=r2)):
                pup[i,j]=1
    return pup



if len(sys.argv)!=4:
    print "Usage: %s lam directory outname"%sys.argv[0]
    sys.exit(0)
lam=float(sys.argv[1])
directory=sys.argv[2]
outname=sys.argv[3]
print "Using %g nm"%lam

# Read in data
ls=os.listdir(directory)#commands.getoutput('ls cov2').split()
ndata=len(ls)
print 'Found',ndata,'files'

tmp=util.FITS.Read(directory+"/"+ls[0])[1]
print tmp.shape
nact=tmp.shape[0]
n=int(numpy.sqrt(nact))

actdata=numpy.zeros((ndata,nact),numpy.float)
for i in range(ndata):
	actdata[i]=util.FITS.Read(directory+"/"+ls[i])[1]
#util.FITS.Write(actdata.reshape((ndata,n,n)),'actdata.fits')


# Subtract means
actmean=actdata.sum(0)/float(ndata)
#util.FITS.Write(actmean.reshape((n,n)),'actmean.fits')
actdata-=actmean
#util.FITS.Write(actdata.reshape((ndata,n,n)),'actdata2.fits')

#actmean2=actdata.sum(0)/float(ndata)
#util.FITS.Write(actmean2.reshape((n,n)),'actmea2n.fits')


# Calculate variance
var=(actdata**2).sum(0)/float(ndata)


# Calculate covariance
cov=numpy.zeros((nact,nact),numpy.float)
for j in range(ndata):
	tmp=actdata[j]	
	for i in range(nact):
		cov[i]+=tmp*tmp[i]
cov/=float(ndata)
cov/=(lam/500.)**2#This is the correct thing to do...

# Write data
#util.FITS.Write(var.reshape((n,n)),'var.fits')
util.FITS.Write(cov,outname)
print "Written to %s"%outname

#os.system('open var.fits')
#os.system('open cov.fits')



