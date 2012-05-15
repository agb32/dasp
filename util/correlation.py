import numpy
"""This code is taken from the rtc code.  It was essentially used for testing of the RTC correlation code - checking that the half conjugate stuff was working okay.  However, it is here, for this, and also because the transformPSF function is used too.  This function has been modified to simulation format.
"""

def transformPSF(psf):
    """Function to transform psf into a form that will be usable by the RTC.
    psf is eg the LGS spot elongation pattern.
    psf should have shape nsubx,nsubx,nimg,nimg. or nimg,nimg.
    """
    if psf==None:
        return None
    psf=psf.copy()
    if len(psf.shape)==4:
        nsubx=psf.shape[0]
        for i in range(nsubx):
            for j in range(nsubx):
                ff=r2hc(numpy.fft.fftshift(psf[i,j]))
                # do the conjugate.
                mm=ff.shape[0]/2+1
                nn=ff.shape[1]/2+1
                ff[:mm,nn:]*=-1
                ff[mm:,:nn]*=-1
                # and now put it into the result.
                psf[i,j]=ff
    elif len(psf.shape)==2:#2D psf.
        ff=r2hc(numpy.fft.fftshift(psf))
        # do the conjugate.
        mm=ff.shape[0]/2+1
        nn=ff.shape[1]/2+1
        ff[:mm,nn:]*=-1
        ff[mm:,:nn]*=-1
        # and now put it into the result.
        psf[:]=ff
    else:
        raise Exception("PSF shape wrong in util/correlation.py")
    return psf.astype(numpy.float32)
        

def r2hc(a):
    """FFT r to hc.
    """
    a=a.copy()
    for i in range(a.shape[0]):
        a[i]=r2hc1d(a[i])
    for i in range(a.shape[1]):
        a[:,i]=r2hc1d(a[:,i])
    return a


def r2hc1d(a):
    """1D version
    """
    res=numpy.fft.fft(a)
    b=a.copy()
    n=a.shape[0]
    b[:n/2+1]=res.real[:n/2+1]
    b[n/2+1:]=res.imag[(n+1)/2-1:0:-1]
    #b[n/2+1:]=res.imag[1:n/2]
    return b

def hc2r1d(a):
    b=numpy.zeros(a.shape,numpy.complex64)
    n=b.shape[0]
    b.real[:n/2+1]=a[:n/2+1]
    b.real[n/2+1:]=a[(n+1)/2-1:0:-1]
    b.imag[(n+1)/2-1:0:-1]=a[n/2+1:]
    b.imag[(n+0)/2+1:]=-a[(n+0)/2+1:]
    res=numpy.fft.ifft(b).real
    return res


def hc2r(a):
    """inverse FFT hc to r
    """
    a=a.copy()
    for i in range(a.shape[0]):
        a[i]=hc2r1d(a[i])
    for i in range(a.shape[1]):#I dont think it matters which axis you go over first.
        a[:,i]=hc2r1d(a[:,i])
    return a

def unpack(a):
    """a is the output from r2hc.  Here, we expand into a standard complex array, so that by inverting this, can test you understand what is going on.
    """
    res=numpy.zeros(a.shape,numpy.complex64)
    m=a.shape[0]
    n=a.shape[1]
    for i in range(1,(m+1)/2):
        res[i,0]=a[i,0]+1j*a[m-i,0]
        res[m-i,0]=a[i,0]-1j*a[m-i,0]
        if n%2==0:
            res[i,n/2]=a[i,n/2]+1j*a[m-i,n/2]
            res[m-i,n/2]=a[i,n/2]-1j*a[m-i,n/2]
        

        for j in range(1,(n+1)/2):
            if i==1:
                res[0,j]=a[0,j]+1j*a[0,n-j]
                res[0,n-j]=a[0,j]-1j*a[0,n-j]
                if m%2==0:
                    res[m/2,j]=a[m/2,j]+1j*a[m/2,n-j]
                    res[m/2,n-j]=a[m/2,j]-1j*a[m/2,n-j]
            res[i,j]=a[i,j]-a[m-i,n-j]+1j*(a[i,n-j]+a[m-i,j])
            res[m-i,n-j]=a[i,j]-a[m-i,n-j]-1j*(a[i,n-j]+a[m-i,j])
            res[i,n-j]=a[i,j]+a[m-i,n-j]+1j*(a[m-i,j]-a[i,n-j])
            res[m-i,j]=a[i,j]+a[m-i,n-j]-1j*(a[m-i,j]-a[i,n-j])
    res[0,0]=a[0,0]+0j
    if n%2==0:
        res[0,n/2]=a[0,n/2]+0j
    else:
        res[0,n/2]=a[0,n/2]+1j*a[0,n/2+1]
        
    if m%2==0:
        res[m/2,0]=a[m/2,0]+0j
    else:
        res[n/2,0]=a[n/2,0]+1j*a[n/2+1,0]
        
    if n%2==0 and m%2==0:
        res[m/2,n/2]=a[m/2,n/2]+0j
    return res

def pack(a):
    """Does the inverse of unpack - converts a 2D FFT to HC format."""
    res=numpy.zeros(a.shape,numpy.float32)
    m=a.shape[0]
    n=a.shape[1]
    for i in range(1,(m+1)/2):
        res[i,0]=a[i,0].real
        res[m-i,0]=a[i,0].imag
        if n%2==0:
            res[i,n/2]=a[i,n/2].real
            res[m-i,n/2]=a[i,n/2].imag

        for j in range(1,(n+1)/2):
            if i==1:
                res[0,j]=a[0,j].real
                res[0,n-j]=a[0,j].imag
                if m%2==0:
                    res[m/2,j]=a[m/2,j].real
                    res[m/2,n-j]=a[m/2,j].imag

            res[i,j]=(a[i,j]+a[i,n-j]).real/2
            res[m-i,n-j]=(a[i,n-j]-a[i,j]).real/2
            res[i,n-j]=(a[i,j]-a[i,n-j]).imag/2
            res[m-i,j]=(a[i,j]+a[i,n-j]).imag/2
    res[0,0]=a[0,0].real
    if n%2==0:
        res[0,n/2]=a[0,n/2].real
    else:
        res[0,n/2]=a[0,n/2].real
        res[0,n/2+1]=a[0,n/2].imag
        
    if m%2==0:
        res[m/2,0]=a[m/2,0].real
        if n%2==0:
            res[m/2,n/2]=a[m/2,n/2].real
    else:
        res[m/2,0]=a[m/2,0].real
        res[m/2+1,0]=a[m/2,0].imag
    return res

def correlate(a,b):
    """correlate a with b"""
    #first do it the numpy way...
    m=a.shape[0]
    n=a.shape[1]
    corr=numpy.fft.ifft2(numpy.fft.fft2(a)*numpy.conjugate(numpy.fft.fft2(numpy.fft.fftshift(b)))).real
    #and now do it using half complex format.
    aa2=r2hc(numpy.fft.fftshift(b))
    #do the conjugate.
    mm=aa2.shape[0]/2+1
    nn=aa2.shape[1]/2+1
    aa2[:mm,nn:]*=-1
    aa2[mm:,:nn]*=-1

    aa1=r2hc(a)
    aa3=aa1*aa2
    #n=aa2.shape[0]
    #now multiply the two together
    res=numpy.zeros(a.shape,numpy.float32)
    res[0,0]=aa1[0,0]*aa2[0,0]
    if n%2==0:
        res[0,n/2]=aa1[0,n/2]*aa2[0,n/2]
    if m%2==0:
        res[m/2,0]=aa1[m/2,0]*aa2[m/2,0]
        if n%2==0:
            res[m/2,n/2]=aa1[m/2,n/2]*aa2[m/2,n/2]
    else:
        pass
        
    for i in range(1,(m+1)/2):
        #a10i=aa1[0,i]+1j*aa1[0,n-i]
        #a20i=aa2[0,i]+1j*aa2[0,n-i]
        #a30i=a10i*a20i
        #res[0,i]=a30i.real
        #res[0,n-i]=a30i.imag

        #a1i0=aa1[i,0]+1j*aa1[n-i,0]
        #a2i0=aa2[i,0]+1j*aa2[n-i,0]
        #a3i0=a1i0*a2i0
        #res[i,0]=a3i0.real
        #res[n-i,0]=a3i0.imag
        res[i,0]=aa1[i,0]*aa2[i,0]-aa1[m-i,0]*aa2[m-i,0]
        res[m-i,0]=aa1[i,0]*aa2[m-i,0]+aa1[m-i,0]*aa2[i,0]

        if n%2==0:
            #a1in=aa1[i,n/2]+1j*aa1[n-i,n/2]
            #a2in=aa2[i,n/2]+1j*aa2[n-i,n/2]
            #a3in=a1in*a2in
            #res[i,n/2]=a3in.real
            #res[n-i,n/2]=a3in.imag
            res[i,n/2]=aa1[i,n/2]*aa2[i,n/2]-aa1[m-i,n/2]*aa2[m-i,n/2]
            res[m-i,n/2]=aa1[i,n/2]*aa2[m-i,n/2]+aa1[m-i,n/2]*aa2[i,n/2]


        for j in range(1,(n+1)/2):
            if i==1:
                res[0,j]=aa1[0,j]*aa2[0,j]-aa1[0,n-j]*aa2[0,n-j]
                res[0,n-j]=aa1[0,j]*aa2[0,n-j]+aa1[0,n-j]*aa2[0,j]
                if m%2==0:
                    #a1ni=aa1[n/2,i]+1j*aa1[n/2,n-i]
                    #a2ni=aa2[n/2,i]+1j*aa2[n/2,n-i]
                    #a3ni=a1ni*a2ni
                    #res[n/2,i]=a3ni.real
                    #res[n/2,n-i]=a3ni.imag
                    res[m/2,j]=aa1[m/2,j]*aa2[m/2,j]-aa1[m/2,n-j]*aa2[m/2,n-j]
                    res[m/2,n-j]=aa1[m/2,j]*aa2[m/2,n-j]+aa1[m/2,n-j]*aa2[m/2,j]



            #first unpack hc
            #a1ij=aa1[i,j]-aa1[n-i,n-j]+1j*(aa1[i,n-j]+aa1[n-i,j])
            #a2ij=aa2[i,j]-aa2[n-i,n-j]+1j*(aa2[i,n-j]+aa2[n-i,j])

            #a1inmj=aa1[i,j]+aa1[n-i,n-j]+1j*(aa1[n-i,j]-aa1[i,n-j])
            #a2inmj=aa2[i,j]+aa2[n-i,n-j]+1j*(aa2[n-i,j]-aa2[i,n-j])
            #now multiply these unpacked values...
            #a3ij=a1ij*a2ij
            #a3inmj=a1inmj*a2inmj
            #and now repack the result.
            #res[i,j]=(a3ij+a3inmj).real/2
            #res[n-i,n-j]=(a3inmj-a3ij).real/2

            #res[i,n-j]=(a3ij-a3inmj).imag/2
            #res[n-i,j]=(a3ij+a3inmj).imag/2

            res[i,j]=aa3[i,j]+aa3[m-i,n-j]-aa3[i,n-j]-aa3[m-i,j]
            res[m-i,n-j]=aa1[i,j]*aa2[m-i,n-j]+aa1[m-i,n-j]*aa2[i,j]+aa1[m-i,j]*aa2[i,n-j]+aa1[i,n-j]*aa2[m-i,j]

            res[i,n-j]=aa1[i,j]*aa2[i,n-j]-aa1[m-i,n-j]*aa2[m-i,j]+aa1[i,n-j]*aa2[i,j]-aa1[m-i,j]*aa2[m-i,n-j]
            res[m-i,j]=aa1[i,j]*aa2[m-i,j]-aa1[m-i,n-j]*aa2[i,n-j]+aa1[m-i,j]*aa2[i,j]-aa1[i,n-j]*aa2[m-i,n-j]
            
    #and then hc2r.
    res=hc2r(res)
    #corr=r2hc(corr)
    return corr,res


## if __name__=="__main__":
##     import sys
##     n=16
##     m=16
##     if len(sys.argv)==3:
##         m=int(sys.argv[1])
##         n=int(sys.argv[2])
##     a=numpy.random.random((m,n)).astype("f")
##     a[3:5,4:6]=20.
##     b=numpy.zeros((m,n),"f")
##     b[m/2-1:m/2+1,n/2-1:n/2+1]=1.
##     #first as a test, write a.
##     open("sh.dat","w").write(a.tostring())
##     tmp=r2hc(a)

##     #corr=numpy.conjugate(numpy.fft.fft2(numpy.fft.fftshift(b)))
##     corr=r2hc(numpy.fft.fftshift(b))#how do I do the conjugate?
##     #Do the conjugation...
##     mm=corr.shape[0]/2+1
##     nn=corr.shape[1]/2+1
##     corr[:mm,nn:]*=-1
##     corr[mm:,:nn]*=-1
##     #Note, this is the same as doing:
##     #corr=r2hc(numpy.fft.ifft2(numpy.conjugate(numpy.fft.fft2(numpy.fft.fftshift(b)))).real)
##     open("corr.dat","w").write(corr.tostring())

##     #Now, do the convolution in python (for comparisons...).
##     c=correlate(a,b)

##     os.system("./corrfftw %d %d"%(m,n))

##     txt=open("ffta.dat").read()
##     ffta=numpy.fromstring(txt,numpy.float32)
##     ffta.shape=m,n

##     #This shows that r2hc works.... (ie my python version gives same as fftw)
##     print "checking r2hc",min(numpy.fabs(ffta.ravel())),max(numpy.fabs(ffta.ravel())),min(numpy.fabs(ffta-tmp).ravel()),max(numpy.fabs(ffta-tmp).ravel())

##     txt=open("fftab.dat").read()
##     fftab=numpy.fromstring(txt,numpy.float32)
##     fftab.shape=m,n

##     pyres=numpy.fft.ifft2(numpy.fft.fft2(a)*numpy.conjugate(numpy.fft.fft2(numpy.fft.fftshift(b)))).real

##     pyfftab=r2hc(pyres)
##     print "fftab diff",min(numpy.fabs(fftab.ravel())),max(numpy.fabs(fftab.ravel())),min(numpy.fabs(fftab-pyfftab).ravel()),max(numpy.fabs(fftab-pyfftab).ravel())
    

##     #and check the result...
##     pyres*=m*n

##     txt=open("res.dat").read()
##     res=numpy.fromstring(txt,numpy.float32)
##     res.shape=m,n

##     print "result diff",min(numpy.fabs(res.ravel())),max(numpy.fabs(res.ravel())),min(numpy.fabs(res-pyres).ravel()),max(numpy.fabs(res-pyres).ravel())

##     #Now compare res with pyres
##     #Fiddle round with creating corr until you get them matching...
