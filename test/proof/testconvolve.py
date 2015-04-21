"""To test whether can bin image before convolution."""

import numpy
import pylab
import util.centroid
import util.elong
import scipy.signal
import science.infScrn
import util.tel

def binimg(img,bin):
    img=img.copy()
    print img.shape
    img.shape=img.shape[0],img.shape[1]/bin,bin
    img=img.sum(2)
    img.shape=img.shape[0]/bin,bin,img.shape[1]
    img=img.sum(1)
    return img


def test(n=32,bin=2):
    #img=util.centroid.createAiryDisc(n,4,1.5,.25)
    #img*=100./img.sum()

    phase=science.infScrn.makeScrnQuick(n/2,0.5,l0=3.,r0=0.1,seed=0)
    npup=phase.shape[0]
    nsubx=1
    fftsize=npup/nsubx*2
    phasesize=npup/nsubx
    nimg=fftsize#phasesize
    ncen=fftsize#phasesize
    c=util.centroid.centroid(nsubx,util.tel.Pupil(npup,npup/2,0,nsubx),fftsize=fftsize,binfactor=None,phasesize=phasesize,nimg=nimg,ncen=ncen)#addPoisson=0,sig=1.
    c.easy()
    c.reformatPhs(phase)
    c.runCalc({"cal_source":0})
    #c.outputData is the slopes
    #c.cmodbimg is the wfs image (divided into subaps).
    #Then, put your data into self.reorderedPhs... 
    img=c.reformatImg()# can be used to get a displayable image.




    #img2,c=util.centroid.compute(scrn,1)
    


    #img+=numpy.random.normal(0,1.,(n,n))#1e readout noise.
    psf=util.elong.make(n,1,n/4,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5)[1]
    conv=scipy.signal.convolve2d(img,psf,"same",boundary="wrap")
    convbin=binimg(conv,bin)
    imgbin=binimg(img,bin)
    psfbin=binimg(psf,bin)
    psfbin2=util.elong.make(n/bin,1,n/4/bin,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5/bin,xoff=-(0.5-0.5/bin),yoff=-(0.5-0.5/bin))[1]
    binconv=scipy.signal.convolve2d(imgbin,psfbin,"same",boundary="wrap")
    binconv2=scipy.signal.convolve2d(imgbin,psfbin2,"same",boundary="wrap")
    #compare conv and binconv2... looks good!  Not identical, but similar, and probably good enough!  At least for simple Gaussian spot profiles.
    return img,psf,conv,convbin,imgbin,psfbin,binconv,phase,psfbin2,binconv2


def test2(n=16):
    fftsize=128
    #nimg=16#gets changed later.
    phase=science.infScrn.makeScrnQuick(n,0.5,l0=3.,r0=0.1,seed=0)
    npup=phase.shape[0]
    nsubx=1
    phasesize=n
    nimg=fftsize#phasesize
    ncen=fftsize#phasesize
    c=util.centroid.centroid(nsubx,util.tel.Pupil(npup,npup/2,0,nsubx),fftsize=fftsize,binfactor=None,phasesize=phasesize,nimg=nimg,ncen=ncen)#addPoisson=0,sig=1.
    nimg=16
    ncen=16
    bin=fftsize/nimg
    c.easy()
    c.reformatPhs(phase)
    c.runCalc({"cal_source":0})
    #c.outputData is the slopes
    #c.cmodbimg is the wfs image (divided into subaps).
    #Then, put your data into self.reorderedPhs... 
    img=c.reformatImg()# can be used to get a displayable image.
    imgpad=numpy.zeros((fftsize*2,fftsize*2),"f")
    imgpad[fftsize/2:fftsize*3/2,fftsize/2:fftsize*3/2]=img

    #psf=util.elong.make(fftsize*2,1,fftsize/2,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5)[1]
    psf=util.elong.make(spotsize=fftsize*2,nsubx=1,wfs_n=n,wfs_nfft=fftsize,wfs_nimg=nimg,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5)[1]
    #first, bin after convolution.
    conv=scipy.signal.convolve2d(imgpad,psf,"same",boundary="wrap")
    conv=conv[fftsize/2:fftsize*3/2,fftsize/2:fftsize*3/2]
    convbin=binimg(conv,bin)

    #Now bin before convolution.
    imgbin=binimg(imgpad,bin)#bin to 16, then pad to 32...
    #psfbin=util.elong.make(n*2,1,n/2,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5/8,xoff=-(0.5-0.5/8),yoff=-(0.5-0.5/8))[1]
    psfbin=util.elong.make(fftsize/bin*2,1,wfs_n=n,wfs_nfft=fftsize/bin,wfs_nimg=nimg,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5,xoff=-(0.5-0.5/bin),yoff=-(0.5-0.5/bin))[1]
    convpostbin=scipy.signal.convolve2d(imgbin,psfbin,"same",boundary="wrap")[nimg/2:nimg*3/2,nimg/2:nimg*3/2]

    #And how about bin, then conv, then bin.
    #bin=2
    imgbin4=binimg(imgpad,bin/2)#bin to 32 then pad to 64.
    #psfbin4=util.elong.make(n*2,1,n/2,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5/bin,xoff=-(0.5-0.5/bin),yoff=-(0.5-0.5/bin))[1]
    psfbin4=util.elong.make(fftsize/(bin/2)*2,1,wfs_n=n,wfs_nfft=fftsize/(bin/2),wfs_nimg=nimg,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5/1,xoff=-(0.5-0.5/(bin/2)),yoff=-(0.5-0.5/(bin/2)))[1]
    convpostbin4=scipy.signal.convolve2d(imgbin4,psfbin4,"same",boundary="wrap")
    convbinbin=binimg(convpostbin4,2)[nimg/2:nimg*3/2,nimg/2:nimg*3/2]

    #and convolve without padding.  No - need the padding to prevent wrapping, so this isn't great.
    imgbin4b=binimg(img,bin/2)#bin to 32.
    psfbin32=util.elong.make(fftsize/(bin/2),1,wfs_n=n,wfs_nfft=fftsize/(bin/2),wfs_nimg=nimg,beacon_depth=5000.,launchDist=20.,launchTheta=30.,telDiam=0.5/1,xoff=-(0.5-0.5/(bin/2)),yoff=-(0.5-0.5/(bin/2)))[1]
    convpostbin4b=scipy.signal.convolve2d(imgbin4b,psfbin32,"same",boundary="wrap")
    convbinbinb=binimg(convpostbin4b,2)

    #Doing binning before convolution and then after seems to give best match.
    return img,psf,conv,convbin,imgbin,psfbin,convpostbin,imgbin4,psfbin4,convbinbin,imgbin4b,convpostbin4b,convbinbinb,psfbin32
    



def testPreBinning(n=16,prebinfactor=4):
    """Tests preBinningFactor in the centroid module"""
    telDiam=4.2
    fftsize=128
    nsubx=7
    npup=n*nsubx
    nimg=16
    ncen=16
    phase=science.infScrn.makeScrnQuick(npup,telDiam,l0=3.,r0=0.15,seed=0)
    phasesize=n
    postbinfactor=8/prebinfactor
    import util.tel
    import util.elong
    lgs=util.elong.make(fftsize/prebinfactor*2,nsubx,wfs_n=n,wfs_nfft=fftsize/prebinfactor,wfs_nimg=nimg,clipsize=fftsize/prebinfactor,telDiam=telDiam,beacon_alt=11000.,beacon_depth=1000.,xoff=(0.5-0.5/prebinfactor),yoff=(0.5-0.5/prebinfactor))

    c=util.centroid.centroid(nsubx,util.tel.Pupil(npup,npup/2,0),fftsize=fftsize,binfactor=None,phasesize=phasesize,nimg=nimg,ncen=ncen,preBinningFactor=prebinfactor,spotpsf=lgs[0],clipsize=fftsize/prebinfactor)
    c.easy()
    c.reformatPhs(phase)
    c.runCalc({"cal_source":0})
    img=c.reformatImg()# can be used to get a displayable image.


    #and now test without prebinning.
    lgs2=util.elong.make(256,nsubx,wfs_n=n,wfs_nfft=fftsize,wfs_nimg=nimg,clipsize=fftsize,telDiam=telDiam,beacon_alt=11000.,beacon_depth=1000.)
    c2=util.centroid.centroid(nsubx,util.tel.Pupil(npup,npup/2,0),fftsize=fftsize,binfactor=None,phasesize=phasesize,nimg=nimg,ncen=ncen,preBinningFactor=1,spotpsf=lgs2[0],clipsize=128)

    c2.easy()
    c2.reformatPhs(phase)
    c2.runCalc({"cal_source":0})
    img2=c2.reformatImg()# can be used to get a displayable image.


    return img,img2

if __name__=="__main__":
    import sys
    n=16
    prebinfactor=4
    if len(sys.argv)>1:
        n=int(sys.argv[1])
    if len(sys.argv)>2:
        prebinfactor=int(sys.argv[2])
    img,img2=testPreBinning(n,prebinfactor)
    import pylab
    pylab.subplot(1,3,1)
    pylab.imshow(img,interpolation="nearest",cmap="gray")
    pylab.subplot(1,3,2)
    pylab.imshow(img2,interpolation="nearest",cmap="gray")
    pylab.subplot(1,3,3)
    pylab.imshow(img-img2,interpolation="nearest",cmap="gray")
    print numpy.max(img),numpy.max(img-img2),numpy.min(img-img2)
    pylab.show()
