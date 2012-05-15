import cmod.sor
import numpy

def fit(centx,centy,output,sorstruct):#mask,avpistDivisor=None,conv=0.01,maxiters=1000):
    """centx,y contain the centroids, mask is the pupil mask, and avpistDivisor is optional"""
    #if avpistDivisor==None:
    #    avpistDivisor=Numeric.sum((mask>0).flat)
    cmod.sor.fit(centx,centy,output,sorstruct)#mask.flat,avpistDivisor,conv,maxiters)

def createSorMask(pupil,wfs_n=1,minarea=0.):
    """typically, pupil will be the telescope pupil (phase points), wfs_n will be number of pixels (of phase) per subaperture in 1 direction, and minarea will be the fraction of pixels required to have a valid subaperture.  If n=1 and minarea=0 (default), then the pupil function will be exactly that given...
    This should be flattened (mask.flat) before being used with util.sor.
    """
    n=wfs_n
    if type(pupil)==numpy.ndarray:#Numeric.ArrayType:
        pass
    else:
        pupil=pupil.fn
    npup=pupil.shape[0]
    nx=npup/n
    pupfn=numpy.zeros((nx+2,nx+2),numpy.int32)
    imirr=numpy.zeros((nx,nx),numpy.int32)
    total=0
    for i in range(nx):
        for j in range(nx):
            subarea=numpy.sum(numpy.sum(pupil[i*n:(i+1)*n,j*n:(j+1)*n]))
            if subarea>=minarea*n*n:
                pupfn[i+1,j+1]=1
                total+=1

    for i in range(nx):
        for j in range(nx):
            k=i+1
            l=j+1
            fl2=numpy.array(pupfn[k-1:k+2,l-1:l+2],copy=1)
            fl2[0,0]=1
            fl2[0,2]=1
            fl2[2,0]=1
            fl2[2,2]=1
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,1,1],[1,1,1],[1,1,1]],numpy.int32) ) ) )) :
                imirr[j][i]=1
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,0,1],[0,1,1],[1,1,1]],numpy.int32) ) ) )) :
                imirr[j][i]=2
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,0,1],[1,1,0],[1,1,1]],numpy.int32) ) ) )) :
                imirr[j][i]=3
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,1,1],[0,1,1],[1,0,1]],numpy.int32) ) ) )) :
                imirr[j][i]=4
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,1,1],[1,1,0],[1,0,1]],numpy.int32) ) ) )) :
                imirr[j][i]=5
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,0,1],[1,1,1],[1,1,1]],numpy.int32) ) ) )) :
                imirr[j][i]=6
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,1,1],[1,1,1],[1,0,1]],numpy.int32) ) ) )) :
                imirr[j][i]=7
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,1,1],[0,1,1],[1,1,1]],numpy.int32) ) ) )) :
                imirr[j][i]=8
            if( numpy.alltrue(numpy.alltrue( numpy.equal( fl2 , numpy.array([[1,1,1],[1,1,0],[1,1,1]],numpy.int32) ) ) )) :
                imirr[j][i]=9
    return imirr

def reconstruct(pupil,centx,centy,minarea=0.1):
    if len(centx.shape)!=2:
        raise Exception("centx/y should be square")
    wfsn=pupil.shape[0]/centx.shape[0]
    output=numpy.zeros(centx.shape,centx.typecode())
    sormask=createSorMask(pupil,wfsn,minarea)
    sorstr=cmod.sor.init(sormask.flat,numpy.sum((sormask>0).flat))
    iter,err=cmod.sor.fit(centx,centy,output,sorstr)
    return output,iter,err
