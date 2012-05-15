import svd,gist
import util.FITS,Numeric,scipy.sparse,numpy,cmod.utils
Adata=Numeric.array(util.FITS.Read("Adata.fits")[1].astype("f"))
Arowind=Numeric.array(util.FITS.Read("Arowind.fits")[1].astype("i"))
Aindptr=Numeric.array(util.FITS.Read("Aindptr.fits")[1].astype("i"))
m=578
n=1536
W=Numeric.zeros((n,),"f")
Ut=Numeric.zeros((m,m),"f")
Vt2=Numeric.zeros((m,n),"f")
Vt=cmod.utils.mmapArray("tmp.mmap",(m,n),"f")#this can be None...

rmx=Numeric.zeros((n,m),"f")#transposed version.
rmxData=Numeric.zeros((90000,),"f")
rmxRowind=Numeric.zeros((90000,),"i")
rmxIndptr=Numeric.zeros((n+1,),"i")
rmx=svd.sparseMatrixCreate(90000,m,n,rmxData,rmxRowind,rmxIndptr)
storefloat=1
sfloat=1
minGIVal=1.
svd.dolibsvd(m,n,Adata,Arowind,Aindptr,Ut,W,Vt,0,0,rmx,0.01,0.,storefloat,sfloat,minGIVal)
gist.window(0);gist.palette("gray.gp");gist.fma();gist.plg(W)
gist.window(1);gist.palette("gray.gp");gist.fma();gist.pli(Ut)
gist.window(2);gist.palette("gray.gp");gist.fma();gist.pli(Vt)
gist.window(3);gist.palette("gray.gp");gist.fma();gist.pli(rmx)


def makeReconmx(ut,w,vt):
    """If used with outputs from other algorithms, u,v must be transposed correctly...
    If calling this from results of dolibsvd, set to zero entries of W smaller than e.g. 0.01...,
    or even better, call with w[:xxx] where xxx is the position beyond which you want to ignore - this is fastest...
    """
    #first dot a matrix with w on the diag with ut.
    ut=ut.copy()#so we don't destroy ut.
    valid=min(w.shape[0],ut.shape[0])
    for i in range(valid):
        if w[i]==0.:
            ut[i]=0.
        else:
            ut[i]*=1/w[i]
    v=Numeric.transpose(vt)
    print v.shape,ut.shape,valid
    #rmx=Numeric.matrixmultiply(v[:,:ut.shape[0]],ut)
    rmx=Numeric.matrixmultiply(v[:,:valid],ut[:valid])
    return rmx


def rotateArray(a,size,x):
    if x == 0:
        return;
    j = start = 0;
    t1 = a[0];
    for i in xrange(size):
        if j>=x:
            n=j-x
        else:
            n=j+size-x
        t2 = a[n];
        a[n] = t1;
        t1 = t2;
        j = n;
        if j == start:
            j+=1
            start = j;
            t1 = a[j];
    
def getRowPositionQuick(s,e,row,rowind):
    """s,e are start and end locations, rowind is the ordered array and row is the position to find in the array."""
    eo=e
    so=s
    dm=(s+e)/2.
    m=int(dm+0.5)
    while (m>so and m<eo and not(rowind[m-1]<row and rowind[m]>=row)):
        if(rowind[m]<row):
            s=m
            dm=(dm+e)/2.
        elif(rowind[m]>row):
            e=m
            dm=(s+dm)/2.
        m=int(dm+0.5)
        print dm,m
    return m

import scipy.sparse
lil=scipy.sparse.lil_matrix((m,n))
for i in range(Adata.shape[0]):
    lil[i/n,i%n]=1

csc=lil.tocsc()
csc.data=Adata
csc.indptr=Aindptr
csc.rowind=Arowind
pmx=csc.todense()

csc=scipy.sparse.csc_matrix((numpy.array(rmxData),numpy.array(rmxRowind),numpy.array(rmxIndptr)),(m,n))
gist.window(3);gist.palette("gray.gp");gist.fma();gist.pli(Numeric.transpose(csc.todense()))
