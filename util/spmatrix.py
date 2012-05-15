"""Additional sparse matrix operations, on top of those from scipy.sparse"""
import numpy as na
import scipy.sparse
def addToDiag(arr,val):
    """val can be a vector or a single value."""
    if arr.format!="csc":
        raise Exception("spmatrix.addToDiag only usable for scipy.sparse csc type matricees")
    v=scipy.sparse.lil_matrix(arr.shape)
    if type(val) in [type(0.),type(0)]:
        v.setdiag(na.ones(min(arr.shape))*val)
        #v=na.zeros((min(arr.shape),),na.float32)
        #v[:]=val
    else:
        if v.shape[0]!=min(arr.shape):
            raise Exception("spmatrix.addToDiag wrong sizes %s %s"%(str(arr.shape),str(val.shape)))
        v.setdiag(val)
    arr=arr+v
    return arr
    #for i in xrange(arr.indptr.shape[0]-1):
    #    if arr.indptr[i]==arr.indptr[i+1]:
    #        print arr.indptr
    #        raise Exception("spmatrix.addToDiag - the %dth entry doesn't have a diagonal element"%i)
    #    pos=arr.rowind[arr.indptr[i]:arr.indptr[i+1]].searchsorted(i)+arr.indptr[i]
    #    print pos,arr.rowind.shape,i
    #    if arr.rowind[pos]==i:#this is the diag
    #        arr.data[pos]+=v[i]
    #    else:
            

def getDiagonal(arr):
    if arr.format!="csc":
        raise Exception("spmatrix.getDiagonal only usable for scipy.sparse csc type matricees")
    diag=na.zeros((min(arr.shape),),arr.dtype)
    for i in xrange(diag.shape[0]):
        pos=arr.rowind[arr.indptr[i]:arr.indptr[i+1]].searchsorted(i)+arr.indptr[i]
        if arr.rowind[pos]==i:#this is the diag
            diag[i]=arr.data[pos]
    return diag
def strip(mx,val):
    """Sets to zero all elements with abs()<=val."""
    if mx.format!="csc":
        raise Exception("spmatrix.getDiagonal only usable for scipy.sparse csc type matricees")
    imax=mx.data.shape[0]
    i=0
    nrem=0
    while i<imax:
        if abs(mx.data[i])<=val:
            mx.data[i:-1]=mx.data[i+1:,].copy()
            mx.data=mx.data[:-1]
            mx.rowind[i:-1]=mx.rowind[i+1:,].copy()
            mx.rowind=mx.rowind[:-1]
            imax-=1
            col = na.searchsorted(mx.indptr, i+1)-1
            mx.indptr[col+1:,]-=1
            nrem+=1
        else:
            i+=1
    mx.rowind=mx.rowind.copy()#this should allow the old version memory to be released
    mx.data=mx.data.copy()
    mx._check()
        
def dotWithSelfTransposed(mx,tol=0.,fmt="csc",fudge=1):
    """
    cmod.svd.csrDotCsc is much faster, and multithreaded...
    Dot mx (a scipy.sparse matrix) with mx transposed.
    use tol to specify the value below which data is ignored.
    fmt is the returned format.  csr is faster than csc, but means the conversion needs to take place later on.
    If fudge is set, will make sure that every row/col has an entry, to avoid the result being singular (necessary but not sufficient), so that it can probably be solved...
    """
    print "util.spmatrix.dotWithSelfTransposed - can take a while (up to 20 minutes for an EAGLE matrix)"
    if mx.format=="csr":#this is probably the faster method...
        lil=scipy.sparse.lil_matrix((mx.shape[0],mx.shape[0]))
        for i in xrange(mx.shape[0]):
            r1=mx.colind[mx.indptr[i]:mx.indptr[i+1]]
            d1=mx.data[mx.indptr[i]:mx.indptr[i+1]]
            n1=mx.indptr[i+1]-mx.indptr[i]
            for j in xrange(i,mx.shape[0]):
                r2=mx.colind[mx.indptr[j]:mx.indptr[j+1]]
                d2=mx.data[mx.indptr[j]:mx.indptr[j+1]]
                n2=mx.indptr[j+1]-mx.indptr[j]
                val=0.
                cnt1=0
                cnt2=0
                while cnt1<n1 and cnt2<n2:#for each element of the first matrix...
                    while cnt2<n2 and r2[cnt2]<r1[cnt1]:
                        cnt2+=1
                    if cnt2<n2 and r2[cnt2]==r1[cnt1]:
                        val+=d1[cnt1]*d2[cnt2]
                    cnt1+=1
                if na.fabs(val)>tol:
                    lil[i,j]=val
                    lil[j,i]=val
    elif mx.format=="csc":#see if this is faster?
        lil=scipy.sparse.lil_matrix((mx.shape[0],mx.shape[0]))
        import time
        t1=time.time()
        for k in xrange(mx.shape[1]):
            col=mx.rowind[mx.indptr[k]:mx.indptr[k+1]]
            data=mx.data[mx.indptr[k]:mx.indptr[k+1]]
            n=mx.indptr[k+1]-mx.indptr[k]
            #print "k %d time %g"%(k,time.time()-t1)
            for l in xrange(n):#for each entry in the row...
                for m in xrange(l,n):#for the entries past this...
                    c1=int(col[l])
                    c2=int(col[m])
                    val=data[l]*data[m]
                    lil[c1,c2]+=val
                    if l!=m:
                        lil[c2,c1]+=val
                            

    elif mx.format=="cscslow":#this is fairly slow - may be faster to swap to a csr type (mx.tocsr)
        lil=scipy.sparse.lil_matrix((mx.shape[0],mx.shape[0]))
        for i in xrange(mx.shape[0]):
            for j in xrange(i,mx.shape[0]):
                val=0.
                for k in xrange(mx.shape[1]):
                    col=mx.rowind[mx.indptr[k]:mx.indptr[k+1]]
                    data=mx.data[mx.indptr[k]:mx.indptr[k+1]]
                    n=mx.indptr[k+1]-mx.indptr[k]
                    pos=col.searchsorted(i)
                    if pos<n and col[pos]==i:#has an entry...
                        pos2=pos+col[pos:,].searchsorted(j)
                        if pos2<n and col[pos2]==j:#also has an entry...
                            val+=data[pos]*data[pos2]
                if na.fabs(val)>tol:
                    lil[i,j]=val
                    lil[j,i]=val
    else:
        raise Exception("spmatrix.dotWithSelfTransposed only usable for scipy.spares csc type matricees")
    if fudge:#now make sure that every row/col has something.
        for ind,row in enumerate(lil.data):
            if len(row)==0:#nothing in this row... make it same as prev or next diag element.
                p=n=0.
                if ind>0:
                    p=lil[ind-1,ind-1]
                if ind<mx.shape[0]-1:
                    n=lil[ind+1,ind+1]
                m=min(p,n)
                if m==0:
                    m=max(p,n)
                if m==0:
                    m=0.1
                lil[ind,ind]=m
                print "util.spmatrix.dotWithSelfTransposed fudging row %d with val %g"%(ind,m)
    if fmt=="csc":
        print "util.spmatrix.dotWithSelfTransposed: converting to csc format (could take up to 10 minutes)"
        mat=lil.tocsc()#liltocsc(lil)
    elif fmt=="csr":#this is fairly fast (2s)
        mat=lil.tocsr()
    else:
        mat=lil
    return mat

def liltocsc(mx):
    """Convert a lil matrix to a csc matrix"""
    nzmax=mx.getnnz()
    data=na.zeros(nzmax,dtype=mx.dtype)
    rowind=na.zeros(nzmax,dtype=na.int32)
    indptr=na.empty(mx.shape[0]+1,dtype=na.int32)
    indptr[:,]=nzmax
    pos=0
    for i in xrange(mx.shape[1]):
        indptr[i]=pos
        for j,col in enumerate(mx.rows):
            if i in col:
                rowind[pos]=j
                data[pos]=mx.data[j][col.index(i)]
                pos+=1
    indptr[-1]=pos
    return scipy.sparse.csc_matrix((data,rowind,indptr),dims=mx.shape,nzmax=nzmax)


###This is copied from scipy.linsolve.  The idea is to make it better for this situation, including the ability to factorise the matrix before doing the computation...

import scipy.linsolve

def spsolve(A, b, permc_spec=2):
    """Taken from scipy.linsolve, and updated to raise an exception if error occurs."""
    if not hasattr(A, 'tocsr') and not hasattr(A, 'tocsc'):
        raise ValueError, "sparse matrix must be able to return CSC format--"\
              "A.tocsc()--or CSR format--A.tocsr()"
    if not hasattr(A, 'shape'):
        raise ValueError, "sparse matrix must be able to return shape" \
                " (rows, cols) = A.shape"
    M, N = A.shape
    if (M != N):
        raise ValueError, "matrix must be square"


    mat, csc = scipy.linsolve.linsolve._toCS_superLU( A )
    if csc:
        index0 = mat.rowind
    else:
        index0 = mat.colind
    ftype, lastel, data, index1 = mat.ftype, mat.nnz, mat.data, mat.indptr
    gssv = eval('scipy.linsolve.linsolve._superlu.' + ftype + 'gssv')
    #print "data-ftype: %s compared to data %s" % (ftype, data.dtype.char)
    #print "Calling _superlu.%sgssv" % ftype
    data,info=gssv(N, lastel, data, index0, index1, b, csc, permc_spec)
    if info!=0:
        raise Exception("util.spmatrix.spsolve failed with info %d"%info)
    return data

def cscselect(csc,yf,yt,xf,xt):
    """select part of a csc matrix."""
    print "I think this might be wrong - and I think the csr version is right... will have to think about this sometime"
    arr=na.zeros((yt-yf,xt-xf),csc.dtype.char)
    for i in range(yf,yt):
	rf=csc.indptr[i]
        rt=csc.indptr[i+1]
        for j in range(rf,rt):
            if csc.rowind[j]>=xf and csc.rowind[j]<xt:
                arr[i-yf,csc.rowind[j]-xf]=csc.data[j]
    return arr
def csrselect(csr,yf,yt,xf,xt):
    """select part of a csr matrix."""
    arr=na.zeros((yt-yf,xt-xf),csr.dtype.char)
    for i in range(yf,yt):
        f=csr.indptr[i]
        t=csr.indptr[i+1]
        for j in range(f,t):
            if csr.colind[j]>=xf and csr.colind[j]<xt:
                arr[i-yf,csr.colind[j]-xf]=csr.data[j]
    return arr
