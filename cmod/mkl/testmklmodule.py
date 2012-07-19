import cmod.mkl
import numpy
import sys,time
#import util.dot as quick

def runtest(size,usesdd=1,typ="d",tst=1,overwriteU=0):
    aa=None
    if (size<=32768 and typ=="f") or (size<=16384 and typ=="d"):
        a=numpy.random.random((size,size)).astype(typ)
        vt=numpy.zeros((size,size),typ)
        if overwriteU:
            u=a
        else:
            u=numpy.zeros((size,size),typ)
            
    else:
        #create a and vt memory mapped.
        a=cmod.utils.mmapArray("a.dat",(size,size),typ)
        for i in range(size):
            a[i]=numpy.random.random((size,)).astype(typ)
        vt=cmod.utils.mmapArray("vt.dat",(size,size),typ)
        if overwriteU:
            u=a
        else:
            #create u memory mapped.
            u=cmod.utils.mmapArray("vt.dat",(size,size),typ)


    evals=numpy.zeros((size,),typ)
    minimum=maximum=0
    t1=t2=inv1=inv2=0
    nsvd=None
    if tst:
        aa=a.copy()
        t1=time.time()
        nsvd=numpy.linalg.svd(a)
        t2=time.time()
        inv1=makeInv(nsvd[0],nsvd[1],nsvd[2],0.1)
    wsize=cmod.mkl.svd(a,u,evals,vt,None,usesdd)
    if (wsize<=32768**2*4+32768*7 and typ=="f") or (wsize<=16384**2*4+16384*7 and typ=="d"):
        work=numpy.zeros((wsize,),typ)
    else:
        work=cmod.utils.mmapArray("work.dat",(wsize,),typ)
    t3=time.time()
    cmod.mkl.svd(a,u,evals,vt,work,usesdd)
    t4=time.time()
    if tst:
        inv2=makeInv(vt,evals,u,0.1)
        d=(inv1-inv2).ravel()
        minimum=min(d)
        maximum=max(d)
    return nsvd,aa,u,evals,vt,inv1,inv2,t2-t1,t4-t3,minimum,maximum,wsize

def makeInv(u,evals,vt,minEig=0.):
    eval2=numpy.where(evals<minEig,0,1/evals)
    ut=u.transpose()
    v=vt.transpose()
    for i in xrange(min(ut.shape[0],eval2.shape[0])):
        ut[i]*=eval2[i]
    if v.shape[0]>ut.shape[0]:#ncents>nmodes...
        inv=quick.dot(v[:,:ut.shape[0]], ut)
    else:
        inv=quick.dot(v,ut[:vt.shape[1]])
    return inv

def loadDat(size=16384,dtype="f",uf="u.dat",vf="v.dat",ef="evals.dat"):
    u=numpy.zeros((size,size),dtype)
    v=numpy.zeros((size,size),dtype)
    e=numpy.zeros((size,),dtype)
    wsize=8
    if dtype=="f":wsize=4
    f=open(uf)
    for i in range(size):
        data=f.read(size*wsize)
        u[i]=numpy.fromstring(data,dtype)
    f.close()
    f=open(vf)
    for i in range(size):
        data=f.read(size*wsize)
        v[i]=numpy.fromstring(data,dtype)
    f.close()
    f=open(ef)
    data=f.read(size*wsize)
    e[:]=numpy.fromstring(data,dtype)
    f.close()
    return u,e,v

        
if __name__=="__main__":
    size=int(sys.argv[1])
    usesdd=int(sys.argv[2])
    typ=sys.argv[3]
    tst=int(sys.argv[4])
    overwriteU=int(sys.argv[5])
    save=int(sys.argv[6])
    res=runtest(size,usesdd,typ,tst,overwriteU)
    txt= "size %d sdd %d typ %s overwriteU %d Time %g (numpy %g) min/max diff: %g %g\n"%(size,usesdd,typ,overwriteU,res[8],res[7],res[9],res[10])
    print txt
    open("timing.txt","a").write(txt)
    if save:
        u=res[2]
        evals=res[3]
        v=res[4]
        f=open("u.dat","w")
        for i in range(size):
            f.write(u[i].tostring())
        f.close()
        f=open("v.dat","w")
        for i in range(size):
            f.write(v[i].tostring())
        f.close()
        f=open("evals.dat","w")
        f.write(evals.tostring())
        f.close()
