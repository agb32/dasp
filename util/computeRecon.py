"""Code to make a reconstructor from a very big (sparse) poke matrix.
See the __main__ part at the bottom for the recommended way to start...

After doing this and finding a typical minimum value to use for the reconstructor, if you need to investigate different svd conditioning, can then use finalDot(rmxtype=1), to give a csr output without needing to make the huge dense reconstructor.  This should be quite a bit faster since paging won't be needed.

Saving files - we by default save them not byteswapped - to they would look like nonsense if read by a normal viewer.  However, util.FITS.Read and associated methods knows about this, and to if it finds an UNORDERD key, it doesn't byteswap even if asked to.
"""
import util.FITS,numpy
import sys,time,os
#import util.dot as quick
#sys.path.insert(0,"/home/ali/c/lapack")
#import svd
acml=0
atlas=0
try:
    import cmod.mkl as mkl
except:
    print "Cannot import mkl - trying acml..."
    try:
        import cmod.acml as mkl
        acml=1
    except:
        print "Cannot import mkl or acml - trying atlas..."
        try:
            import cmod.atlas as mkl
            atlas=1
        except:
            print "Warning: in computeRecon.py - cmod.mkl, cmod.acml and cmod.atlas not found - this may cause problems depending what functions you want - but continuing anyway"

        
import cmod.svd,scipy.sparse

class sparseHolder:
    def __init__(self):
        self.data=None
        self.shape=None
        self.colind=None
        self.indptr=None

class makeRecon:
    def __init__(self,pmxfname,rcond=0.1,regularisation=0.):
        self.pmxfname=pmxfname
        if regularisation!=0:
            pmxfname=pmxfname[:-5]+"reg%g.fits"%regularisation
        self.rcond=rcond
        self.dottedname=pmxfname[:-5]+"_dotted.fits"
        self.uname=pmxfname[:-5]+"_u.fits"
        self.vtname=pmxfname[:-5]+"_vt.fits"
        self.ename=pmxfname[:-5]+"_evals.fits"
        self.invname=pmxfname[:-5]+"_inv%g.fits"%rcond
        self.rmxdensename=pmxfname[:-5]+"_rmxden%g.fits"%rcond
        self.rmxcsrname=pmxfname[:-5]+"_rmxcsr%g.fits"%rcond
        self.rmxcscname=pmxfname[:-5]+"_rmxcsc%g.fits"%rcond
        self.timename=pmxfname[:-5]+"_timing.txt"
        self.cntname=pmxfname[:-5]+"_cnt.txt"
        self.regularisation=regularisation
    def log(self,txt):
        f=open(self.timename,"a")
        f.write(txt)
        f.close()
    def dotdot(self,csc=None):
        """Do the sparse dot product with self transposed..."""
        if csc==None:
            pmxfname=self.pmxfname
            csc=util.FITS.loadSparse(pmxfname,doByteSwap=1)
        if isinstance(csc,numpy.ndarray):
            shape=csc.shape
            if atlas:
                res=numpy.zeros((shape[0],shape[0]),csc.dtype,order="C")
            else:
                res=numpy.zeros((shape[0],shape[0]),csc.dtype,order="F")
            lines=open("/proc/meminfo").readlines()
            for line in lines:
                if "MemTotal:" in line:
                    mem=int(line.split()[1])
                    multiplier={"kB":1024,"b":1,"B":1,"mB":1024*1024,"MB":1024**2,"GB":1024**3,"gB":1024**3}
                    if line.split()[2] in multiplier.keys():
                        mem*=multiplier[line.split()[2]]
                    else:
                        print "WARNING - multiplier %s not known for memory"
                    print "Total system memory %d bytes"%mem
                    break

            #compute the number of blocks to divide the multiplication into...
            nblock=1
            while (shape[0]/nblock*shape[1]*2+(shape[0]/nblock)**2)*csc.itemsize>mem:
                nblock+=1
            print "Doing gemm in %d quadrants, shape %s"%(nblock*nblock,str(csc.shape))
            size=shape[0]/nblock
            t0=time.time()
            for i in range(nblock):
                vstart=i*size
                vend=(i+1)*size
                if i==nblock-1:
                    vend=shape[0]
                a=csc[vstart:vend]#c contig still.
                for j in range(nblock):
                    ustart=j*size
                    uend=(j+1)*size
                    if j==nblock-1:
                        uend=shape[0]
                    b=csc[ustart:uend].transpose()#f contig
                    print "Starting gemm %d %d"%(i,j)
                    t1=time.time()
                    mkl.gemm(a,b,res[vstart:vend,ustart:uend])
                    t2=time.time()
                    print "GEMM time %gs"%(t2-t1)
                    del(b)
                del(a)
            t2=time.time()
            print "Total GEMM time %gs"%(t2-t0)
            open(self.timename,"a").write("mx dot mx dense %gs\n"%(t2-t1))
            extraHeader=None
            if self.regularisation!=0.:
                print "Doing regularisation %g"%self.regularisation
                for i in range(res.shape[0]):
                    res[i,i]+=self.regularisation
                print "Regularisation done"
                extraHeader="REGULARD= %g"%self.regularisation
            util.FITS.Write(res,self.dottedname,doByteSwap=0,extraHeader=extraHeader)
            del(csc)
        else:
            csr=csc.tocsr()
            csc=csr.transpose()
            #resindptr=numpy.zeros((csr.shape[0]+1,),numpy.int32)
            #rescolind=numpy.zeros((csr.nnz*2,),numpy.int32)
            #resdata=numpy.zeros((csr.nnz*2,),numpy.float32)
            print "Doing sparse dot %s %s"%(str(csr.shape),str(csc.shape))
            t1=time.time()
            res=cmod.svd.csrDotCsc(csr.data,csr.colind,csr.indptr,csc.data,csc.rowind,csc.indptr,8)#resdata,rescolind,resindptr,8)
            t2=time.time()
            print "Time for mx dot mx %gs"%(t2-t1)
            open(self.timename,"a").write("mx dot mx %gs\n"%(t2-t1))
            res=scipy.sparse.csr_matrix(res,dims=(csr.shape[0],csr.shape[0]))
            util.FITS.saveSparse(res,self.dottedname,doByteSwap=0)
            del(csc)
            del(csr)
        return res

    def dosvd(self,issparse=1,dtype=numpy.float32,a=None,usesdd=0):
        """Here, a is the densified result from dotdot, ie a square matrix"""
        fname=self.dottedname
        #res=util.FITS.loadSparse(fname,doByteSwap=1)#scipy.sparse.csr_matrix() fails for large sparse matrix...
        if a==None:
            if issparse:
                res=util.FITS.Read(fname,doByteSwap=1,savespace=1)
                shape=eval(res[0]["parsed"]["SHAPE"])
                if res[0]["parsed"].has_key("REGULARD"):
                    print "Using regularised matrix (%g added to diag)"%float(res[0]["parsed"]["REGULARD"])
                if res[0]["parsed"]["FORMAT"]!="csr":
                    raise Exception("Must be CSR matrix")
                #print "Loaded sparse",res.shape,res.data.shape,res.data.dtype.char,type(res),res.indptr.shape
                a=numpy.zeros(shape,res[1].dtype)
                sp=sparseHolder()
                sp.shape=shape
                sp.data=res[1]
                sp.colind=res[3].view(numpy.uint32)
                sp.indptr=res[5].view(numpy.uint32)
                cmod.svd.densifyCsr(a,sp)
                #a=res.todense()#seems to fail if a==18GB.
                print "done todense"
                del(sp)
                del(res)
            else:
                res=util.FITS.Read(fname,doByteSwap=1,savespace=1)
                a=res[1]
                if res[0]["parsed"].has_key("REGULARD"):
                    print "Using regularised matrix (%g added to diag)"%float(res[0]["parsed"]["REGULARD"])
        if not a.flags.c_contiguous:
            print "Transposing a before svd to make c contiguous"
            a=a.transpose()#should now be c-contiguous.
        if not a.flags.c_contiguous:
            raise Exception("A not contiguous")
        if a.dtype!=dtype:
            a=a.astype(dtype)
        print "SVD of matrix with shape",a.shape
        u=a
        #usesdd=0
        vt=numpy.zeros(a.shape,dtype)
        evals=numpy.zeros((a.shape[0],),dtype)
        print "Computing lwork"
        lwork=mkl.svd(a,u,evals,vt,None,usesdd)
        # there is an issue here, that when lwork is large, it isn't accurately represented by a 32 bit float, so is wrong in the conversion.  So, need to work out a larger value that it could be.  This is a bug with intel MKL
        i=0
        aa=numpy.zeros((2,),dtype)
        aa[:]=lwork
        while aa.astype("l")[0]==aa.astype("l")[1]:
            aa[1]=aa[0]+i
            i+=1
        lwork=aa.astype("l")[1]
        print "Got lwork (both) of %d"%lwork
        work=numpy.zeros((lwork,),dtype)
        print "Doing SVD"
        t1=time.time()
        mkl.svd(a,u,evals,vt,work,usesdd)#may take a bit of time of paging at the start, but then should run okay.
        t2=time.time()
        print "SVD time %g"%(t2-t1)
        open(self.timename,"a").write("SVD %gs\n"%(t2-t1))
        # note, u and vt are swapped around compared with what I'm used to from numpy.svd, so I swap them back here...
        util.FITS.Write(evals,self.ename,doByteSwap=0)
        util.FITS.Write(u,self.vtname,doByteSwap=0)
        util.FITS.Write(vt,self.uname,doByteSwap=0)
        return vt,evals,u

    def makeInv(self):
        """Use the SVD result to make the inverse"""
        rcond=self.rcond
        evals=util.FITS.Read(self.ename,savespace=1,doByteSwap=1)[1]
        ievals=numpy.where(evals<evals[0]*rcond,0.,1/evals).astype(evals.dtype)
        neig=numpy.nonzero(ievals)[0][-1]+1
        print "neig %d"%neig
        print "Normalised cancelled evals:",evals[neig:]/evals[0]
        del(evals)
        #Now work out how many stages to do the multiply in...
        #ie optimal memory use...
        lines=open("/proc/meminfo").readlines()
        for line in lines:
            if "MemTotal:" in line:
                mem=int(line.split()[1])
                multiplier={"kB":1024,"b":1,"B":1,"mB":1024*1024,"MB":1024**2,"GB":1024**3,"gB":1024**3}
                if line.split()[2] in multiplier.keys():
                    mem*=multiplier[line.split()[2]]
                else:
                    print "WARNING - multiplier %s not known for memory"
                print "Total system memory %d bytes"%mem
                break

        invmem=ievals.shape[0]**2*ievals.itemsize#memory to store the inverse - and also the memory to store u and vt as well, since all same size.
        memleft=mem-invmem
        if memleft>0:#compute the number of blocks to divide u and vt into.
            nblock=int(numpy.ceil(invmem*2./memleft))
        else:#??? just try doing the whole thing in memory!
            nblock=1
        print "Doing gemm in %d quadrants"%(nblock*nblock)
        size=ievals.shape[0]/nblock
        veut=None
        t0=time.time()
        for i in range(nblock):
            vstart=i*size
            vend=(i+1)*size
            if i==nblock-1:
                vend=ievals.shape[0]
            vt=util.FITS.Read(self.vtname,savespace=1,doByteSwap=1)[1]#c contig
            #tmp=vt[:,vstart:vend].copy()
            #del(vt)
            #vt=tmp
            #del(tmp)
            vt=vt[:neig,vstart:vend].transpose()#not contig
            v=vt.copy()#c contig
            del(vt)
            for k in xrange(neig):
                v[:,k]*=ievals[k]#multiply column by eval.
            for j in range(nblock):
                ustart=j*size
                uend=(j+1)*size
                if j==nblock-1:
                    uend=ievals.shape[0]
                u=util.FITS.Read(self.uname,savespace=1,doByteSwap=1)[1]
                tmp=u[ustart:uend,:neig].copy()#c contig
                del(u)
                ut=tmp.transpose()#fcontig
                del(tmp)

                if type(veut)==type(None):
                    veut=numpy.zeros((ievals.shape[0],ievals.shape[0]),ut.dtype,order="F")
                print "Starting gemm %d %d"%(i,j)
                t1=time.time()
                mkl.gemm(v,ut,veut[vstart:vend,ustart:uend])
                t2=time.time()
                print "GEMM time %gs"%(t2-t1)
                del(ut)
            del(v)
        t2=time.time()
        #if veut.dtype!=numpy.float32:
        #    veut=veut.astype(numpy.float32)
        util.FITS.Write(veut,self.invname,doByteSwap=0)
        print "Total GEMM time %gs"%(t2-t0)
        open(self.timename,"a").write("makeInv %gs\n"%(t2-t0))
        self.log("rcond=%g"%self.rcond)
        return veut

    def makeLUInv(self,a=None,dtype=numpy.float32):
        """Do LU decomposition to compute inverse..."""
        if type(a)!=numpy.ndarray:
            if type(a)==type(""):
                fname=a
            else:
                fname=self.dottedname
            res=util.FITS.Read(fname,doByteSwap=1,savespace=1)
            a=res[1]
            if res[0]["parsed"].has_key("REGULARD"):
                print "Using regularised matrix (%g added to diag)"%float(res[0]["parsed"]["REGULARD"])
        else:
            print "Will overwrite a with inv(a)"
        if not a.flags.c_contiguous:
            print "Transposing a before LU inversion to make c contiguous"
            a=a.transpose()#should now be c-contiguous.
        if not a.flags.c_contiguous:
            raise Exception("A not contiguous")
        if a.dtype!=dtype:
            a=a.astype(dtype)
        if acml or atlas:
            ipiv=numpy.zeros((min(a.shape),),numpy.int32)
        else:
            ipiv=numpy.zeros((min(a.shape),),numpy.int64)
        t1=time.time()
        mkl.ludecomp(a,ipiv)
        if acml or atlas:
            lw=1
        else:
            lw=mkl.luinv(a,ipiv,None)
        work=numpy.zeros((lw,),a.dtype)
        mkl.luinv(a,ipiv,work)
        t2=time.time()
        util.FITS.Write(a,self.invname,doByteSwap=0)
        self.log("makeLUInv time %g\n"%(t2-t1))
        return a


    def denseDotDense(self,a=None,b=None,transA=0,transB=0,resname=None):
        """Dot product of the poke matrix with the gen inv of pmx dot pmxT to give the reconstructor, assuming inputs and output all dense.
        a and b must be the filenames
        Can be used instead of finalDot() giving a dense result...
        """
        if a==None:
            pmxfname=self.pmxfname
        else:
            pmxfname=a
        #pmx=util.FITS.Read(pmxfname,savespace=1,doByteSwap=1)[1]
        pmxh=util.FITS.ReadHeader(pmxfname)["parsed"]
        bytepix=-int(pmxh["BITPIX"])/8
        print "Using %d byte itemsize for denseDotDense"%bytepix
        if pmxh["NAXIS"]!="2":
            raise Exception("A Matrix must be 2d")
        if transA:
            pmxshape=int(pmxh["NAXIS1"]),int(pmxh["NAXIS2"])
        else:
            pmxshape=int(pmxh["NAXIS2"]),int(pmxh["NAXIS1"])
        #if transA:
        #    pmx=pmx.transpose()
        #if not pmx.flags.c_contiguous:
        #    pmx=pmx.copy()

        if b==None:
            invname=self.invname
        else:
            invname=b
        #veut=util.FITS.Read(invname,savespace=1,doByteSwap=1)[1]
        veuth=util.FITS.ReadHeader(invname)["parsed"]
        if veuth["NAXIS"]!="2":
            raise Exception("B matrix must be 2d")
        if transB:
            veutshape=int(veuth["NAXIS1"]),int(veuth["NAXIS2"])
        else:
            veutshape=int(veuth["NAXIS2"]),int(veuth["NAXIS1"])
        #if not veut.flags.c_contiguous:
        #    veut=veut.copy()
        if resname==None:
            resname=self.rmxdensename
        mem=self.getMem()*0.9#get available memory - problems when required is just less then available - swapping - so reduce total by x0.9... 
        ansmem=pmxshape[0]*veutshape[1]*bytepix#memory to store the result
        memleft=mem-ansmem
        if memleft>0:
            nblock=int(numpy.ceil((pmxshape[0]*pmxshape[1]*bytepix+veutshape[0]*veutshape[1]*bytepix)/float(memleft)))
        else:
            print "Warning - the multiply result will fill entire memory plus - this could take a while... (%d,%d)"%(pmxshape[0],veutshape[1])
            nblock=4
        print "Doing gemm in %d quadrants.  Total shape %s x %s"%(nblock*nblock,str(pmxshape),str(veutshape))
        asize=pmxshape[0]/nblock
        bsize=veutshape[1]/nblock
        res=None
        t0=time.time()
        for i in range(nblock):
            astart=i*asize
            aend=(i+1)*asize
            if i==nblock-1:
                aend=pmxshape[0]
            try:
                a=util.FITS.Read(pmxfname,savespace=1,doByteSwap=1,memmap="r")[1]
            except:
                print "Unable to memmap %s, reading..."%pmxfname
                a=util.FITS.Read(pmxfname,savespace=1,doByteSwap=1)[1]

            if transA:
                a2=a[:,astart:aend].transpose()#not contig
                a2=a2.copy()#c contig
            else:
                a2=a[astart:aend]#c contig
            del(a)
            a=a2
            del(a2)
            for j in range(nblock):
                bstart=j*bsize
                bend=(j+1)*bsize
                if j==nblock-1:
                    bend=veutshape[1]
                try:
                    b=util.FITS.Read(invname,savespace=1,doByteSwap=1,memmap="r")[1]
                except:
                    print "Unable to memmap %s, reading..."%invname
                    b=util.FITS.Read(invname,savespace=1,doByteSwap=1)[1]
                if transB:
                    b2=b[:,bstart:bend].copy()#c contig.
                    b2=b2.transpose()#now f contig and transposed, what we want.
                else:
                    b2=b[:,bstart:bend].transpose().copy()#c contig, transposed
                    b2=b2.transpose()#now f contig... which is what we want.
                del(b)
                b=b2
                del(b2)
                if type(res)==type(None):
                    if atlas:
                        res=numpy.empty((pmxshape[0],veutshape[1]),a.dtype,order="C")
                    else:
                        res=numpy.empty((pmxshape[0],veutshape[1]),a.dtype,order="F")
                print "Starting gemm %d %d"%(i,j)
                t1=time.time()
                mkl.gemm(a,b,res[astart:aend,bstart:bend])
                t2=time.time()
                print "GEMM time %gs"%(t2-t1)
                del(b)
            del(a)
        util.FITS.Write(res,resname,doByteSwap=0)
        t1=time.time()
        print "Total GEMM time %gs"%(t2-t0)
        open(self.timename,"a").write("denseDotDense %gs\n"%(t2-t0))
        return res
        
    def dot(self,a,b,res=None,order=None):
        """A simple dot product using mkl"""
        if res==None:
            if order==None:
                if atlas:
                    order="C"
                else:
                    order="F"
            res=numpy.empty((a.shape[0],b.shape[1]),b.dtype,order=order)
        t1=time.time()
        mkl.gemm(a,b,res)
        t2=time.time()
        print "dot GEMM time %gs"%(t2-t1)
        #if order=="C":
        #    res=res.T
        return res

    def finalDot(self,veut=None,rmxtype=2,valmin=0.,save=1,maxelements=2**30,minelements=0):
        """Dot product of the poke matrix with the generalised inverse of pmx dot pmxT, to give the reconstructor.
        rmxtype can be 0 (csc), 1 (csr) or 2 (dense).  For testing purposes, dense is recommended, until you know a suitable valmin value, and is only slightly slower than a non-paged csr/csc type.
        maxelements is the max number of elements allowed.  The default is 2**30, which corresponds to 8GB (2 arrays of this size, each elemnt is 4 bytes).
        mn and mx are the min and max number of elements (fraction or whole) allowed in the final sparse matrix.  If not specified, anything is allowed.
        """
        self.log("finalDot max %g min %g"%(maxelements,minelements))
        pmxfname=self.pmxfname
        if type(veut)==type(None):
            veut=util.FITS.Read(self.invname,savespace=1,doByteSwap=1)[1]
        if not veut.flags.c_contiguous:
            veut=veut.copy()
        csc=util.FITS.loadcsc(pmxfname,doByteSwap=1)
        if maxelements!=None and maxelements<=1:
            maxelements=int(maxelements*veut.shape[0]*csc.shape[1])
        if minelements!=None and minelements<1:
            minelements=int(minelements*veut.shape[0]*csc.shape[1])
        nthreads=8

        done=0
        vallist=[]#this list will contain failed values of valmin.
        donecnt=0
        while done==0:
            print "Starting denseDotCsc with valmin %g"%valmin
            self.log("Trying denseDotCsc with valmin %g"%valmin)
            vallist.append(valmin)
            t1=time.time()
            rmx=cmod.svd.denseDotCsc(veut,csc,valmin,rmxtype,nthreads,maxelements)
            t2=time.time()
            open(self.timename,"a").write("denseDotCsc %gs\n"%(t2-t1))
            print "denseDotCsc time %gs"%(t2-t1)
            if type(rmx)==type(0):#failed when reached maxelements...
                #the value of rmx is the row reached - which will help us guess a new value for valmin.
                #Need to increase valmin.
                if rmxtype==0:
                    frac=float(rmx)*nthreads/csc.shape[1]
                else:
                    frac=float(rmx)*nthreads/veut.shape[0]
                tmp=0
                for v in vallist:
                    if v>valmin:#can't go larger than v
                        if v<tmp or tmp==0:
                            tmp=v
                if tmp==0:
                    if valmin==0:
                        valmin=frac
                    else:
                        valmin/=frac#*=1.5
                else:
                    valmin=(valmin+tmp)/2.

                print "denseDotCsc failed - return value %d - increasing valmin to %g and trying again."%(rmx,valmin)
            elif type(rmx)==type(()):#csc or csr...
                if rmx[2][-1]<minelements or rmx[2][-1]>maxelements:
                    if rmx[2][-1]<minelements:#valmin too high...
                        tmp=0
                        for v in vallist:
                            if v<valmin:#can't go lower than v
                                if v>tmp:
                                    tmp=v
                        if tmp==0:#no previous elements
                            valmin/=4.
                        else:
                            valmin=(valmin+tmp)/2.
                    else:#valmin too low (too many values)
                        tmp=0
                        for v in vallist:
                            if v>valmin:#can't go larger than v
                                if v<tmp or tmp==0:
                                    tmp=v
                        if tmp==0:
                            if valmin==0:
                                valmin=0.1
                            else:
                                valmin*=1.5
                        else:
                            valmin=(valmin+tmp)/2.
                    print "finalDot - sparse matrix created with %d elements, not in allowed range of %d to %d - trying again with valmin %g"%(rmx[2][-1],minelements,maxelements,valmin)
                else:
                    done=1
            else:
                done=1
            donecnt+=1
        print "finalDot done in %d iterations"%donecnt
        self.log("finalDot done in %d iters"%donecnt)
        cscshape=csc.shape
        veutshape=veut.shape
        del(csc)
        del(veut)

        if rmxtype==2:
            if save:
                util.FITS.Write(rmx,self.rmxdensename,doByteSwap=0)
        elif rmxtype==1:
            rmx=scipy.sparse.csr_matrix(rmx,(veutshape[0],cscshape[1]))
            if save:
                util.FITS.saveSparse(rmx,self.rmxcsrname,doByteSwap=0)
        elif rmxtype==0:
            rmx=scipy.sparse.csc_matrix(rmx,(veutshape[0],cscshape[1]))
            if save:
                util.FITS.saveSparse(rmx,self.rmxcscname,doByteSwap=0)
        return rmx

    def count(self,vals,rmx=None):
        """vals is an array (or a list) of values to be counted.  This method returns a count array, with the values equal to the number of times a number greater or equal to the corresponding vals occurs in rmx.
        """
        vals=numpy.array(vals).astype(numpy.float32)
        cnt=numpy.zeros(vals.shape,numpy.int64)
        if type(rmx)==type(None):
            rmx=util.FITS.Read(self.rmxdensename,doByteSwap=1)[1]
        t1=time.time()
        cmod.svd.countInstances(rmx,vals,cnt)
        t2=time.time()
        print "Count time %gs"%(t2-t1)
        open(self.timename,"a").write("count time %gs\n"%(t2-t1))
        f=open(self.cntname,"w")
        for i in range(cnt.shape[0]):
            f.write("%g\t%d\n"%(vals[i],cnt[i]))
        f.close()
        return cnt

    def sparsify(self,rmx,csr,val):
        """Sparsify the rmx to a csr matrix.  Here, the csr matrix shopuld be large enough to hold all values in rmx with an absolute value greater or equal to val.  This can be determined by previously calling the count method."""
        t1=time.time()
        cmod.svd.sparsifyCsr(rmx,csr,val)
        t2=time.time()
        print "Sparsify time %gs"%(t2-t1)
        open(self.timename,"a").write("sparsify time %gs\n"%(t2-t1))
        util.FITS.saveSparse(csr,self.rmxcsrname,doByteSwap=0)

    def autoSparsify(self,rmx=None,frac=0.1,fracmin=0.,vals=None):
        """Sparsify the rmx to frac sparsity (or there abouts)"""
        if rmx==None:
            rmx=self.rmxdensename
        if type(rmx)==type(""):
            print "Loading rmx"
            try:
                rmx=util.FITS.Read(rmx,memmap="r")[1]
            except:
                print "Couldn't load memmapped - trying normally..."
                rmx=util.FITS.Read(rmx,memmap="r")[1]
        size=min(int(rmx.size*frac),2**31-1)
        sizemin=min(int(rmx.size*fracmin),size)
        print "Looking for count between",sizemin,size
        val=0
        if vals==None:
            vals=numpy.array([1e-6,1e-5,1e-4,1e-3,1e-2,0.05,0.1,0.2,0.5,1.,2.,5.,10.,100.,1e3,1e4,1e6]).astype("f")
        else:
            vals=numpy.array(vals).astype("f")
        found=0
        while found==0:
            valmin=1e-100
            valmax=1e100
            print "Counting for",vals
            cnt=self.count(vals,rmx)
            print cnt
            for i in range(cnt.shape[0]):
                if cnt[i]<=size:
                    if cnt[i]>=sizemin:#have found the correct one
                        val=vals[i]
                        ndata=cnt[i]
                        found=1
                        break
                    else:
                        valmax=vals[i]
                        if i==0:
                            valmin=0.
                        else:
                            valmin=vals[i-1]
                        break
            if found==0:
                #set up a new vals array...
                print "Setting up new vals ranging from %g to %g"%(valmin,valmax)
                vals=numpy.arange(10)*(valmax-valmin)/9.+valmin
        print "Got val %g cnt %d sparsity %g (looking for %d)"%(val,ndata,ndata/float(rmx.size),size)
        if val>0:
            #csr=scipy.sparse.csr_matrix(rmx.shape,nzmax=ndata,dtype=numpy.float32)
            tdata=numpy.zeros((ndata,),numpy.float32)
            tcolind=numpy.zeros((ndata,),numpy.uint32)
            tindptr=numpy.zeros((rmx.shape[0]+1,),numpy.uint32)
            class dummy:
                data=tdata
                colind=tcolind
                indptr=tindptr
                shape=rmx.shape
                format="csr"
            d=dummy()
            #csr=scipy.sparse.csr_matrix((data,colind,indptr),shape=rmx.shape)
            #csr.data=data
            #csr.colind=colind
            #csr.indptr=indptr
            self.sparsify(rmx,d,val)
            try:
                csr=scipy.sparse.csr_matrix((tdata,tcolind,tindptr),shape=rmx.shape)
            except:
                csr=scipy.sparse.csr_matrix((tdata,tcolind,tindptr),dims=rmx.shape)
            try:
                csr.check_format()
            except:
                print "WARNING: csr.check_format() failed - method may not exist (scipy version?) or csr may be invalid."
            return csr,cnt
        else:
            print "Unable to autoSparsify"
            return None,cnt
        
        
    def getMem(self):
        lines=open("/proc/meminfo").readlines()
        for line in lines:
            if "MemTotal:" in line:
                mem=int(line.split()[1])
                multiplier={"kB":1024,"b":1,"B":1,"mB":1024*1024,"MB":1024**2,"GB":1024**3,"gB":1024**3}
                if line.split()[2] in multiplier.keys():
                    mem*=multiplier[line.split()[2]]
                else:
                    print "WARNING - multiplier %s not known for memory"
                print "Total system memory %d bytes"%mem
                break
        return mem


    def test(self,rcond,valmin=0.):
        """Only use this on small matricees"""
        csc=util.FITS.loadcsc(self.pmxfname,doByteSwap=1)
        data=csc.todense()
        dd=quick.dot(data,data.transpose())
        idd=numpy.linalg.pinv(dd,rcond)
        rmx=quick.dot(idd,data)
        rmx=numpy.where(numpy.fabs(rmx)<valmin,0,rmx).astype("f")
        return data,dd,idd,rmx

    def compress(self,rmx=None,level=None):
        """Compress the rmx.
        This can be useful for simulations where rmx is too big to fit in memory, but compressed version isn't - the sim can then uncompress each iteration (may be faster than paging...), and means that only one machine is required, rather than splitting the mvm by MPI over more than one.
        TODO later maybe.  Use compression to fit the RMX in memory, then decompress bits at a time in tomoRecon.py.
        """
        

def getMem():
    lines=open("/proc/meminfo").readlines()
    for line in lines:
        if "MemTotal:" in line:
            mem=int(line.split()[1])
            multiplier={"kB":1024,"b":1,"B":1,"mB":1024*1024,"MB":1024**2,"GB":1024**3,"gB":1024**3}
            if line.split()[2] in multiplier.keys():
                mem*=multiplier[line.split()[2]]
            else:
                print "WARNING - multiplier %s not known for memory"
            print "Total system memory %d bytes"%mem
            break
    return mem


def compress(rmx,mbits=15,outfile=None):
    """Compress a float32 array into truncated floating point consisting of bits bits per element
    Can be used to compress large reconstructors so that they fit in memory - they can then be uncompressed on the fly, a bit at a time.
    """

    import cmod.utils
    if type(rmx)==type(""):
        rmx=util.FITS.Read(rmx)[1]
    else:
        print "computeRecon.compress: Destroying input"
    r=rmx.ravel()
    mn,mx=cmod.utils.compressFloatArrayAll(r,mbits)
    offset=mx-mn+1
    ebits=1
    while (offset>>ebits)>0:
        ebits+=1
    bits=mbits+ebits+1
    words=(rmx.size*bits+31)/32
    print "Got ebits=%d, total=%d, giving %d words (shape %s)"%(ebits,bits,words,str(rmx.shape))
    r=r[:words]
    if outfile!=None:
        #don't bother byteswapping because this data isn't in a standard format anyway (eg 24 bit float???)
        util.FITS.Write(r,outfile,extraHeader=["COMPBITS= %d"%mbits,"SHAPE   = %s"%str(rmx.shape),"EXPMIN  = %d"%mn,"EXPMAX  = %d"%mx],doByteSwap=0)
    return r

if __name__=="__main__":
    if len(sys.argv)>=2:
        fname=sys.argv[1]
    else:
        fname="/var/ali/spmx42mb0.fits"
    if len(sys.argv)>=3:
        rcond=float(sys.argv[2])
    else:
        rcond=0.01
    if len(sys.argv)>=4:
        frac=float(sys.argv[3])
    else:
        frac=0.01
    t1=time.time()
    mr=makeRecon(fname,rcond)
    res=mr.dotdot()
    del(res)
    u,e,vt=mr.dosvd()
    del(u)
    del(e)
    del(vt)
    inv=mr.makeInv()
    del(inv)
    rmx=mr.finalDot(rmxtype=1,valmin=0.01,save=1,maxelements=0.1,minelements=0.05)
    #csr,cnt=mr.autoSparsify(rmx,frac)
    
    t2=time.time()
    print "Completed in %gs"%(t2-t1)
    open(mr.timename,"a").write("Total time %gs\n"%(t2-t1))
