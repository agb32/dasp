import numpy
#import util.dot as quick
class BlockMatrix:
    """A class for implementing block diagonal matricees"""
    def __init__(self,sizeList=None,dtype=numpy.float32):
        self.sizeList=sizeList
        self.dtype=dtype
        self.blockList=[]
        
        if sizeList!=None:
            self.allocate(sizeList)
        
    def allocate(self,sizeList):
        """Allocate storage for block diagonal matrix with square blocks of size specified in sizeList."""
        self.sizeList=numpy.array(sizeList)
        if len(self.sizeList.shape)==1:
            tmp=numpy.zeros((self.sizeList.shape[0],2),self.sizeList.dtype)
            tmp.transpose()[:]=self.sizeList
            self.sizeList=tmp
        self.sizeCumList=numpy.zeros((self.sizeList.shape[0]+1,2),numpy.int32)
        for i in range(1,self.sizeList.shape[0]+1):
            self.sizeCumList[i,0]=self.sizeCumList[i-1,0]+self.sizeList[i-1,0]
            self.sizeCumList[i,1]=self.sizeCumList[i-1,1]+self.sizeList[i-1,1]
        #self.data=numpy.zeros((numpy.sum(self.sizeList*self.sizeList),),self.dtype)#a 1D holder for the data.
        #self.size=numpy.sum(self.sizeList)
        self.shape=(self.sizeList[:,0].sum(),self.sizeList[:,1].sum())
        self.blockList=[]
        for y,x in self.sizeList:
            self.blockList.append(numpy.zeros((y,x),self.dtype))
        self.dtype=self.blockList[0].dtype

    def assign(self,blocklist):
        """assign blocklist to be the memory for this...
        """
        self.blockList=blocklist[:]
        self.sizeList=numpy.zeros((len(self.blockList),2),numpy.int32)
        self.sizeCumList=numpy.zeros((self.sizeList.shape[0]+1,2),numpy.int32)
        self.dtype=self.blockList[0].dtype
        for i in range(len(self.blockList)):
            self.sizeList[i,0]=self.blockList[i].shape[0]
            if len(self.blockList[i].shape)>1:
                self.sizeList[i,1]=self.blockList[i].shape[1]
            else:
                self.sizeList[i,1]=1
            self.sizeCumList[i+1,0]=self.sizeCumList[i,0]+self.sizeList[i,0]
            self.sizeCumList[i+1,1]=self.sizeCumList[i,1]+self.sizeList[i,1]
        #self.size=numpy.sum(self.sizeList)
        self.shape=(self.sizeList[:,0].sum(),self.sizeList[:,1].sum())
            
    def todense(self):
        """convert to a (large?) dense matrix"""
        if self.sizeList==None:
            return None
        arr=numpy.zeros(self.shape,self.dtype)
        totx=0
        toty=0
        for i in range(len(self.sizeList)):
            sizey,sizex=self.sizeList[i]
            arr[toty:toty+sizey,totx:totx+sizex]=self.blockList[i]
            toty+=sizey
            totx+=sizex
        return arr

    def getBlock(self,i):
        return self.blockList[i]

    def __len__(self):
        return self.shape[0]
    def __repr__(self):
        txt=""
        for block in self.blockList:
            txt+=block.__repr__()+"\n"
        return txt
    def __getitem__(self,index,withBlockIndex=0):
        #print index
        
        if type(index) in [type(1),type(1.)]:#just a single value-return a row.
            index=slice(index,index+1)
        if type(index)==type(slice(1)):#a slice...
            #need to make sure the slice doesn't span multiple blocks...
            start,stop,step=index.indices(self.shape[0])
            #print start,stop,step
            tot=0
            blockno=0
            m=min(start,stop)
            while tot<=m:
                tot+=self.sizeList[blockno,0]
                blockno+=1
            block=self.blockList[blockno-1]
            start-=tot-self.sizeList[blockno-1,0]
            stop-=tot-self.sizeList[blockno-1,0]
            if start<0 or start>block.shape[0] or stop<0 or stop>block.shape[0]:
                raise Exception("blockMatrix.py - Invalid slices - must be sliced from within a block %d %d"%(start,stop))
            if withBlockIndex:
                return block[start:stop:step],blockno-1
            else:
                return block[start:stop:step]
            
        elif type(index)==type(()):#a tuple
            if len(index)!=2:
                raise Exception("blockMatrix.py - Invalid slice index - 2D only")

            if type(index[0]) in [type(1),type(1.)]:
                ind1=slice(index[0],index[0]+1)
            else:
                ind1=index[0]
                
            if type(index[1]) in [type(1),type(1.)]:
                ind2=slice(index[1],index[1]+1)
            else:
                ind2=index[1]
            #if ind1 is valid, use it to get a sub block then apply ind2.
            #Otherwise, apply ind2 first, and then apply ind1 to that.
            transpose=0
            try:
                block,bindex=self.__getitem__(ind1,withBlockIndex=1)
            except:
                transpose=1
            if transpose:
                start,stop,step=ind2.indices(self.shape[0])
                tot=0
                blockno=0
                while tot<=start:
                    tot+=self.sizeList[blockno,0]
                    blockno+=1
                bindex=blockno-1
                block=self.blockList[bindex]
                diff=stop-start
                start-=tot-self.sizeList[bindex,0]
                block=block[:,start:start+diff:step]
                #now slice indx1 on this block.
                #start,stop,step=ind1.indices(block.shape[0])
                block=block[ind1]
            else:
                start,stop,step=ind2.indices(block.shape[1])
                block=block[:,ind2]
            return block
        elif type(index)==type(""):
            return self.blockList[int(index)]
        else:
            raise Exception("blockMatrix.py - Invalid slice index")



    def fill(self):
        """test fill a block matrix"""
        s=0
        for b in self.blockList:
            e=b.shape[0]*b.shape[0]
            b.ravel()[:]=numpy.arange(s,s+e).astype(b.dtype)
            s+=e

    def dot(self,m):
        """dot with a matrix or vector... (dense)
        Result will be dense.
        ie returns self dot m
        """
        if m.shape[0]!=self.shape[1]:
            raise Exception("Dot product of wrong shape")
        if len(m.shape)==2:
            out=numpy.zeros((self.shape[0],m.shape[1]),max(self.dtype,m.dtype))
        else:
            out=numpy.zeros((self.shape[0],),max(self.dtype,m.dtype))

        cnt=0
        for i in range(len(self.sizeList)):
            size=self.sizeList[i,1]
            out[cnt:cnt+self.sizeList[i,0]]=quick.dot(self.blockList[i],m[cnt:cnt+size])
            cnt+=size
        return out

    def dottedWith(self,m):
        """dot a matrix or vector with this... ie m dot self.
        """
        if m.shape[1]!=self.shape[0]:
            raise Exception("Dot product of wrong shape %s %s"%(str(m.shape),str(self.shape)))
        if len(m.shape)==2:
            out=numpy.zeros((m.shape[0],self.shape[1]),max(self.dtype,m.dtype))
        else:
            out=numpy.zeros((self.shape[1],),max(self.dtype,m.dtype))
        cnt=0
        for i in range(len(self.sizeList)):
            size=self.sizeList[i,0]
            #print cnt,cnt+size,m[:,cnt:cnt+size].shape,self.blockList[0].shape,out[cnt:cnt+size].shape
            out[:,cnt:cnt+self.sizeList[i,1]]=quick.dot(m[:,cnt:cnt+size],self.blockList[i])
            cnt+=size
        return out

    def dotBlock(self,m):
        """dot self with another block matrix.
        """
        if m.shape[0]!=self.shape[1]:
            raise Exception("Dot product of wrong shape")
        if numpy.alltrue(self.sizeList==m.sizeList) and numpy.alltrue(self.sizeList[:,0]==self.sizeList[:,1]):#note the second test may be unneccessary - I just odn't have time to think it through at the moment...
            #all blocks same size - just dot the blocks...
            out=BlockMatrix(self.sizeList)
            for i in range(len(self.sizeList)):
                out.blockList[i]=quick.dot(self.blockList[i],m.blockList[i])
            return out
        else:
            #output could be any size... will probably be block, but maybe not square blocks.
            out=numpy.zeros((self.shape[0],m.shape[1]),max(self.dtype,m.dtype))
            x=0
            y=0
            for i in xrange(len(m.sizeList)):
                mblock=m.blockList[i]
                for j in xrange(mblock.shape[1]):
                    for y in xrange(self.shape[0]):
                        #print self.getSlice([slice(y,y+1),slice(int(m.sizeCumList[i]),int(m.sizeCumList[i+1]))]).shape
                        #print mblock[:,j].shape
                        out[y,x]=quick.dot(self.getSlice([slice(y,y+1),slice(int(m.sizeCumList[i,1]),int(m.sizeCumList[i+1,1]))]),mblock[:,j])
                    x+=1
            return out

##             prodList=[]
##             mi=0
##             for i in range(len(self.sizeList)):
##                 b1=self.blockList[i]
##                 sfrom=self.sizeCumList[i]
##                 sto=self.sizeCumList[i+1]
##                 #Now find out which blocks of m this intersects...
##                 while m.sizeCumList[mi]<sto:#while we're still looking at parts of m in the current block...
##                     while m.sizeCumList[mi+1]<sfrom:#make sure that we're in the current block to start with (should be)...
##                         mi+=1
##                     mfrom=m.sizeCumList[mi]
##                     mto=m.sizeCumList[mi+1]
##                     if mfrom<sfrom:
##                         if mto<sto:
##                             part1=b1[:mto-sfrom,:mto-sfrom]
##                             part2=m.blockList[mi][sfrom-mfrom:,sfrom-mfrom:]
##                         else:
##                             part1=b1
##                             part2=m.blockList[mi][sfrom-mfrom:sto-mfrom,sfrom-mfrom:sto-mfrom]
##                     else:#sfrom<mfrom
##                         if mto<sto:
##                             part1=b1[mfrom-sfrom:mto-sfrom,mfrom-sfrom:mto-sfrom]
##                             part2=m.blockList[mi]
##                         else:#sto<mto
##                             part1=b1[mfrom-sfrom:,mfrom-sfrom:]
##                             part2=m.blockList[mi][:sto-mfrom,:sto-mfrom]
##                     prodList.append(quick.dot(part1,part2))
##                     mi+=1
##             out=BlockMatrix()
##             out.assign(prodList)
##             return out
    def getSlice(self,indx):
        """This is equivalent to doing:
        self.todense()[indx]
        where indx can be slices, or a single index.
        It returns a copy of the data.
        Can be called using python slice() object to represent slices.
        We don't use the normal slice method to do this, because that may
        involve copying the data, which may be unintended...
        """
        #Have failed to get the data without copying, so now need to copy the
        #data... Some off-block data has been requested (ie empty matrix). 
        single=0
        if type(indx) in [type(1),type(1.)]:
            #single row requested.
            indx=slice(indx,indx+1,1)
            single=1
            
        if type(indx)==type(slice(1)):
            start,stop,step=indx.indices(self.shape[0])
            diff=(int(numpy.fabs(stop-start))+step-1)/step
            out=numpy.zeros((diff,self.shape[1]),self.dtype)
            #print "shape",out.shape,start,stop,step
            #fill the matrix here...
            rows=range(start,stop,step)
            i=0
            for row in rows:
                #first find out which block we're in, and then which row of this block.
                j=0
                while row>=self.sizeCumList[j+1,0]:
                    j+=1
                out[i,self.sizeCumList[j,1]:self.sizeCumList[j+1,1]]=self.blockList[j][row-self.sizeCumList[j,0]]
                i+=1
                    
        elif type(indx) in [type(()),type([])]:
            if len(indx)!=2:
                raise Exception("Must be 1 or 2D slice")
            tmp=[]
            coord=[]
            diff=[]
            single=0
            for i in range(len(indx)):
                if type(indx[i]) in [type(1),type(1.)]:
                    tmp.append(slice(indx[i],indx[i]+1,1))
                    single+=1
                else:#already a slice...
                    tmp.append(indx[i])
                coord.append(tmp[i].indices(self.shape[i]))
                diff.append((int(numpy.fabs(coord[i][1]-coord[i][0]))+coord[i][2]-1)/coord[i][2])
            indx=tmp
            if single!=2:#if they are both single values, requesting a single value.  Otherwise, will be an array.
                single=0
            out=numpy.zeros((diff[0],diff[1]),self.dtype)
            #now fill the matrix.
            rows=range(coord[0][0],coord[0][1],coord[0][2])
            i=0
            cfrom,cto,cstep=map(int,coord[1])
            cols=range(cfrom,cto,cstep)
            for row in rows:
                #get block, and row of this block.
                j=0
                while row>=self.sizeCumList[j+1,0]:
                    j+=1
                tmp=self.blockList[j][row-self.sizeCumList[j,0]]
                nstart=None
                #now select the part of this row - is the block used at all?
                if cfrom<self.sizeCumList[j,1]:
                    if cto<=self.sizeCumList[j,1]:#block not used
                        pass
                    else:#block is used (if step allows)...
                        #Find the first index of use:
                        #ie sizeCumList[j]+index==cfrom+n*cstep and index<cstep
                        indfrom=(self.sizeCumList[j,1]-cfrom)%cstep
                        if indfrom!=0:
                            indfrom=cstep-indfrom
                        nstart=(self.sizeCumList[j,1]+indfrom-cfrom+cstep-1)/cstep
                        #now find the stopping points...
                        nused=(min(cto,self.sizeCumList[j+1,1])-(self.sizeCumList[j,1]+indfrom)+cstep-1)/cstep
                        indto=indfrom+cstep*nused
                        nend=nstart+nused
                        if indfrom>=tmp.shape[0]:#step misses the block!
                            nstart=None
                        
                elif cfrom<self.sizeCumList[j+1,1]:#block used
                    indfrom=cfrom-self.sizeCumList[j,1]
                    nused=(min(cto,self.sizeCumList[j+1,1])-cfrom+cstep-1)/cstep
                    indto=indfrom+cstep*nused
                    nstart=0
                    nend=nused
                else:#block not used
                    pass
                if nstart!=None:
                    out[i,nstart:nend]=tmp[indfrom:indto:cstep]
                i+=1

        if single==2:
            out=out[0,0]
        elif single==1:
            out=out[0]
        return out

    def copy(self):
        """Create a copy"""
        bl=[]
        for b in self.blockList:
            bl.append(b.copy())
        b=BlockMatrix()
        b.assign(bl)
        return b
    def inv(self):
        """Inverse and return as new matrix."""
        if not numpy.alltrue(self.sizeList[:,0]==self.sizeList[:,1]):
            raise Exception("Cannont inverse a non-square block matrix - please use todense and do it manually")
        newBM=BlockMatrix()
        bl=[]
        for b in self.blockList:
            bl.append(numpy.linalg.inv(b))
        newBM.assign(bl)
        return newBM
    def pinv(self,rcond=1e-15):
        """Inverse and return as new matrix."""
        if not numpy.alltrue(self.sizeList[:,0]==self.sizeList[:,1]):
            raise Exception("Cannont inverse a non-square block matrix - please use todense and do it manually")
        newBM=BlockMatrix()
        bl=[]
        for b in self.blockList:
            bl.append(numpy.linalg.pinv(b,rcond))
        newBM.assign(bl)
        return newBM
    def rand(self):
        """fill with random values."""
        for b in self.blockList:
            b[:]=numpy.random.random(b.shape)
    def transpose(self):
        """Return the transpose of self."""
        if not numpy.alltrue(self.sizeList[:,0]==self.sizeList[:,1]):
            raise Exception("Not yet implemented - transpose of block matrix with non-square blocks - should be easy to do if you need it...")
        newBM=self.copy()
        for b in newBM.blockList:
            b[:]=b.copy().transpose()
            
def dot(mx1,mx2):
    """Take the dot prod of 2 matrices, which may or may not be block
    """
    if isinstance(mx1,BlockMatrix):
        if isinstance(mx2,BlockMatrix):
            return mx1.dotBlock(mx2)
        else:
            return mx1.dot(mx2)
    else:
        if isinstance(mx2,BlockMatrix):
            return mx2.dottedWith(mx1)
        else:
            return quick.dot(mx1,mx2)
