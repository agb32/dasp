import svd,Numeric,random,gist
rows=10
cols=10
ndata=50
a=Numeric.zeros((rows,cols),"d")
data=Numeric.zeros((ndata,),"d")
rowind=Numeric.zeros((ndata,),"i")
indptr=Numeric.zeros((cols+1,),"i")
sp=svd.sparseMatrixCreate(ndata,rows,cols,data,rowind,indptr)
b=Numeric.zeros((rows,cols),"d")
for y in xrange(rows):
    for x in xrange(cols):
        r=random.random()*2-1
        a[y,x]=r
        tmp=svd.sparseMatrixInsert(sp,y,x,r)
        print "Inserted (%d) %g at (%d, %d), got %g"%(tmp,r,y,x,svd.sparseMatrixGet(sp,y,x))
        print data
        print rowind
        print indptr
        print svd.sparseMatrixInfo(sp),"(ndata,rows,cols,min,rowmin,indmin,cnt,alloced)"
        raw_input()

for y in xrange(rows):
    for x in xrange(cols):
        b[y,x]=svd.sparseMatrixGet(sp,y,x)


val=Numeric.sort(Numeric.fabs(a.flat))[-ndata]
c=Numeric.where(Numeric.fabs(a)<val,0,a)

gist.window(0);gist.fma();gist.pli(a)
gist.window(1);gist.fma();gist.pli(b)
gist.window(2);gist.fma();gist.pli(c-b)
#At this stage, b and c should be identical.

for i in xrange(1000):
    x=random.randint(0,cols-1)
    y=random.randint(0,rows-1)
    r=random.random()*2-1
    a[y,x]=r
    tmp=svd.sparseMatrixInsert(sp,y,x,r)
    print "Inserted (%d) %g at (%d, %d), got %g"%(tmp,r,y,x,svd.sparseMatrixGet(sp,y,x))
    print data
    print rowind
    print indptr
    print svd.sparseMatrixInfo(sp),"(ndata,rows,cols,min,rowmin,indmin,cnt,alloced)"

#now repeat the plotting of the b/c data... note, they may not be the same any more.

#now test setting all to zeros.
for y in xrange(rows):
    for x in xrange(cols):
        tmp=svd.sparseMatrixInsert(sp,y,x,0.)
        print "Inserted (%d) %g at (%d, %d), got %g"%(tmp,0,y,x,svd.sparseMatrixGet(sp,y,x))
        print data
        print rowind
        print indptr
        print svd.sparseMatrixInfo(sp),"(ndata,rows,cols,min,rowmin,indmin,cnt,alloced)"

#now repeat the plotting of the b data - should all be zeros.

