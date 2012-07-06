"""Module to create theoretical poke matrix for Fried geometry.  Can forsee several different methods here, for different types of pokes.  eg poking individual actuators, poking as zernikes etc."""
import numpy as na
#import util.dot as quick
#na.float32="f"#hack for Numeric - not needed for numpy.
import scipy.sparse
def straightPokes(nact,vlist=[1.],pupfn=None,actmap=None,returnfull=0):
    """Create poke matrix based on individual pokes.
    This assumes fried geometry, and nsubx==nact-1.  If you want other
    geometries, you'll have to implement them... (sorry - won't be so trivial!)
    nact is number of actuators in 1d.
    vlist is a list of values to place in the pokematrix.  The first value
    in the list is the value to give to centroids immediately next to an
    actuator, the second value is for centroids one away from the actuator
    etc.  At the moment, only the first element of vlist is used, ie only
    subapertures immediately next to the centroid are affected.  This is
    probably a reasonable approximation unless you want a really detailed
    mirror model, inwhich case, you should add to this.
    pupfn is the pupil function, array size (nact-1,nact-1), which defines
    which centroid values are used and ignored (if None, all used).
    actmap is array size (nact,nact) defining which actuators are poked.
    """
    if actmap==None:#note, numpy gets the comparison with None correct.
        actmap=na.ones((nact,nact))
    if pupfn==None:
        pupfn=na.ones((nact-1,nact-1))
    ncent=na.sum(pupfn.flat)*2
    nacts=na.sum(actmap.flat)
    pupcnt=pupfn.copy()
    pos=0
    #Index the pupcnt array... (with centroid count).
    for i in range(nact-1):
        for j in range(nact-1):
            if pupcnt[i,j]==1:
                pupcnt[i,j]=pos
                pos+=1
            else:
                pupcnt[i,j]=-1
    pokemx=na.zeros((nacts,ncent),na.float32)
    #probably easiest to create a whole pokemx, and strip out the actuators
    #and centroids that aren't needed.
    pokemxfull=na.zeros((nact*nact,(nact-1)*(nact-1)*2),na.float32)
    pokemxtmp=na.zeros((nacts,(nact-1)*(nact-1)*2),na.float32)
    nsubx=nact-1
    nc=nsubx*nsubx
    for i in range(nact):
        for j in range(nact):
            actpos=i*nact+j
            if i>0:
                if j>0:
                    pokemxfull[actpos,(i-1)*nsubx+j-1]=-vlist[0]
                    pokemxfull[actpos,(i-1)*nsubx+j-1+nc]=-vlist[0]
                if j<nsubx:
                    pokemxfull[actpos,(i-1)*nsubx+j]=vlist[0]
                    pokemxfull[actpos,(i-1)*nsubx+j+nc]=-vlist[0]
            if i<nsubx:
                if j>0:
                    pokemxfull[actpos,i*nsubx+j-1]=-vlist[0]
                    pokemxfull[actpos,i*nsubx+j-1+nc]=vlist[0]
                if j<nsubx:
                    pokemxfull[actpos,i*nsubx+j]=vlist[0]
                    pokemxfull[actpos,i*nsubx+j+nc]=vlist[0]
    #now get rid of unwanted actuators
    actpos=0
    for i in range(nact):
        for j in range(nact):
            if actmap[i,j]==1:
                pokemxtmp[actpos]=pokemxfull[i*nact+j]
                actpos+=1
            else:
                pokemxfull[i*nact+j]=0
    #and get rid of unwanted centroids
    centpos=0
    for i in range(nact-1):
        for j in range(nact-1):
            if pupfn[i,j]==1:
                pokemx[:,centpos]=pokemxtmp[:,i*nsubx+j]
                pokemx[:,centpos+ncent/2]=pokemxtmp[:,i*nsubx+j+nc]
                centpos+=1
            else:
                pokemxfull[:,i*nsubx+j]=0
                pokemxfull[:,i*nsubx+j+nc]=0
    if returnfull:
        pokemx=pokemxfull
    return pokemx


class sparsePMX:
    """A class to represent a poke matrix (fried geometry, zonal) with a
    sparse representation.  You can check that your poke matrix has been
    converted to sparse format okay by using the bulkify method - this should
    give you back your poke matrix"""
    def __init__(self,nact,pmx=None):
        self.nact=nact
        self.nact2=nact*nact
        self.nsubx=nact-1
        self.ncent=self.nsubx**2*2
        self.storage=None
        self.val=None
        if type(pmx)!=type(None):
            self.sparsify(pmx)
    
    def sparsify(self,pmx):
        """converts a pokematrix to sparse format, and returns such an object.
        This object has a dot method which can be used for dotting with
        things."""
        self.storage=na.zeros((8,self.nsubx,self.nsubx),pmx.dtype.char)
        for i in range(self.nsubx):#block...
            for j in range(self.nsubx):#element in block...
                self.storage[0,i,j]=pmx[i*self.nact+j,j+i*self.nsubx]
                self.storage[1,i,j]=pmx[i*self.nact+j+1,j+i*self.nsubx]
                self.storage[2,i,j]=pmx[i*self.nact+j+self.nact,j+i*self.nsubx]
                self.storage[3,i,j]=pmx[i*self.nact+j+self.nact+1,j+i*self.nsubx]
                self.storage[4,i,j]=pmx[i*self.nact+j,j+i*self.nsubx+self.ncent/2]
                self.storage[5,i,j]=pmx[i*self.nact+j+1,j+i*self.nsubx+self.ncent/2]
                self.storage[6,i,j]=pmx[i*self.nact+j+self.nact,j+i*self.nsubx+self.ncent/2]
                self.storage[7,i,j]=pmx[i*self.nact+j+self.nact+1,j+i*self.nsubx+self.ncent/2]

    def create(self,val=1.,dtype="f"):
        """Create a sparse poke matrix.  If this has worked, self.bulkify() should return the same as a straightPokes() poke matrix.
        This avoids creating the poke matrix, which sometimes can be too large for memory.
        Assumes that all actuators and subaps are used (the default case for straightpokes).
        """
        self.val=val
        self.storage=na.zeros((8,self.nsubx,self.nsubx),dtype)
        self.storage[0]=val
        self.storage[1]=-val
        self.storage[2]=val
        self.storage[3]=-val
        self.storage[4]=val
        self.storage[5]=val
        self.storage[6]=-val
        self.storage[7]=-val
        

    def bulkify(self,res=None):
        """converts a sparse matrix into a full matrix..."""
        if type(res)==type(None):
            res=na.zeros((self.nact2,self.ncent),self.storage.dtype.char)
        if res.shape!=(self.nact2,self.ncent):
            raise Exception("result for bulkify is wrong shape")
        for i in range(self.nsubx):#block...
            for j in range(self.nsubx):#element in block...
                res[i*self.nact+j,j+i*self.nsubx]=self.storage[0,i,j]
                res[i*self.nact+j+1,j+i*self.nsubx]=self.storage[1,i,j]
                res[i*self.nact+j+self.nact,j+i*self.nsubx]=self.storage[2,i,j]
                res[i*self.nact+j+self.nact+1,j+i*self.nsubx]=self.storage[3,i,j]
                res[i*self.nact+j,j+i*self.nsubx+self.ncent/2]=self.storage[4,i,j]
                res[i*self.nact+j+1,j+i*self.nsubx+self.ncent/2]=self.storage[5,i,j]
                res[i*self.nact+j+self.nact,j+i*self.nsubx+self.ncent/2]=self.storage[6,i,j]
                res[i*self.nact+j+self.nact+1,j+i*self.nsubx+self.ncent/2]=self.storage[7,i,j]
        return res
    def dot(self,mx,res=None):
        """dot product of self with mx (1d or 2d).  Optionally, put the result
        in res... (stops a new array being allocated).
        Result should be the same as numpy.dot(self.bulkify(),mx)"""
        if mx.shape[0]!=self.ncent:
            raise Exception("mx has wrong shape for multiplication with sparse poke matrix")
        dim=1
        nsubx=self.nsubx
        nact=self.nact
        nc=self.ncent/2
        resshape=(self.nact2,)
        if len(mx.shape)==2:
            dim=2
            resshape=(self.nact2,mx.shape[1])
        if type(res)==type(None):
            res=na.zeros(resshape,self.storage.dtype.char)
        else:
            res[:,]=0
        if res.shape[0]!=self.nact2:
            raise Exception("result matrix has wrong shape for sparse poke matrix multiplication")
        if dim==1:
            self._dot1d(mx,res)
        elif dim==2:
            self._dot2d(mx,res)
        return res

    def _dot2d(self,mx,res):
        """internal work function - does a 2d dot product (matrix.matrix)"""
        for i in range(mx.shape[1]):#iter over columns...
            self._dot1d(mx[:,i],res[:,i])
        
    
    def _dot1d(self,mx,res):
        """internal work function... does a 1D dot product (matrix.vector)"""
        nsubx=self.nsubx
        nact=self.nact
        nc=self.ncent/2
        for i in range(nsubx):
            s=i*nsubx
            e=s+nsubx
            rs=i*nact
            re=rs+nsubx
            s2=s+nc
            e2=e+nc
            res[rs:re]+=self.storage[0,i]*mx[s:e]+self.storage[4,i]*mx[s2:e2]
            res[rs+1:re+1]+=self.storage[1,i]*mx[s:e]+self.storage[5,i]*mx[s2:e2]
            res[rs+nact:re+nact]+=self.storage[2,i]*mx[s:e]+self.storage[6,i]*mx[s2:e2]
            res[rs+nact+1:re+nact+1]+=self.storage[3,i]*mx[s:e]+self.storage[7,i]*mx[s2:e2]

    def dotWithSparsePmx(self,spmxT):
        """Basically does pmx dot pmxT.  Except that pmx is self, and pmxT is a sparse PMX instance (may also be self) which hasn't really been transposed.
        Returns a sparse matrix in the form of a dictionary of 1D arrays, which are the diagonals of the sparse system. (The sparsity is such that only some of the diagonals have data).
        """
        prod={0:na.zeros((self.nact2,),na.float32)}#the diagonal
        prod[1]=na.zeros((self.nact2,),na.float32)#one away from the diagonal
        prod[-1]=na.zeros((self.nact2,),na.float32)#one away from the diagonal
        prod[self.nact]=na.zeros((self.nact2,),na.float32)
        prod[-self.nact]=na.zeros((self.nact2,),na.float32)
        prod[self.nact-1]=na.zeros((self.nact2,),na.float32)
        prod[-(self.nact-1)]=na.zeros((self.nact2,),na.float32)
        prod[self.nact+1]=na.zeros((self.nact2,),na.float32)
        prod[-(self.nact+1)]=na.zeros((self.nact2,),na.float32)
        nsubx=self.nsubx
        m=self.storage
        n=spmxT.storage
        for y in range(self.nact):
            for x in range(self.nact):
                i=y*self.nact+x
                if y<nsubx and x<nsubx:
                    prod[0][i]+=m[0,y,x]*n[0,y,x]
                    prod[0][i]+=m[4,y,x]*n[4,y,x]
                    prod[1][i]+=m[0,y,x]*n[1,y,x]
                    prod[1][i]+=m[4,y,x]*n[5,y,x]
                    prod[self.nact][i]+=m[0,y,x]*n[2,y,x]
                    prod[self.nact][i]+=m[4,y,x]*n[6,y,x]
                    prod[self.nact+1][i]+=m[0,y,x]*n[3,y,x]
                    prod[self.nact+1][i]+=m[4,y,x]*n[7,y,x]
                if y<nsubx and x>0:
                    prod[0][i]+=m[1,y,x-1]*n[1,y,x-1]
                    prod[0][i]+=m[5,y,x-1]*n[5,y,x-1]
                    prod[-1][i]+=m[1,y,x-1]*n[0,y,x-1]
                    prod[-1][i]+=m[5,y,x-1]*n[4,y,x-1]
                    prod[self.nact][i]+=m[1,y,x-1]*n[3,y,x-1]
                    prod[self.nact][i]+=m[5,y,x-1]*n[7,y,x-1]
                    prod[self.nact-1][i]+=m[1,y,x-1]*n[2,y,x-1]
                    prod[self.nact-1][i]+=m[5,y,x-1]*n[6,y,x-1]
                if y>0 and x<nsubx:
                    prod[0][i]+=m[2,y-1,x]*n[2,y-1,x]
                    prod[0][i]+=m[6,y-1,x]*n[6,y-1,x]
                    prod[1][i]+=m[2,y-1,x]*n[3,y-1,x]
                    prod[1][i]+=m[6,y-1,x]*n[7,y-1,x]
                    prod[-self.nact][i]+=m[2,y-1,x]*n[0,y-1,x]
                    prod[-self.nact][i]+=m[6,y-1,x]*n[4,y-1,x]
                    prod[-(self.nact-1)][i]+=m[2,y-1,x]*n[1,y-1,x]
                    prod[-(self.nact-1)][i]+=m[6,y-1,x]*n[5,y-1,x]
                if y>0 and x>0:
                    prod[0][i]+=m[3,y-1,x-1]*n[3,y-1,x-1]
                    prod[0][i]+=m[7,y-1,x-1]*n[7,y-1,x-1]
                    prod[-1][i]+=m[3,y-1,x-1]*n[2,y-1,x-1]
                    prod[-1][i]+=m[7,y-1,x-1]*n[6,y-1,x-1]
                    prod[-self.nact][i]+=m[3,y-1,x-1]*n[1,y-1,x-1]
                    prod[-self.nact][i]+=m[7,y-1,x-1]*n[5,y-1,x-1]
                    prod[-(self.nact+1)][i]+=m[3,y-1,x-1]*n[0,y-1,x-1]
                    prod[-(self.nact+1)][i]+=m[7,y-1,x-1]*n[4,y-1,x-1]
                
        prod[1]=prod[1][:self.nact2-1]
        prod[-1]=prod[-1][1:,]
        prod[self.nact]=prod[self.nact][:self.nact2-self.nact]
        prod[-self.nact]=prod[-self.nact][self.nact:,]
        prod[self.nact+1]=prod[self.nact+1][:self.nact2-(self.nact+1)]
        prod[-(self.nact+1)]=prod[-(self.nact+1)][(self.nact+1):,]
        prod[self.nact-1]=prod[self.nact-1][:self.nact2-(self.nact-1)]
        prod[-(self.nact-1)]=prod[-(self.nact-1)][(self.nact-1):,]
        return prod
    
def tomographicPokeMatrix(ngsList,lgsList,dmList,telDiam,target):
    """Creates a tomographic poke matrix, ie one using several DMs, and several guide stars etc.
    ngsList is a list of NGS objects, each of which includes information about direction, wfs order.
    lgsList is a list of LGS objects, each of which includes information about direction and height, wfs order.
    dmList is a list of DM objects, each of which includes information about conjugate height (and maybe a virtual flag?), actuator number.
    telDiam is the telescope diameter.
    target includes information about target direction (poke matrix will be created for this direction).
    """
    nacts=0
    ncent=0
    for gs in ngsList+lgsList:
        ncent+=2*(gs.nsubx**2)
    for dm in dmList:
        nacts+=dm.nact**2
    
    pokemxfull=na.zeros((nacts,ncent),na.float32)#create whole pokemx, then strip out unwanted parts later.
    actcnt=0
    for dm in dmList:
        dmDiam=na.tan(dm.fov/3600./180*na.pi)*dm.height*2+telDiam#width of the physical dm.
        actDiam=dmDiam/(dm.nact-1)
        dmGradient=(1-dm.coupling)/actDiam
        print dmGradient,dmDiam,actDiam
        centCnt=0
        for gs in ngsList+lgsList:
            subapDiam=telDiam/gs.nsubx
            r=dm.height*na.tan(gs.theta/3600./180.*na.pi)
            wsx=r*na.cos(gs.phi/180.*na.pi)#position of centre of phase on dm.
            wsy=r*na.sin(gs.phi/180.*na.pi)
            wsx+=dmDiam/2.-telDiam/2.#wavefront sensor start coordinates.
            wsy+=dmDiam/2.-telDiam/2.
            print wsx,wsy,dmDiam,dm.height,gs.theta,gs.phi
            if wsx<0 or wsy<0 or wsx+telDiam>dmDiam or wsy+telDiam>dmDiam:
                raise Exception("WFS is outside DM fov %g %g %g"%(dmDiam,wsx,wsy))
            for i in xrange(dm.nact):
                ay=i*actDiam#actuator y position
                wy=ay-wsy#position on wfs - 0 to subapDiam means within first subap etc... subapDiam/2 means the centre of the subap.
                wymin=wy-actDiam#the possible range of influence of this actuator (will be less if coupling<0.
                wymax=wy+actDiam
                if wymax>=-0.5*subapDiam and wymin<=(gs.nsubx+0.5)*subapDiam:
                    # poke may be measurable...
                    for j in xrange(dm.nact):
                        # work out which centroids for this gs will be affected, assume first order only.
                        # is the actuator in the field of view of the wfs?
                        # Assume that the centre of each DM is on-axis.  Note, that actuators are assumed to start at the edge of the DM.
                        #act[i,j] is being poked here of deformable mirror dm.

                        ax=j*actDiam#actuator x position
                        # compute position on wfs...
                        wx=ax-wsx
                        wxmin=wx-actDiam
                        wxmax=wx+actDiam
                        centList=[]
                        if wxmax>=-0.5*subapDiam and wxmin<=(gs.nsubx+0.5)*subapDiam:
                            #poke will have an affect on the wfs.
                            for sy in range(int(na.floor(wymin/subapDiam)),int(na.ceil(wymax/subapDiam))):
                                for sx in range(int(na.floor(wxmin/subapDiam)),int(na.ceil(wxmax/subapDiam))):
                                    #the affected centroids...
                                    csy=(sy+0.5)*subapDiam#centre of subap coordinates.
                                    csx=(sx+0.5)*subapDiam
                                    dy=csy-wy
                                    dx=csx-wx
                                    valy=1-dmGradient*na.fabs(dy)
                                    valx=1-dmGradient*na.fabs(dx)
                                    if valy<0:
                                        valy=None
                                    if valx<0:
                                        valx=None
                                    if dy<0 and valy!=None:
                                        valy*=-1
                                    if dx<0 and valx!=None:
                                        valx*=-1
                                    centList.append((sy,sx,valy,valx))

                        for (sy,sx,valy,valx) in centList:
                            if sy>=0 and sy<gs.nsubx and sx>=0 and sx<gs.nsubx:
                                actnumber=actcnt+i*dm.nact+j
                                centnumber=centCnt+sy*gs.nsubx+sx
                                if valx!=None:
                                    pokemxfull[actnumber,centnumber]=valx
                                if valy!=None:
                                    pokemxfull[actnumber,centnumber+ncent/2]=valy
                                #print i,j,sy,sx,actnumber,centnumber,valx,valy,wxmin,wxmax,wymin,wymax
            #finished with this wfs for now...
            centCnt+=gs.nsubx**2
        #and now we're finished with this DM.
        actcnt+=dm.nact**2
    return pokemxfull
        
class DM:
    """This is now depreciated - use util.dm.physicalDM instead."""
    def __init__(self,nact,height,fov,coupling=0.1):
        self.fov=fov#field of view in arcsec (diameter, not radius).  Probably equal to GS.theta*2.
        self.height=height#conjugate height in m.
        self.nact=nact#number of actuators.
        self.coupling=coupling#coupling between actuators (rudimentary)
        print "createPokeMx.DM - please use util.dm.physicalDM instead"
    def computeCoords(self,telDiam,fname=None):
        """Creates an array containing coords of each actuator.
        0,0 means on axis.
        You can use gnuplot to plot the relative coordinates of the actuators and sub-aperture centres, e.g. using:
        set style line 2 lt 2#green
        plot "filename.csv" every :1::0::0 using 1:2,"" every :1::1::1 using 1:2 ls 2
        This assumes that there are 2 sets of coordinates in the file (e.g. one DM and one LGS).
        """
        self.coords=na.zeros((self.nact,self.nact,2),na.float32)
        dmDiam=na.tan(self.fov/3600./180*na.pi)*self.height*2+telDiam#width of the physical dm.
        actDiam=dmDiam/(self.nact-1)
        if fname!=None:
            f=open(fname,"a")
            f.write("#DM height %g, fov %g, nact %g\n"%(self.height,self.fov,self.nact))
        for i in range(self.nact):
            for j in range(self.nact):
                self.coords[i,j,0]=-dmDiam/2+j*actDiam
                self.coords[i,j,1]=-dmDiam/2+i*actDiam
                if fname!=None:
                    f.write("%g\t%g\n"%(self.coords[i,j,0],self.coords[i,j,1]))
        if fname!=None:
            f.write("\n")
            f.close()

class LGS:
    def __init__(self,nsubx,theta,phi):
        self.theta=theta#angle in arcsec
        self.phi=phi#angle in degrees
        self.nsubx=nsubx
        print "createPokeMx.LGS - please use util.guideStar.LGS instead"
    def computeCoords(self,telDiam,height,fname=None):
        """Creates an array containing coords of centre of each subap.
        Height is the height of interest (eg the DM conjugate height).
        0,0 means on axis."""
        self.coords=na.zeros((self.nsubx,self.nsubx,2),na.float32)
        r=height*na.tan(self.theta/3600./180.*na.pi)
        wsx=r*na.cos(self.phi/180.*na.pi)#position of centre of phase on dm.
        wsy=r*na.sin(self.phi/180.*na.pi)
        subapDiam=telDiam/self.nsubx
        if fname!=None:
            f=open(fname,"a")
            f.write("#GS height %g, nsubx %g, theta %g, phi %g\n"%(height,self.nsubx,self.theta,self.phi))
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                self.coords[i,j,0]=wsx-telDiam/2.+j*subapDiam+subapDiam/2
                self.coords[i,j,1]=wsy-telDiam/2.+i*subapDiam+subapDiam/2
                if fname!=None:
                    f.write("%g\t%g\n"%(self.coords[i,j,0],self.coords[i,j,1]))
        if fname!=None:
            f.write("\n")
            f.close()
class NGS:
    def __init__(self,nsubx,theta,phi):
        self.theta=theta#angle in arcsec
        self.phi=phi#angle in degrees
        self.nsubx=nsubx
        print "createPokeMx.NGS - please use util.guideStar.NGS instead"
    def computeCoords(self,telDiam,height,fname=None):
        """Creates an array containing coords of centre of each subap.
        Height is the height of interest (eg the DM conjugate height).
        0,0 means on axis."""
        self.coords=na.zeros((self.nsubx,self.nsubx,2),na.float32)
        r=height*na.tan(self.theta/3600./180.*na.pi)
        wsx=r*na.cos(self.phi/180.*na.pi)#position of centre of phase on dm.
        wsy=r*na.sin(self.phi/180.*na.pi)
        subapDiam=telDiam/self.nsubx
        if fname!=None:
            f=open(fname,"a")
            f.write("#GS height %g, nsubx %g, theta %g, phi %g\n"%(height,self.nsubx,self.theta,self.phi))
        for i in range(self.nsubx):
            for j in range(self.nsubx):
                self.coords[i,j,0]=wsx-telDiam/2.+j*subapDiam+subapDiam/2
                self.coords[i,j,1]=wsy-telDiam/2.+i*subapDiam+subapDiam/2
                if fname!=None:
                    f.write("%g\t%g\n"%(self.coords[i,j,0],self.coords[i,j,1]))
        if fname!=None:
            f.write("\n")
            f.close()

def testTomoSingle():
    """A single DM conjugated to ground layer"""
    dm=DM(9,0,0,0.1)
    gs=LGS(8,0,0)
    pmx=tomographicPokeMatrix([gs],[],[dm],42.,None)
    return pmx
def testTomo2():
    dm=DM(9,0,0,0.1)
    dm2=DM(9,2000.,1200.,0.1)
    gs=LGS(8,0,0)
    gs2=LGS(8,500.,0)
    pmx=tomographicPokeMatrix([gs,gs2],[],[dm,dm2],42.,None)
    return pmx
def testTomo3():
    """Poke matrix for a large system"""
    dm=[DM(9,0,0)]
    dm.append(DM(9,200.,150.))
    dm.append(DM(9,2000.,150.))
    dm.append(DM(9,10000.,150.))
    gs=[LGS(8,0.,0.)]
    gs.append(LGS(8,60.,0.))
    gs.append(LGS(8,60.,120.))
    gs.append(LGS(8,60.,240.))
    pmx=tomographicPokeMatrix(gs,[],dm,42.,None)
    return pmx
def testTomo4(sparse=0):
    dm=[DM(9,0,0)]
    dm.append(DM(5,200.,150.))
    gs=[LGS(4,0.,0.)]
    gs.append(LGS(8,60.,0.))
    gs.append(LGS(8,60.,120.))
    gs.append(LGS(8,60.,240.))
    if sparse:
        pmx=sparseTomo(ngsList=gs,lgsList=[],dmList=dm)
    else:
        pmx=tomographicPokeMatrix(gs,[],dm,42.,None)
    return pmx
    
def testTomoEagle(sparse=0):
    """Possibly poke matrix for EAGLE. 40GB in size."""
    dm=[DM(128,0,0)]
    dm.append(DM(128,2000.,5.*60))
    lgs=[]
    for i in range(9):
        lgs.append(LGS(128,2*60.,i*360./9))
    ngs=[]
    for i in range(3):
        ngs.append(LGS(4,1*60.,i*360./3))
    if sparse:
        pmx=sparseTomo(ngsList=ngs,lgsList=lgs,dmList=dm)
    else:
        pmx=tomographicPokeMatrix(ngs,lgs,dm,42.,None)
    return pmx
def testTomoEagleSmall(sparse=0,nact=32):
    """Possibly poke matrix for EAGLE. 40GB in size."""
    dm=[DM(nact,0,0)]
    dm.append(DM(nact,2000.,5.*60))
    lgs=[]
    for i in range(9):
        lgs.append(LGS(nact,2*60.,i*360./9))
    ngs=[]
    for i in range(3):
        ngs.append(LGS(4,1*60.,i*360./3))
    if sparse:
        pmx=sparseTomo(ngsList=ngs,lgsList=lgs,dmList=dm)
    else:
        pmx=tomographicPokeMatrix(ngs,lgs,dm,42.,None)
    return pmx
#To solve for the wavefront you can try:
#import util.createPokeMx,numpy.random,numpy,util.spmatrix,scipy.linsolve,time
#t1=time.time()
#pmx=util.createPokeMx.testTomoEagle(1)
#t2=time.time()
#cents=numpy.random.random(pmx.ncent)
#t3=time.time()
#pTc=pmx.dot(cents)#this is very quick.
#t4=time.time()
#pTp=util.spmatrix.dotWithSelfTransposed(pmx.mx)#this takes 1-2000 seconds, but is once only.
#t5=time.time()
#solution=scipy.linsolve.spsolve(pTp,pTc)#this takes about 150 seconds.
#or
#solution=util.spmatrix.spsolve(pTp,pTc)#gives error warning if LU decomp fails...
#or
#ludecomp=scipy.linsolve.splu(pTp)
#ludecomp.solve(pTc)
#t6=time.time()
#print t2-t1,t3-t2,t4-t3,t5-t4,t6-t5

def getTomoEagle():
    dm=[DM(128,0,0)]
    dm.append(DM(128,2000.,5.*60))
    lgs=[]
    for i in range(9):
        lgs.append(LGS(128,2*60.,i*360./9))
    ngs=[]
    for i in range(3):
        ngs.append(LGS(4,1*60.,i*360./3))
    return ngs, lgs, dm
def test(nact=9):
    """test, including performance.  Current tests show that:
    For 9x9 actuators, sparse is slower (10x).
    For 17x17 act, sparse is slightly faster (2x).
    For 33x33 act, sparse is much faster (10-33x).
    Obviously, coding in c would speed up the sparse implementation."""
    import time
    pmx=straightPokes(nact)
    spmx=sparsePMX(nact,pmx)
    c1=na.random.rand((nact-1)**2*2)#vector
    c2=na.random.rand(((nact-1)**2*2)**2)#matrix
    c2.shape=((nact-1)**2*2,(nact-1)**2*2)
    t1=time.time()
    r=quick.dot(pmx,c1)
    t2=time.time()
    print "numpy vector",t2-t1
    t1=time.time()
    r=quick.dot(pmx,c2)
    t2=time.time()
    print "numpy matrix",t2-t1
    t1=time.time()
    r=spmx.dot(c1)
    t2=time.time()
    print "sparse vector",t2-t1
    t1=time.time()
    r=spmx.dot(c2)
    t2=time.time()
    print "sparse matrix",t2-t1

class tomo:
    def __init__(self,ngsList,lgsList,dmList,telDiam,target):
        self.pmx=tomographicPokeMatrix(ngsList,lgsList,dmList,telDiam,target)
        self.ngsList=ngsList
        self.lgsList=lgsList
        self.dmList=dmList
        self.telDiam=telDiam
        self.target=target
        
class sparseTomo:
    def __init__(self,tpmx=None,ngsList=[],lgsList=[],dmList=[],telDiam=42.):
        """tpmx should be a tomo class object, or None.
        ngsList/lgsList are lists of instances of
        util.guideStar.LGS/NGS and dmList is a list of
        util.dm.physicalDM objects.
        """
        self.ngsList=ngsList
        self.lgsList=lgsList
        self.dmList=dmList
        self.telDiam=telDiam
        self.tpmx=tpmx
        if tpmx!=None:
            self.sparsify(tpmx)
        else:
            self.create()
    def sparsify(self,tpmx):
        """Sparsify a tomographic poke matrix."""
        #Each *diag* of the pmx can be defined by which dm/wfs it
        #corresponds to, and a number representing the position away
        #from the origin to start.  +ve numbers are along x axis, -ve
        #numbers are up y axis.
        #The length of these will be (nact-1)**2
        #The number of these will depend on over/under sampling of wfs.
        #print "TODO - sparseTomo.sparsify()"
        self.mx=scipy.sparse.csc_matrix(tpmx)


    def create(self):
        """Create a tomographic poke matrix stored in sparse format.  Full storage is never required here."""
        ngsList=self.ngsList
        lgsList=self.lgsList
        dmList=self.dmList
        telDiam=self.telDiam
        nacts=0
        ncent=0
        for gs in ngsList+lgsList:
            ncent+=2*(gs.nsubx**2)
        for dm in dmList:
            nacts+=dm.nact**2
        self.nacts=nacts
        self.ncent=ncent
        #pokemxfull=na.zeros((nacts,ncent),na.float32)#create whole pokemx, then strip out unwanted parts later.
        nelem=0
        actcnt=0
        #centList=[]
        centList=scipy.sparse.lil_matrix((nacts,ncent),dtype="f")
        for dm in dmList:
            dmDiam=na.tan(dm.fov/3600./180*na.pi)*dm.height*2+telDiam+1e-6#width of the physical dm. 1micron is for float rounding errors.
            actDiam=dmDiam/(dm.nact-1)
            dmGradient=(1-dm.coupling)/actDiam
            #print dmGradient,dmDiam,actDiam
            centCnt=0
            for gs in ngsList+lgsList:
                subapDiam=telDiam/gs.nsubx
                r=dm.height*na.tan(gs.theta/3600./180.*na.pi)
                wsx=r*na.cos(gs.phi/180.*na.pi)#position of centre of phase on dm.
                wsy=r*na.sin(gs.phi/180.*na.pi)
                wsx+=dmDiam/2.-telDiam/2.#wavefront sensor start coordinates.
                wsy+=dmDiam/2.-telDiam/2.
                #print wsx,wsy,dmDiam,dm.height,gs.theta,gs.phi
                if wsx<0 or wsy<0 or wsx+telDiam>dmDiam or wsy+telDiam>dmDiam:
                    raise Exception("WFS is outside DM fov %g %g %g %g"%(dmDiam,wsx,wsy,dm.fov))
                for i in xrange(dm.nact):
                    ay=i*actDiam#actuator y position
                    wy=ay-wsy#position on wfs - 0 to subapDiam means within first subap etc... subapDiam/2 means the centre of the subap.
                    wymin=wy-actDiam#the possible range of influence of this actuator (will be less if coupling<0.
                    wymax=wy+actDiam
                    if wymax>=-0.5*subapDiam and wymin<=(gs.nsubx+0.5)*subapDiam:
                        # poke may be measurable...
                        for j in xrange(dm.nact):
                            # work out which centroids for this gs will be affected, assume first order only.
                            # is the actuator in the field of view of the wfs?
                            # Assume that the centre of each DM is on-axis.  Note, that actuators are assumed to start at the edge of the DM.
                            #act[i,j] is being poked here of deformable mirror dm.

                            ax=j*actDiam#actuator x position
                            # compute position on wfs...
                            wx=ax-wsx
                            wxmin=wx-actDiam
                            wxmax=wx+actDiam
                            if wxmax>=-0.5*subapDiam and wxmin<=(gs.nsubx+0.5)*subapDiam:
                                #poke will have an affect on the wfs.
                                for sy in range(int(na.floor(wymin/subapDiam)),int(na.ceil(wymax/subapDiam))):
                                    for sx in range(int(na.floor(wxmin/subapDiam)),int(na.ceil(wxmax/subapDiam))):
                                        #the affected centroids...
                                        csy=(sy+0.5)*subapDiam#centre of subap coordinates.
                                        csx=(sx+0.5)*subapDiam
                                        dy=csy-wy
                                        dx=csx-wx
                                        valy=1-dmGradient*na.fabs(dy)
                                        valx=1-dmGradient*na.fabs(dx)
                                        if valy<0:
                                            valy=None
                                        if valx<0:
                                            valx=None
                                        if dy<0 and valy!=None:
                                            valy*=-1
                                        if dx<0 and valx!=None:
                                            valx*=-1
                                        if sy>=0 and sy<gs.nsubx and sx>=0 and sx<gs.nsubx:
                                            actnumber=actcnt+i*dm.nact+j
                                            centnumber=centCnt+sy*gs.nsubx+sx
                                            if valx==None:
                                                valx=0
                                            if valy==None:
                                                valy=0
                                            if valx!=0:
                                                centList[actnumber,centnumber]=valx
                                            if valy!=0:
                                                centList[actnumber,centnumber+ncent/2]=valy
                                            #if valx!=0 or valy!=0:
                                            #    centList.append((actnumber,centnumber,valx,valy))
                            #nelem+=len(centList)
                            #for (sy,sx,valy,valx) in centList:
                            #    if sy>=0 and sy<gs.nsubx and sx>=0 and sx<gs.nsubx:
                            #        actnumber=actcnt+i*dm.nact+j
                            #        centnumber=centCnt+sy*gs.nsubx+sx
                                    #if valx!=None:
                                    #    pokemxfull[actnumber,centnumber]=valx
                                    #if valy!=None:
                                    #    pokemxfull[actnumber,centnumber+ncent/2]=valy
                                    #print i,j,sy,sx,actnumber,centnumber,valx,valy,wxmin,wxmax,wymin,wymax
                #finished with this wfs for now...
                centCnt+=gs.nsubx**2
            #and now we're finished with this DM.
            actcnt+=dm.nact**2
        #centList.sort()
        #centList=na.array(centList).astype(na.float32)
        #self.centList=centList
        self.mx=centList.tocsc()
        #return centList#pokemxfull
                
    def expand(self,coords=None,inarr=None):
        """
        Expands the sparse matrix into a dense (full) matrix.  
        If specified, coords is a tuple of (xmin,ymin,xmax,ymax).  This can allow parts of a huge matrix to be viewed."""
        if type(inarr)==type(None):
            inarr=self.mx
        if type(coords)==type(()) and len(coords)==4:
            xmin,ymin,xmax,ymax=coords
            tmp=scipy.sparse.lil_matrix((ymax-ymin+1,xmax-xmin+1))
            for i in xrange(inarr.data.shape[0]):
                row,col=inarr.rowcol(i)
                if row>=ymin and row<=ymax and col>=xmin and col<=xmax:
                    tmp[row-ymin,col-xmin]=inarr.data[i]
            inarr=tmp.tocsc()
        elif type(coords)!=type(None):
            print "createPokeMx: coords should be a 4-tuple, (xmin,ymin,xmax,ymax) if specified"
        return inarr.todense()
##         if inarr==None:
##             inarr=self.centList
##             ncent=self.ncent
##         else:
##             ncent=(max(inarr[:,1])+1)*(inarr.shape[1]-2)
##         if xmin==None:
##             xmin=0
##         if ymin==None:
##             ymin=0
##         if xmax==None:
##             xmax=(max(inarr[:,1])+1)*(inarr.shape[1]-2)#x and y centroids... or possibly just position if inarr is specified!  The factor inarr.shape should be either 3 or 4.
##         if ymax==None:
##             ymax=max(inarr[:,0])+1#actuators...
##         arr=na.zeros((ymax-ymin,xmax-xmin),na.float32)
##         for i in xrange(inarr.shape[0]):
##             x=inarr[i,1]
##             x2=x+ncent/2
##             y=inarr[i,0]
##             if y>=ymin and y<ymax:
##                 if x>=xmin and x<xmax:
##                     arr[y-ymin,x-xmin]=inarr[i,2]
##                 if inarr.shape[1]==4 and x2>=xmin and x2<xmax:
##                     arr[y-ymin,x2-xmin]=inarr[i,3]
##         return arr

    def transpose(self):
        return self.mx.transpose()
        
    def dot(self,mx,res=None,sparseRes=1):
        """dot product of self with mx (1d or 2d).  Optionally, put the result
        in res... (stops a new array being allocated).
        Result should be the same as numpy.dot(self.expand(),mx)

        If sparseRes, the result will be sparse too.

        If this doesn't work, try self.mx.matvec(vec) for dotting with a vector.
        """
        mx=na.array(mx)

        if type(mx)==na.ndarray and len(mx.shape)==1:#a vector...
            res=self.mx.matvec(mx)
        else:
            tres=self.mx.matmat(mx)
            if sparseRes==0:
                if type(res)!=type(None):
                    res[:,]=tres.todense()
                else:
                    res=tres.todense()
            else:
                res=tres
        return res

    def strip(self,val,mx=None):
        """Sets to zero all elements with abs()<=val."""
        if type(mx)==type(None):
            mx=self.mx
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
##         dim=1
##         nc=self.ncent/2
##         resshape=(self.nacts,)
##         if len(mx.shape)==2:
##             dim=2
##             resshape=(self.nacts,mx.shape[1])
##         if type(res)==type(None):
##             res=na.zeros(resshape,self.centList.dtype.char)
##         else:
##             res[:,]=0
##             if res.shape!=resshape:
##                 raise Exception("result matrix has wrong shape for sparse poke matrix multiplication.  Should be %s, is %s"%(str(resshape),str(res.shape)))
##         if dim==1:
##             self._dot1d(mx,res)
##         elif dim==2:
##             self._dot2d(mx,res)
##         return res

##     def _dot1d(self,mx,res):
##         for i in xrange(self.centList.shape[0]):#for each centroid measurement...
##             x=self.centList[i,1]
##             x2=x+self.ncent/2
##             y=self.centList[i,0]
##             xc=self.centList[i,2]
##             yc=self.centList[i,3]
##             res[y]+=mx[x]*xc+mx[x2]*yc
##     def _dot2d(self,mx,res):
##         for i in xrange(self.centList.shape[0]):#for each centroid measurement...
##             x=self.centList[i,1]
##             x2=x+self.ncent/2
##             y=self.centList[i,0]
##             xc=self.centList[i,2]
##             yc=self.centList[i,3]
##             res[y]+=mx[x]*xc+mx[x2]*yc
##             #for j in xrange(mx.shape[1]):
##             #    res[y,j]+=mx[x,j]*xc+mx[x2,j]*yc
            
                            
##     def dotWithSparseMx(self,spmx,sparse=1,transpose=1,ignore=0.,res=None):
##         """Basically dots self with a sparse matrix.  Result can be sparse or full.
##         If abs(element)<=ignore, the element will be ignored.
##         The format of spmx can either be an array of x,y,val or a sparseTomo object.  The latter case
##         corresponds to a sparse tomographic poke matrix (eg self).  The former is for a general sparse matrix.
##         If supplied, res is the full matrix, not sparse (or a dictionary, which is then converted to sparse).
##         Only an array is returned (not a sparse class object).
##         If transpose is set, the matrix will be transposed before being dotted.  If spmx is a sparseTomo object,
##         this should (usually) be set.
##         """
##         if type(spmx)==na.ndarray:#should be 2D array, with x,y,val, ie shape=(n,3).
##             nc=0
##         else:#is a sparseTomo object.
##             nc=spmx.ncent/2
##             spmx=spmx.centList
##         xindx=1-transpose
##         yindx=transpose

##         if res==None:
##             if sparse:
##                 res={}
##             else:
##                 spmxshape1=max(spmx[:,xindx])+1
##                 res=na.zeros((self.nacts,spmxshape1),na.float32)
##         if type(res)==type({}):
##             restype=1
##         else:
##             restype=0
##             res[:,]=0.
##         for i in xrange(self.centList.shape[0]):
##             for j in xrange(spmx.shape[0]):
##                 val=0.
##                 val2=0.
##                 if self.centList[i,1]==spmx[j,yindx]:
##                     val+=self.centList[i,2]*spmx[i,2]
##                 elif self.centList[i,1]+self.ncent/2==spmx[j,yindx]:
##                     val+=self.centList[i,3]*spmx[i,2]
##                 if nc>0:
##                     if self.centList[i,1]==spmx[j,yindx]:
##                         val2+=self.centList[i,2]*spmx[i,3]
##                     elif self.centList[i,1]+self.ncent/2==spmx[j,yindx]:
##                         val2+=self.centList[i,3]*spmx[i,3]
##                 if abs(val)>ignore:
##                     if restype==1 and not res.has_key((self.centList[i,0],spmx[j,xindx])):
##                         res[self.centList[i,0],spmx[j,xindx]]=0.
##                     res[self.centList[i,0],spmx[j,xindx]]+=val
##                 if abs(val2)>ignore:
##                     if restype==1 and not res.has_key((self.centList[i,0],spmx[j,xindx]+nc)):
##                         res[self.centList[i,0],spmx[j,xindx]+nc]=0.
##                     res[self.centList[i,0],spmx[j,xindx]+nc]+=val2

##         if restype==1:
##             res2=na.zeros((len(res),3),na.float32)
##             i=0
##             for j,k in res.keys():
##                 res2[i]=j,k,res[j,k]
##                 i+=1
##             res=res2
##         return res
