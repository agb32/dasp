#$Id: tel.py,v 1.32 2010/02/17 06:16:30 ali Exp $
"""
Tel.py : functions and classes to define telescope pupil geometry
"""


#import UserArray
import numpy.lib.user_array as user_array
import numpy
na=numpy
import types

def makeCircularGrid(nxpup,nypup=None,natype=na.float64,dosqrt=1,xoff=0,yoff=0):
    """
    returns a npup*npup numpy array with a grid of distance of pixels
    from the center of the screen

    Distance definition not suitable with FFT !!!

    default return type : Float64
    """
    if nypup==None:
        nypup=nxpup
    tabx=na.arange(nxpup)-float(nxpup/2.)+0.5-xoff ##RWW's convention
    taby=na.arange(nypup)-float(nypup/2.)+0.5-yoff ##RWW's convention
    grid=tabx[na.newaxis,:,]**2+taby[:,na.newaxis]**2
    if dosqrt:
        na.sqrt(grid,grid)
    return grid.astype(natype)

### Anular pupil function #############################################
class Pupil(user_array.container):#UserArray.UserArray):
    """
    Defines telescope pupil geometry

    Class variables (important to simulation programmer):
     - npup : number of pixels of the fn array
     - r1 : radius in PIXELS of the primary mirror
     - r2 : radius in PIXELS of the secondary mirror
     - area : area in PIXELS of the pupil
     - fn : numpy npup*npup array storing the pupil geometry
    @cvar npup: number of pixels of the function array
    @type npup: Int
    @cvar area: Area in pixels of the pupil
    @type area: Int
    @cvar r1: Radius of primary mirror in Pixels
    @type r1: Int
    @cvar r2: Radius of secondary mirror in Pixels
    @type r2: Int
    
    """

    def __init__(self,npup,r1=None,r2=0,nsubx=None,minarea=0.5,apoFunc=None,nAct=None,dmminarea=None,spider=None):
        """ Constructor for the Pupil class

        Parameters: 
         - r1 : radius in PIXELS of the primary mirror
         - r2 : radius in PIXELS of the secondary mirror
         - apoFunc : function defining pupil function in the case of apodised pupils
        @param npup: number of pixels of the function array
        @type npup: Int
        @param apoFunc: Function determining pupil fucntion for apodised pupils
        @type apoFunc: Function
        @param r1: Radius of primary mirror in Pixels
        @type r1: Int
        @param r2: Radius of secondary mirror in Pixels
        @type r2: Int
        @param spider: Definition of spiders.
        @type spider: Tuple of narms, thickness/degrees or "elt"
        """
##         print "creating"
##         inarr=None
##         if type(npup)!=type(1):#assume its an array...
##             inarr=npup
##             npup=inarr.shape[0]
        if nAct!=None:
            raise Exception("nAct in util.tel.Pupil() not used")
        self.npup=npup
        self.checkerboard=None
        self.checkerboardargs=None
        self.area=0.
        if r1==None:
            r1=npup/2
        self.r1=r1
        self.r2=r2
        self.nsubx=nsubx
        self.minarea=minarea
        self.apoFunc=apoFunc
        if dmminarea==None:
            self.dmminarea=minarea
        else:
            self.dmminarea=dmminarea
        ## we create a grid of x and y lines (to avoid for loops)
        grid=makeCircularGrid(npup)

        if type(apoFunc)==types.NoneType:
            self.fn=na.logical_and((grid<=r1),(grid>=r2))
            self.area=na.sum(na.sum(self.fn))
        elif type(apoFunc)==na.ndarray:#ArrayType:
            self.fn=apoFunc*na.logical_and((grid<=r1),(grid>=r2))
            self.area=na.sum(na.sum(self.fn))
        else:
            self.fn=apoFunc(grid)*na.logical_and((grid<=r1),(grid>=r2))
            self.area=na.sum(na.sum(na.logical_and((grid<=r1),(grid>=r2))))
##         if type(inarr)!=type(None):
##             self.fn=inarr
        #UserArray.UserArray.__init__(self,self.fn,copy=0)
        #self.shape=self.fn.shape
        if spider=="elt":
            self.makeELTSpider()
        elif spider!=None:
            self.makeSpider(spider[0],spider[1])
        self.calcSubaps()
        user_array.container.__init__(self,self.fn,copy=0)
    # END of __init__

    def calcSubaps(self):
        """ To be used only from tel.py::Pupil::__init__ """
        self.sum=na.sum(na.sum(self.fn))
        nsubx=self.nsubx
        npup=self.npup
        minarea=self.minarea
        if nsubx!=None:
            #if nAct==None:
            nAct=nsubx+1
            self.nAct=nAct
            self.ndata=0
            self.subflag=na.zeros((nsubx,nsubx),na.int32)
            self.subarea=na.zeros((nsubx,nsubx),na.float64)
            self.dmflag=na.zeros((nAct,nAct),na.int32)
            n=npup/nsubx
            self.pupsub=na.zeros((nsubx,nsubx,n,n),na.float64)
            self.dmpupil=na.zeros((npup,npup),na.float64)
            for i in range(nsubx):
                for j in range(nsubx):
                    self.pupsub[i,j]=self.fn[i*n:(i+1)*n,j*n:(j+1)*n]
                    self.subarea[i,j]=na.sum(na.sum(self.pupsub[i,j]))
                    if self.subarea[i,j]>=minarea*n*n:#flag non-vignetted subaps
                        self.subflag[i,j]=1
                        self.ndata+=2#number of centroids that will be computed (note, 2== 1 for x, 1 for y).
                    if self.subarea[i,j]>self.dmminarea*n*n:#this is only valid for nact==nsubx+1.
                        self.dmflag[i,j]=self.dmflag[i+1,j]=self.dmflag[i,j+1]=self.dmflag[i+1,j+1]=1
                        self.dmpupil[i*n:(i+1)*n,j*n:(j+1)*n]=1.
        
    def getSubapFlag(self,nsubx,minarea=None):
        """Compute the subap flags for a given nsubx"""
        if minarea==None:
            minarea=self.minarea
        subflag=na.zeros((nsubx,nsubx),na.int32)
        n=self.npup/nsubx
        minarea*=n*n
        for i in range(nsubx):
            for j in range(nsubx):
                subarea=na.sum(na.sum(self.fn[i*n:(i+1)*n,j*n:(j+1)*n]))
                if subarea>=minarea:
                    subflag[i,j]=1
        return subflag
        
    def getSubarea(self,nsubx):
        """compute subarea"""
        subarea=na.zeros((nsubx,nsubx),na.float64)
        n=self.npup/nsubx
        #might be faster to use cmod.binimg here?
        for i in xrange(nsubx):
            for j in xrange(nsubx):
                subarea[i,j]=na.sum(self.fn[i*n:(i+1)*n,j*n:(j+1)*n])
        return subarea

    def _rc(self, a):
        if len(na.shape(a)) == 0:
            return a
        else:
            p=self.__class__(self.npup,self.r1,self.r2,self.nsubx,self.minarea,self.apoFunc)
            p.fn=a
            p.array=a
            return p#self.__class__(a)

    def perSubap(self,nsubx=None,vectorAlign=1):
        """Rearrange the data so that memory is sequential for subaps - this is used by the cell...
        If vectorAlign is set, it makes each row of the subap pupils to be 16 bytes in size."""
        if nsubx==None:
            nsubx=self.nsubx
        if nsubx==None:
            raise Exception("Error: Number of x subaps must be defined")
        nphs=self.npup/nsubx
        if vectorAlign:
            nphs_v=(nphs+3)&~3
        else:
            nphs_v=nphs
        a=na.zeros((nsubx,nsubx,nphs,nphs_v),na.float32)
        for i in range(nsubx):
            for j in range(nsubx):
                a[i,j,:,:nphs]=self.fn[i*nphs:(i+1)*nphs,j*nphs:(j+1)*nphs].astype(na.float32)
        return a

    def asDoubleSized(self):
        """return a pupil fn array double the size..."""
        return Pupil(self.npup*2,self.r1*2,self.r2*2,self.nsubx,self.minarea)
    
    def makeSpider(self,narms,thickness):
        """narms is number of arms, thickness is their thickness in degrees.
        Spiders here are purely radial.
        """
        thickness=thickness/180.*na.pi/2#get half-width in radians.
        theta=na.fromfunction(lambda x,y:na.arctan2(x-self.npup/2.+0.5,y-self.npup/2.+0.5),(self.npup,self.npup))
        theta[:]=na.where(theta<0,theta+2*na.pi,theta)
        for k in range(narms):
            tmp=na.abs(theta-k*2.*na.pi/narms)
            arm=na.where((tmp<thickness) + (2*na.pi-tmp<thickness),0,1)
            self.fn*=arm

#        for i in xrange(self.npup):
#            for j in xrange(self.npup):
#                #If pxl i,j falls in a spider, mask it out.
#                theta=na.arctan2(i-self.npup/2.+0.5,j-self.npup/2.+0.5)
#                if theta<0:
#                    theta+=2*na.pi
#                #theta ranges from 0 to 2pi
#                for k in range(narms):
#                    if na.abs(theta-k*2.*na.pi/narms)<thickness or 2*na.pi-na.abs(theta-k*2.*na.pi/narms)<thickness:
#                        self.fn[i,j]=0
            
    def makeELTSpider(self,theta=0.):
        """Attempts to make the appropriate spider corresponding to the ELT pupil design document

        This was send from Clelia.Robert@onera.fr on 29/5/9 to Myers.  The file appears to have been created 16/9/8 according to the footer.
        
        To sum up, there are 4 vanes from centre at 30, 150, 210, 330 degrees.  These are 0.5m wide.  There is also 2 vanes offset at +-4.95m from centre, that span a 141 degree arc.
        The central obscuration is 12.43m diameter.
        Theta is the angle by which the pupil is rotated...
        """
        grid=makeCircularGrid(self.npup)
        #first force the central obscuration to correct size.
        self.fn=na.logical_and((grid<=self.r1),(grid>=self.r1/42.*12.43))
        self.area=na.sum(na.sum(self.fn))
        #now add the spiders...
        #r1 corresponds to 42m in pixels.
        #So 0.5m is r1/42.*0.5
        vt=self.r1/42.*0.5#vane thickness
        ht=vt/2.#vane half-thickness.
        tgrid=na.fromfunction(lambda x,y:na.arctan2(x-self.npup/2.+0.5,y-self.npup/2.+0.5),(self.npup,self.npup))#grid of angles.
        phi=(30+theta)*na.pi/180.#angle of vane
        self.fn*=na.where(na.abs(na.sin(phi-tgrid)*grid)<ht,0,1)
        phi=(-30+theta)*na.pi/180.#angle of vane
        self.fn*=na.where(na.abs(na.sin(phi-tgrid)*grid)<ht,0,1)
        #Now do the off-centre arms...
        xpos=grid*na.sin(theta*na.pi/180.+(na.pi/2-tgrid))#The x coord of a pixel in the rotated (by theta) frame.
        doff=self.r1*2/42.*4.95#distance off centre in pixels
        xoff=-doff*na.sin(theta*na.pi/180)
        yoff=doff*na.cos(theta*na.pi/180)
        grid2=makeCircularGrid(self.npup,xoff=xoff,yoff=yoff)
        tgrid2=na.fromfunction(lambda x,y:na.arctan2(x-self.npup/2.+0.5-yoff,y-self.npup/2.+0.5-xoff),(self.npup,self.npup))#grid of angles.
        #top right spider...
        phi=(19.5+theta)*na.pi/180.#angle of vane
        self.fn*=na.where(na.logical_and(na.abs(na.sin(phi-tgrid2)*grid2)<ht,xpos>=0),0,1)
        #top left spider...
        phi=(-19.5+theta)*na.pi/180.
        self.fn*=na.where(na.logical_and(na.abs(na.sin(phi-tgrid2)*grid2)<ht,xpos<=0),0,1)
        #bottom right spider...
        doff=-self.r1*2/42.*4.95#distance off centre in pixels
        xoff=-doff*na.sin(theta*na.pi/180)
        yoff=doff*na.cos(theta*na.pi/180)
        grid2=makeCircularGrid(self.npup,xoff=xoff,yoff=yoff)
        tgrid2=na.fromfunction(lambda x,y:na.arctan2(x-self.npup/2.+0.5-yoff,y-self.npup/2.+0.5-xoff),(self.npup,self.npup))#grid of angles.
        self.fn*=na.where(na.logical_and(na.abs(na.sin(phi-tgrid2)*grid2)<ht,xpos>=0),0,1)
        #bottom left spider
        phi=(19.5+theta)*na.pi/180.
        self.fn*=na.where(na.logical_and(na.abs(na.sin(phi-tgrid2)*grid2)<ht,xpos<=0),0,1)
        
    def atHeight(self,height,fov,telDiam):
        """Computes a new pupil conjugate at height, with fov.
        height in m
        fov in arcsec (actually the half-fov).
        """
        pxlperm=float(self.r1)/telDiam*2.
        w=2*height*na.tan(fov/3600./180.*na.pi)
        newrad=(telDiam+w)*pxlperm/2.
        npup=self.npup*newrad/self.r1

        r2=self.r2-height*na.tan(fov/3600./180.*na.pi)*pxlperm
        if r2<0.:
            r2=0.
        pup=Pupil(int(npup+0.5),newrad,r2)
        return pup

    def makeCheckerboard(self,sym="hex",n=8,seq=None,rms=1.):
        """Make a checkerboard pattern - eg for pupil segmentation studies...
        Arguments - sym can be"hex" or "square".
        n is the number of checkerboards across the pupil.
        seq is the sequence of patterns to use (eg [-1,1] for sym=="square"
        rms is the RMS of the checkerboard phase.
        """
        if self.checkerboard!=None and self.checkerboardargs==(sym,n,seq,rms):
            return self.checkerboard
        self.checkerboardargs=(sym,n,seq,rms)
        self.checkerboard=None
        arr=numpy.zeros((self.npup,self.npup),numpy.float32)
        if sym=="hex":
            if seq==None:
                seq=[-1,0,1]
            #get the coord of centres...
            #Then assign pixels to nearest to these...
            #increase the sequence length...
            seq=seq*int(numpy.ceil((n*2.)/len(seq)))
            coord=[]
            ysep=self.npup/float(n)
            y=-int(self.npup/float(n)/2.)
            rowno=0
            while y<self.npup+ysep:
                s=seq[:]
                if rowno%2==0:#even row
                    x=-int(self.npup/float(n)/2.)
                    sep=self.npup/float(n)
                    s.append(s.pop(0))
                else:
                    x=0
                    sep=self.npup/float(n)
                #for i in range(rowno%len(seq)):
                #    s.append(s.pop(0))
                while x<self.npup+sep:
                    coord.append([y,x,s.pop(0)])
                    x+=sep
                y+=ysep
                rowno+=1
            coord.reverse()
            coord=numpy.array(coord)
            #now assign the pixels...
            r=range(self.npup)
            for y in r:
                for x in r:
                    dist=(y-coord[:,0])**2+(x-coord[:,1])**2
                    mini=numpy.argmin(dist)
                    arr[y,x]=coord[mini,2]
        elif sym=="square":
            if seq==None:
                seq=[-1,1]
            coord=[]
            seq=seq*int(numpy.ceil((n*2.)/len(seq)))
            sep=self.npup/float(n)
            y=0
            rowno=0
            while y<self.npup+sep:
                s=seq[:]
                if rowno%2==0:
                    s.append(s.pop(0))
                x=0
                while x<self.npup+sep:
                    coord.append([y,x,s.pop(0)])
                    x+=sep
                y+=sep
                rowno+=1
            coord.reverse()
            coord=numpy.array(coord)
            #now assign the pixels...
            r=range(self.npup)
            for y in r:
                for x in r:
                    dist=(y-coord[:,0])**2+(x-coord[:,1])**2
                    mini=numpy.argmin(dist)
                    arr[y,x]=coord[mini,2]
            

        if rms!=None:
            arr*=self.fn
            n=self.fn.sum()
            av=arr.sum()/n
            av2=(arr*arr).sum()/n
            stdev=numpy.sqrt(av2-av*av)
            scale=rms/stdev
            arr*=scale
            

        self.checkerboard=arr
        return arr
                    

class RectangularPupil:
    """Defines telescope pupil geometry, can allow for oval pupil shape.  Note, the pupil is still circular, but the array on which it belongs is rectangular.
    Class variables (important to simulation programmer):
     - nxpup : number of pixels of the fn array in x direction
     - nypup : number of pixels of the fn array in x direction
     - r1 : radius in PIXELS of the primary mirror
     - r2 : radius in PIXELS of the secondary mirror
     - area : area in PIXELS of the pupil
     - fn : numpy npup*npup array storing the pupil geometry
    @cvar nxpup: number of pixels of the function array
    @type nxpup: Int
    @cvar nypup: number of pixels of the function array
    @type nypup: Int
    @cvar area: Area in pixels of the pupil
    @type area: Int
    @cvar r1: Radius of primary mirror in Pixels
    @type r1: Int
    @cvar r2: Radius of secondary mirror in Pixels
    @type r2: Int
    """
    def __init__(self,nxpup,nypup,r1,r2,apoFunc=None):
        """ Constructor for the Pupil class

        Parameters: 
         - r1 : radius in PIXELS of the primary mirror
         - r2 : radius in PIXELS of the secondary mirror
         - apoFunc : function defining pupil function in the case of apodised pupils
        @param nxpup: number of pixels of the function array
        @type nxpup: Int
        @param nypup: number of pixels of the function array
        @type nypup: Int
        @param apoFunc: Function determining pupil fucntion for apodised pupils
        @type apoFunc: Function
        @param r1: Radius of primary mirror in Pixels
        @type r1: Int
        @param r2: Radius of secondary mirror in Pixels
        @type r2: Int
        """
        self.nxpup=nxpup
        self.nypup=nypup
        self.area=0.
        self.r1=r1
        self.r2=r2

        ## we create a grid of x and y lines (to avoid for loops)
        grid=makeCircularGrid(nxpup,nypup)

        if type(apoFunc)==types.NoneType:
            self.fn=na.logical_and((grid<=r1),(grid>=r2))
            self.area=na.sum(na.sum(self.fn))
        else:
            self.fn=apoFunc(grid)*na.logical_and((grid<=r1),(grid>=r2))
            self.area=na.sum(na.sum(na.logical_and((grid<=r1),(grid>=r2))))
    

def rms(arr):
    """compute the RMS of an array..."""
    arr=na.array(arr).flat
    std=na.sqrt(na.average(arr*arr)-na.average(arr)**2)
    return std
