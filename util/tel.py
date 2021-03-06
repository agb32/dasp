#dasp, the Durham Adaptive optics Simulation Platform.
#Copyright (C) 2004-2016 Alastair Basden and Durham University.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU Affero General Public License as
#published by the Free Software Foundation, either version 3 of the
#License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Affero General Public License for more details.

#You should have received a copy of the GNU Affero General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Tel.py : functions and classes to define telescope pupil geometry
"""
#ELT:
#pup=util.tel.Pupil(1600.,800.*0.97,800/3.5*0.97,hexDiam=59.6)
#or for a smaller pupil:
#pup=util.tel.Pupil(1280.,640.*0.97,640/3.5*0.97,hexDiam=59.6*1280/1600.)

#ELT with spiders (39m spiders):
#pup=util.tel.Pupil(1280.,640.*0.97,640/3.5*0.97,hexDiam=59.6*1280/1600.,spider=(6,16.,30.),symmetricHex=0)


#import UserArray
import numpy.lib.user_array as user_array
import numpy
na=numpy
import types
import util.FITS
import util.dist

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

    def __init__(self,npup,r1=None,r2=0,nsubx=None,minarea=0.5,apoFunc=None,nAct=None,dmminarea=None,spider=None,hexDiam=0,hexAreaInner=0.,hexAreaOuter=0.,hexEllipseFact=1.,symmetricHex=0,pupilMap=None,height=0.,fov=0.,outside=0.,rotation=0.):
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
        @type spider: Tuple of (narms, thickness/degrees) or (narms, half-thickness/pixels,offset) or "elt" (not recommended)
        @param hexDiam: If >0, will use hexagons for primary mirror segments.  The diameter (in pixels) of hexagon from point to point (not side to side).
        @type pupType: float.
        @param hexAreaInner: Minarea of vignetted hex mirror segments by obscuration
        @type hexAreaInner: float.
        @param hexAreaOuter: Minarea of vignetted hex mirror segments by primary mirror diameter
        @type hexAreaOuter:  float.
        @param hexEllipseFact: A value to scale the pupil circumscribing circle in the x direction, in case this helps with fitting the hexagons.
        @type hexEllipseFact: float.
        @param symmetricHex: If 1, the resulting pupil function will be symmetric in x and y.  If 0, it might not quite be.
        @type symmetricHex: Int
        @param pupilMap: A filename.fits or 2d array or None.  If not None, will be used to define the pupil function.
        @type pupilMap: None, string, array
        For use with MultiPup object:
          height: The height of this pupil in m.
          fov: fov of this pupil (radius, not diam), in arcsec
          outside: value of pupil outside the fov.
        @type rotation: float, the angle to rotate the pupil.  Only useful for hexagonal pupils obviously, since rotation of a circle is a null operation!
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
        self.height=float(height)
        self.fov=fov
        self.outside=outside
        self.diamIncrease=2*height*numpy.tan(fov*numpy.pi/(180*3600.))
        if r1==None:
            r1=npup/2
        self.r1=r1
        self.r2=r2
        self.spider=spider
        self.nsubx=nsubx
        self.minarea=minarea
        self.hexAreaInner=hexAreaInner
        self.hexAreaOuter=hexAreaOuter
        self.hexDiam=hexDiam
        self.hexEllipseFact=hexEllipseFact
        self.symmetricHex=symmetricHex
        self.apoFunc=apoFunc
        self.rotation=float(rotation)
        if dmminarea==None:
            self.dmminarea=minarea
        else:
            self.dmminarea=dmminarea
        ## we create a grid of x and y lines (to avoid for loops)
        if pupilMap is not None:
            if type(pupilMap)==type(""):
                pupilMap=util.FITS.Read(pupilMap)[1]
            if pupilMap.shape!=(npup,npup):
                raise Exception("pupilMap wrong shape (%s) in tel.py - should be %s"%(str(pupilMap.shape),str((npup,npup))))
            self.fn=pupilMap.astype(numpy.uint8)
            self.area=self.fn.sum()
        elif hexDiam==0:
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
        else:
            self.fn=self.makeHexGridSegmented()
            self.area=self.fn.sum()
##         if type(inarr)!=type(None):
##             self.fn=inarr
        #UserArray.UserArray.__init__(self,self.fn,copy=0)
        #self.shape=self.fn.shape
        if spider=="elt":
            self.makeELTSpider()
        elif spider!=None:
            if len(spider)==2:
                self.makeSpider(spider[0],spider[1])
            elif len(spider)==3:
                armoffset=spider[2]
                narm=spider[0]
                armlist=[((x+armoffset/360.*narm)*360./narm*numpy.pi/180.) for x in range(narm)]
                self.makeSpiderPhysical(armlist,spider[1])
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
        

    def calcHexShape(self,npup,d):               
        r=d/2.
        sx=1.5*r
        sy=r*numpy.sqrt(3)
        #compute number of hexagons across pupil
        nhexx=int(numpy.ceil(npup/sx))
        nhexy=int(numpy.ceil(npup/sy))
        #make an odd number...
        if nhexx%2==0:
            print "Even x"
        if nhexy%2==0:
            print "Even y"
        nhexx+=1-nhexx%2
        nhexy+=1-nhexy%2
        print "nhex y,x:",nhexy,nhexx
        return nhexx,nhexy


    def makeHexGridSegmented(self):
        """Computes a pupil function for a primary made of hexagons.  If these are vignetted by more than hexAreaInner or hexAreaOuter, they are not used.
        Defined from centre of telescope.
        Assumes flat bottomed hexagons.
        Assumes a centred hexagon in the middle.
        Requires the polygon package...
        This way, we divide the circle up into 12 segments and use the distance from centre for these.

        """
        d=self.hexDiam
        r=d/2.
        sx=1.5*r
        sy=r*numpy.sqrt(3)
        #compute number of hexagons across pupil
        nhexx=int(numpy.ceil(self.npup/sx))
        nhexy=int(numpy.ceil(self.npup/sy))
        #make an odd number...
        #nhexx+=1-nhexx%2
        nhexy+=1-nhexy%2
        #print "nhex (for ELT should be 31,36) y,x:",nhexy,nhexx
        pup=numpy.zeros((self.npup+int(d),self.npup+int(d)),numpy.int32)
        #Now, for each hexagon, compute its position, and its vignetting 
        if 1:
            xarr,yarr=numpy.meshgrid(numpy.arange(nhexx)-nhexx//2,numpy.arange(nhexy)-nhexy//2)
            cy=sy*yarr
            cy[:,1::2]+=sy/2.
            cx=sx*xarr
            rcentre=numpy.sqrt(cy*cy+cx*cx)
            theta=numpy.arctan2(cy,cx)
            theta=numpy.where(theta<0,theta+2*numpy.pi,theta)
            seg=((theta+15/180.*numpy.pi)/(30*numpy.pi/180.)).astype(numpy.int32)
            rdist=rcentre*numpy.cos(theta-seg*30*numpy.pi/180.)
            hexmap=numpy.logical_and(rdist<self.r1, rdist>self.r2)
            for i in range(nhexy):
                for j in range(nhexx):
                    if hexmap[i,j]:
                        hexplot,offx,offy=self.drawHex(cx[i,j],cy[i,j],d)
                        offx+=self.npup/2
                        offy+=self.npup/2
                        #print hexplot.shape,offx,offy,j,i,cx[i,j],cy[i,j],rdist[i,j],theta[i,j]-seg[i,j]*30*numpy.pi/180.,seg[i,j],theta[i,j]
                        pup[offy:offy+hexplot.shape[0],offx:offx+hexplot.shape[1]]|=hexplot
        else:#slower
            hexmap=numpy.zeros((nhexy,nhexx),numpy.int32)
            for i in range(nhexy):
                cyreal=sy*(i-nhexy//2)
                for j in range(nhexx):
                    cx=sx*(j-nhexx//2)
                    if j%2==1:#shift y position by half a hex
                        cy=cyreal+sy/2.
                    else:
                        cy=cyreal
                    rcentre=numpy.sqrt(cy*cy+cx*cx)
                    #if rcentre>self.r2*1.3:
                    #    continue
                    theta=numpy.arctan2(cy,cx)
                    if theta<0:
                        theta+=2*numpy.pi
                    seg=int((theta+15/180.*numpy.pi)/(30*numpy.pi/180.))
                    rdist=rcentre*numpy.cos(theta-seg*30*numpy.pi/180.)
                    if rdist<self.r1 and rdist>self.r2:
                        hexmap[i,j]=1
                    if hexmap[i,j]:
                        hexplot,offx,offy=self.drawHex(cx,cy,d)
                        offx+=self.npup/2
                        offy+=self.npup/2
                        #print hexplot.shape,offx,offy,j,i,cx,cy,rdist,theta-seg*30*numpy.pi/180.,seg,theta
                        pup[offy:offy+hexplot.shape[0],offx:offx+hexplot.shape[1]]|=hexplot
                        
        # Now check for filling in of neighbours.
        sd=pup[2:,2:]+pup[:-2,2:]+pup[:-2,:-2]+pup[2:,:-2]
        pup[1:-1,1:-1]|=sd
        if self.symmetricHex:
            pup+=pup[::-1]+pup[:,::-1]+pup[::-1,::-1]
        pup=(pup>0)
        if self.rotation!=0:
            pup=pup.astype("f")
            import cmod.utils
            cmod.utils.rotateArray(pup,self.rotation)
            pup=(pup>=0.3).astype(numpy.uint8)
        #print nhexy,nhexx
        pup=pup[int(d)/2:-int(d)/2,int(d)/2:-int(d)/2]
        #self.puptmp=Pupil(self.npup,self.r1,self.r2).fn.astype("i")
        #self.puptmp+=pup
        return pup.copy()#to make it contiguous.


    def writeZemaxPupil(self,telDiam,fout=None,fullSpider=0):
        """Writes a zemax uda file (user defined aperture)
        Currently, not general - will only write the E-ELT pupil 39m.

        """
        #r is computed from E-TRE-ESO-313-1000_2-63709.pdf
        #using the central obscuration smallest diameter... (and some trig)
        r=9.417/13.*numpy.cos(numpy.arctan(numpy.sqrt(3)/13))
        sx=1.5*r
        sy=r*numpy.sqrt(3)
        d=r*2
        coords=[]
        for i in range(9):
            #rhs.  Centre of these is at 17 hexes across.
            coords.append((17*sx+r,(i-4.5)*sy))
            coords.append((17*sx+r*numpy.sin(30*numpy.pi/180),(i-4)*sy))
        for i in range(4):
            #2nd section...Starts at 17*sr+r, 4.5*sy
            coords.append(((17-i)*sx+r,(4.5+1.5*i)*sy))
            coords.append(((17-i)*sx+r*numpy.sin(30*numpy.pi/180),(5+1.5*i)*sy))
            coords.append(((17-i)*sx+r-r-r*numpy.cos(60*numpy.pi/180),(5+1.5*i)*sy))
            coords.append(((17-i)*sx+r-r-r*numpy.cos(60*numpy.pi/180)-r*numpy.sin(30*numpy.pi/180),(5+1.5*i)*sy+r*numpy.cos(30*numpy.pi/180)))
        i=4
        coords.append(((17-i)*sx+r,(4.5+1.5*i)*sy))
            
        c1=numpy.array(coords)
        c2=numpy.zeros((c1.shape[0]*6,2),c1.dtype)
        rmat=numpy.zeros((2,2),numpy.float64)
        for i in range(6):
            rmat[0,0]=numpy.cos(i*60*numpy.pi/180)
            rmat[0,1]=numpy.sin(i*60*numpy.pi/180)
            rmat[1,0]=-numpy.sin(i*60*numpy.pi/180)
            rmat[1,1]=numpy.cos(i*60*numpy.pi/180)
            c2[i*c1.shape[0]:(i+1)*c1.shape[0]]=numpy.dot(c1,rmat)
        
        #and now the central obs.
        coords=[]
        for i in range(4):
            coords.append((4*sx+r,(i-2)*sy))
            coords.append((4*sx+r*numpy.sin(30*numpy.pi/180),(i-1.5)*sy))
        i=4
        coords.append((4*sx+r,(i-2)*sy))
        c1=numpy.array(coords)
        c3=numpy.zeros((c1.shape[0]*6,2),c1.dtype)
        for i in range(6):
            rmat[0,0]=numpy.cos(i*60*numpy.pi/180)
            rmat[0,1]=numpy.sin(i*60*numpy.pi/180)
            rmat[1,0]=-numpy.sin(i*60*numpy.pi/180)
            rmat[1,1]=numpy.cos(i*60*numpy.pi/180)
            c3[i*c1.shape[0]:(i+1)*c1.shape[0]]=numpy.dot(c1,rmat)
        
        #and the spiders - note, these can't go to zero because of zemax limitations.
        c4=None
        if len(self.spider)==3:
            
            #narm, thickness, offset.
            #Need to define it as triangles (or segments).
            #Do it for one segment, then rotate to get all.
            #On a unit circle first, then scale later.
            narm,thickness,offset=self.spider#actually, half thickness.
            sep=360./narm
            r=telDiam*0.6
            thickness*=telDiam/self.npup
            x0=y0=0.
            x1=r*numpy.cos(offset*numpy.pi/180)
            y1=r*numpy.sin(offset*numpy.pi/180)
            x2=r*numpy.cos((offset+sep)*numpy.pi/180)
            y2=r*numpy.sin((offset+sep)*numpy.pi/180)
            #Now reduce by thickness.
            #Special case for the central one...
            c=numpy.cos
            s=numpy.sin
            a=offset*numpy.pi/180.
            x0,y0=numpy.dot([[c(a),-s(a)],[s(a),c(a)]],[thickness/numpy.tan(sep/2*numpy.pi/180),thickness])
            #and the others...
            x1-=thickness*s(a)
            y1+=thickness*c(a)
            x2+=thickness*c(numpy.pi/2-a-sep*numpy.pi/180)
            y2-=thickness*s(numpy.pi/2-a-sep*numpy.pi/180)
            c1=numpy.array([[x0,y0],[x1,y1],[x2,y2]])
            if fullSpider:
                c1=numpy.concatenate([c1,[[x0,y0]]])
            c4=numpy.zeros((c1.shape[0]*narm,2),c1.dtype)
            for i in range(narm):
                rmat[0,0]=numpy.cos(i*sep*numpy.pi/180)
                rmat[0,1]=numpy.sin(i*sep*numpy.pi/180)
                rmat[1,0]=-numpy.sin(i*sep*numpy.pi/180)
                rmat[1,1]=numpy.cos(i*sep*numpy.pi/180)
                c4[i*c1.shape[0]:(i+1)*c1.shape[0]]=numpy.dot(c1,rmat)
        if fout!=None:
            if fout[-4:]!=".uda":
                fout+=".uda"
            fd=open(fout[:-4]+"pupil.uda","w")
            for i in range(c2.shape[0]):
                fd.write("LIN %g, %g\n"%(c2[i,0],c2[i,1]))
            fd.write("BRK\n")
            for i in range(c3.shape[0]):
                fd.write("LIN %g, %g\n"%(c3[i,0],c3[i,1]))
            fd.write("BRK\n")
            fd.close()
            if c4!=None:
                fd=open(fout[:-4]+"spider.uda","w")
                for i in range(c4.shape[0]):
                    if fullSpider==0 or i%4!=3:
                        fd.write("LIN %g, %g\n"%(c4[i,0],c4[i,1]))
                    if i%(3+fullSpider)==2+fullSpider:
                        fd.write("BRK\n")
                fd.close()
        return c2,c3,c4
            

    def makeHexGridGeometric(self):
        """Computes a pupil function for a primary made of hexagons.  If these are vignetted by more than hexAreaInner or hexAreaOuter, they are not used.
        Defined from centre of telescope.
        Assumes flat bottomed hexagons.
        Assumes a centred hexagon in the middle.
        Requires the polygon package...
        Might be useful, but not for the E-ELT... 
        Not fully finished - data might have some holes in it.
        """
        import Polygon
        import Polygon.Shapes
        circOuter=Polygon.Shapes.Circle(self.r1,(0,0),points=int(2*numpy.pi*self.r1))
        circOuter.scale(self.hexEllipseFact,1.)
        circInner=Polygon.Shapes.Circle(self.r2,(0,0),points=int(2*numpy.pi*self.r2))
        d=self.hexDiam
        hexagon=Polygon.Polygon([(d*numpy.cos(i),d*numpy.sin(i)) for i in numpy.arange(6)*60.*numpy.pi/180.])
        hexarea=hexagon.area()
        r=d/2.
        sx=1.5*r
        sy=r*numpy.sqrt(3)
        #compute number of hexagons across pupil
        nhexx=int(numpy.ceil(self.npup/sx))
        nhexy=int(numpy.ceil(self.npup/sy))
        #make an odd number...
        nhexx+=1-nhexx%2
        nhexy+=1-nhexy%2
        print "nhex (for ELT should be 31,36) y,x:",nhexy,nhexx
        pup=numpy.zeros((self.npup+int(d),self.npup+int(d)),numpy.int32)
        hexmap=numpy.zeros((nhexy,nhexx),numpy.int32)
        #Now, for each hexagon, compute its position, and its vignetting 
        for i in range(nhexy):
            cyreal=sy*(i-nhexy//2)
            for j in range(nhexx):
                cx=sx*(j-nhexx//2)
                if j%2==0:#shift y position by half a hex
                    cy=cyreal+sy/2.
                else:
                    cy=cyreal
                rcentre=numpy.sqrt(cy*cy+cx*cx)
                if rcentre<self.r1*.85:
                    continue
                if rcentre-r<=self.r1:#within the pupil
                    if rcentre+r>=self.r2:#and not wholly within the central obscuration
                        if rcentre-r<self.r2:#partially within central obs
                            #compute overlap
                            hexagon.shift(cx,cy)
                            ol=(hexarea-(hexagon&circInner).area())/hexarea
                            if ol>=self.hexAreaInner:
                                hexmap[i,j]=1
                            hexagon.shift(-cx,-cy)
                        elif rcentre+r>self.r1:#partially outside pupil
                            #compute overlay
                            hexagon.shift(cx,cy)
                            ol=(hexagon&circOuter).area()/hexarea
                            if ol>=self.hexAreaOuter:
                                hexmap[i,j]=1
                            hexagon.shift(-cx,-cy)
                        else:
                            hexmap[i,j]=1
                        if hexmap[i,j]:
                            hexplot,offx,offy=self.drawHex(cx,cy,d)
                            offx+=self.npup/2
                            offy+=self.npup/2
                            print hexplot.shape,offx,offy,j,i,cx,cy
                            pup[offy:offy+hexplot.shape[0],offx:offx+hexplot.shape[1]]+=hexplot
        pup=pup[int(d)/2:-int(d)/2,int(d)/2:-int(d)/2]
        self.puptmp=Pupil(self.npup,self.r1,self.r2).fn.astype("i")
        self.puptmp+=pup
        return pup
    
    def drawHex(self,xpos,ypos,d):
        """Draws a hexagon, returns the array, and the offset of this to place is at x,y."""
        r=d/2.
        sx=1.5*r
        sy=d*numpy.sqrt(3)
        #First, move x,y into the range 0<=x<1.
        x=xpos-numpy.floor(xpos)
        y=ypos-numpy.floor(ypos)
        #how much have we moved?
        xshift=xpos-x
        yshift=ypos-y
        #max size of the hex
        nx=int(numpy.ceil(d+x))
        ny=int(numpy.ceil(d+y))
        #centre of the hex
        xc=x+nx//2
        yc=y+ny//2
        xarr,yarr=numpy.meshgrid(numpy.arange(nx)-xc,numpy.arange(ny)-yc)
        #xarr=xarr-xc
        #yarr=yarr-yc
        theta=numpy.arctan2(yarr,xarr)
        rr=numpy.sqrt(xarr*xarr+yarr*yarr)
        arr=(rr<=numpy.sin(60*numpy.pi/180.)*r/numpy.sin(2*numpy.pi/3.-(theta%(numpy.pi*60/180.)))).astype(numpy.int32)
        #arr=numpy.zeros((ny,nx),numpy.int32)
        #for i in range(ny):
        #    for j in range(nx):
        #        xx=j-xc
        #        yy=i-yc
        #        theta=numpy.arctan2(yy,xx)
        #        rr=numpy.sqrt(xx*xx+yy*yy)
        #        if rr<=numpy.sin(60*numpy.pi/180.)*r/numpy.sin(2*numpy.pi/3.-(theta%(numpy.pi*60/180.))):
        #            arr[i,j]=1
        return arr,xshift,yshift

    
    def drawHexSlow(self,xpos,ypos,d):
        """Draws a hexagon, returns the array, and the offset of this to place is at x,y."""
        r=d/2.
        sx=1.5*r
        sy=d*numpy.sqrt(3)
        #First, move x,y into the range 0<=x<1.
        x=xpos-numpy.floor(xpos)
        y=ypos-numpy.floor(ypos)
        #how much have we moved?
        xshift=xpos-x
        yshift=ypos-y
        #max size of the hex
        nx=int(numpy.ceil(d+x))
        ny=int(numpy.ceil(d+y))
        #centre of the hex
        xc=x+nx//2
        yc=y+ny//2
        arr=numpy.zeros((ny,nx),numpy.int32)
        for i in range(ny):
            for j in range(nx):
                xx=j-xc
                yy=i-yc
                theta=numpy.arctan2(yy,xx)
                rr=numpy.sqrt(xx*xx+yy*yy)
                if rr<=numpy.sin(60*numpy.pi/180.)*r/numpy.sin(2*numpy.pi/3.-(theta%(numpy.pi*60/180.))):
                    arr[i,j]=1
        return arr,xshift,yshift

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
        self.area=self.fn.sum()

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
    
    def makeSpiderPhysical(self,armAngleList,thickness):
        """armAngleList is a list/array of arm angles.  Thickness is in units of npup.  It is actually the radius, not diameter of the beams.
        """
        #xarr,yarr=numpy.meshgrid(numpy.arange(self.npup)-self.npup/2.,numpy.arange(self.npup)-self.npup/2.)
        theta=numpy.fromfunction(lambda y,x:numpy.arctan2(y-self.npup/2+0.5,x-self.npup/2.+0.5),(self.npup,self.npup))
        theta[:]=na.where(theta<0,theta+2*na.pi,theta)
        r=numpy.fromfunction(lambda y,x:numpy.sqrt((x-self.npup/2.+0.5)**2+(y-self.npup/2.+0.5)**2),(self.npup,self.npup))
        arr=numpy.zeros((self.npup,self.npup),numpy.int32)
        for thetaArm in armAngleList:
            arr[:]+=numpy.abs(r*numpy.sin(theta-thetaArm))<thickness
        mask=(arr==0).astype(numpy.int32)
        self.spiderMask=mask
        self.fn*=mask.astype(self.fn.dtype)#different type
        self.area=self.fn.sum()

                             
    def makeELTSpider(self,theta=0.):
        """Attempts to make the appropriate spider corresponding to the ELT pupil design document

        This was send from Clelia.Robert@onera.fr on 29/5/9 to Myers.  The file appears to have been created 16/9/8 according to the footer.
        
        To sum up, there are 4 vanes from centre at 30, 150, 210, 330 degrees.  These are 0.5m wide.  There is also 2 vanes offset at +-4.95m from centre, that span a 141 degree arc.
        The central obscuration is 12.43m diameter.
        Theta is the angle by which the pupil is rotated...
        """
        grid=makeCircularGrid(self.npup)
        #first force the central obscuration to correct size.
        self.fn=na.logical_and((grid<=self.r1),(grid>=self.r1/42.*12.43)).astype("d")
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
        self.area=na.sum(na.sum(self.fn))
        
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
    def rescale(self,npupNew,thresh=0.5):
        """Rescales by spline interpolation to new size."""
        import scipy.interpolate
        x=numpy.arange(self.fn.shape[1])/(self.fn.shape[1]-1.)
        y=numpy.arange(self.fn.shape[0])/(self.fn.shape[0]-1.)
        interpfn=scipy.interpolate.RectBivariateSpline(x,y,self.fn)
        nx=numpy.arange(npupNew)/(npupNew-1.)
        ny=numpy.arange(npupNew)/(npupNew-1.)
        newpup=numpy.array(interpfn(nx,ny)>=thresh).astype(numpy.int32)
        sf=float(npupNew)/self.npup
        pfnNew=Pupil(npupNew,self.r1*sf,self.r2*sf,self.nsubx,self.minarea,pupilMap=newpup)
        return pfnNew
        

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
def smooth(data,degree=5):  
    """Smooths 1D data"""
    window=degree*2-1  
    weight=numpy.ones((window,),numpy.float32)  
    weight=weight/numpy.exp((4*(numpy.arange(window)-degree+1.)/window)**2)
    smoothed=numpy.zeros((len(data)-window,),numpy.float32)
    sw=weight.sum()
    for i in range(len(smoothed)):
        smoothed[i]=(numpy.array(data[i:i+window])*weight).sum()/sw
    return smoothed


class MultiPupil:
    """Class to generate pupil functions with different scalings, and where there are non-ground conjugate parts, so that looking in a different direction gives a different pupil.
    """
    def __init__(self,pupList,diam):
        """pupList: is a list of Pupil objects.
        diam: the diameter of the telescope pupil in m"""
        self.pupList=pupList
        self.diam=diam

    def generatePupil(self,theta,phi,n,tol=0.5,kind="linear"):
        """theta: is the off-axis direction in arcsec.
        phi: is the polar coord in degrees.
        n: the dimensions of the returned pupil.
        """
        import scipy.interpolate
        theta*=numpy.pi/(180*3600.)
        phi*=numpy.pi/180.
        outarr=numpy.ones((n,n),numpy.int8)
        for pup in self.pupList:
            h=pup.height
            diam=pup.diamIncrease+self.diam
            scale=pup.npup/diam#pixels per m.
            w=scale*self.diam#Number of pixels required for the sub-pupil.
            r=h*numpy.tan(theta)
            x=r*numpy.cos(phi)*scale
            y=r*numpy.sin(phi)*scale
            #Select the correct part of this pupil and interpolate
            xstart=x-w/2.+pup.npup/2.
            xoff=xstart%1
            xend=xstart+w
            xstart=int(numpy.floor(xstart))
            xend=int(numpy.ceil(xend))
            ystart=y-w/2.+pup.npup/2.
            yoff=ystart%1
            yend=ystart+w
            ystart=int(numpy.floor(ystart))
            yend=int(numpy.ceil(yend))
            inarr=pup.fn[ystart:yend+1,xstart:xend+1]
            print inarr.shape
            ifn=scipy.interpolate.interp2d(numpy.arange(inarr.shape[1]),numpy.arange(inarr.shape[0]),inarr,kind=kind,copy=False,fill_value=pup.outside)
            #Now interpolate where we want it...
            print xstart,xend,ystart,yend,xoff,yoff
            arr=ifn((numpy.arange(n)+xoff)*(xend-xstart)/float(n),(numpy.arange(n)+yoff)*(yend-ystart)/float(n))
            outarr&=(arr>tol)
        return outarr
            
class CoordPupil(Pupil):

    def setCoords(self,coordsList,telDiam):
        """coordsList is a list of arrays.  Each array contains a boundary inside which is valid pupil.
        If 2 valid pupils overlap, they become invalid. e.g. to do a central obscuration.
        Coords are defined in metres, a 2D array with [N,2], each entry being x,y.
        Last point wraps to first.

        e.g.
        [[(1,1), (5,1), (5,9),(3,2),(1,1)]]

        """
        self.coordsList=coordsList
        self.telDiam=telDiam
    def generatePoly(self,coords):
        from matplotlib.path import Path
        npup=self.npup
        x,y=numpy.meshgrid((numpy.arange(npup)-(npup-1)/2.)/((npup)/2.),(numpy.arange(npup)-(npup-1)/2.)/((npup)/2.))
        x=x.flatten()
        y=y.flatten()
        points=numpy.vstack((x,y)).T
        path=Path(coords/(self.telDiam/2.))
        grid=path.contains_points(points)
        grid.shape=(self.npup,self.npup)
        return grid

    def generate(self):
        fn=0
        for coords in self.coordsList:
            grid=self.generatePoly(coords)
            if fn is None:
                fn=grid.astype(numpy.int8)
            else:
                fn+=grid
        self.fn[:]=fn%2
        self.area=self.fn.sum()
        #compute r2, the inner diameter.
        r2=util.dist.dist(self.fn.shape[0],sqrt=0).ravel()
        nnRad=numpy.argsort(r2) ##we flatten the grid of distances and sort the values
        r2r=numpy.take(r2,nnRad) ## we sort the grid of distances
        xRad=numpy.nonzero(r2r[1:]-r2r[:-1])[0]

        tabNbPointsRad=xRad[1:]-xRad[:-1]

        
        psfSort=numpy.take(self.fn.ravel(),nnRad)
        #### We compute the encircled energy of the PSF
        tabEncircledEnergy=numpy.take(numpy.cumsum(psfSort),xRad)
        #### We compute the radial profile
        tabEnergyPerPixBin=tabEncircledEnergy[1:]-tabEncircledEnergy[:-1]
        #tabEnergyPerPixBin.savespace(1)#prevent conversion to double.
        profil=tabEnergyPerPixBin/tabNbPointsRad
        #pos of first element.
        pos=numpy.argmax(profil)
        self.r2=numpy.sqrt(xRad[pos])/2.


    def makeELTPupil(self,telDiam=38.542):
        """An example to make the E-ELT pupil function.  Data taken from the ESO telescope data package, and also the MOSAIC Baseline Simulation document v0.4 (Gendron et al), E-MOS-TEC-ANR-0015"""
        y=[0,0,614,1229,1229,1843,2458,2458,3071,3686,3686,4297,4913,4913,5523,6138,6144,6144,6759,6759,7368,7983,7990,7990,8605,8605,9212,9826,9834,9834,10448,10448,11053,11666,11675,11675,12287,12287,12890,13500,13510,13510,14120,14129,14129,14739,14748,14748,15356,15365,15365,15972,15980,15980,16587,16594,16594,17199,17206,17206,17810,17815,17815,18418,18423,18423,19023,19028,18431,18431,18434,18434,19034,19036,18438,18438,18439,18439,19038,19038,18439,18439,18438,18438,19036,19034,18434,18434,18431,18431,19028,19023,18423,18423,18418,17815,17815,17810,17206,17206,17199,16594,16594,16587,15980,15980,15972,15365,15365,15356,14748,14748,14739,14129,14129,14120,13510,13510,13500,12890,12287,12287,11675,11675,11666,11053,10448,10448,9834,9834,9826,9212,8605,8605,7990,7990,7983,7368,6759,6759,6144,6144,6138,5523,4913,4913,4297,3686,3686,3071,2458,2458,1843,1229,1229,614,0,0,-614,-1229,-1229,-1843,-2458,-2458,-3071,-3686,-3686,-4297,-4913,-4913,-5523,-6138,-6144,-6144,-6759,-6759,-7368,-7983,-7990,-7990,-8605,-8605,-9212,-9826,-9834,-9834,-10448,-10448,-11053,-11666,-11675,-11675,-12287,-12287,-12890,-13500,-13510,-13510,-14120,-14129,-14129,-14739,-14748,-14748,-15356,-15365,-15365,-15972,-15980,-15980,-16587,-16594,-16594,-17199,-17206,-17206,-17809,-17815,-17815,-18418,-18423,-18423,-19023,-19027,-18431,-18431,-18434,-18434,-19034,-19036,-18438,-18438,-18439,-18439,-19038,-19038,-18439,-18439,-18438,-18438,-19036,-19034,-18434,-18434,-18431,-18431,-19027,-19023,-18423,-18423,-18418,-17815,-17815,-17809,-17206,-17206,-17199,-16594,-16594,-16587,-15980,-15980,-15972,-15365,-15365,-15356,-14748,-14748,-14739,-14129,-14129,-14120,-13510,-13510,-13500,-12890,-12287,-12287,-11675,-11675,-11666,-11053,-10448,-10448,-9834,-9834,-9826,-9212,-8605,-8605,-7990,-7990,-7983,-7368,-6759,-6759,-6144,-6144,-6138,-5523,-4913,-4913,-4297,-3686,-3686,-3071,-2458,-2458,-1843,-1229,-1229,-614]
        x=[18452,18452,18798,18451,18451,18796,18448,18448,18792,18443,18443,18786,18436,18436,18778,18427,17735,17735,17383,17383,17724,17372,16677,16677,16323,16323,16665,16310,15614,15614,15258,15258,15599,15243,14545,14545,14188,14188,14529,14172,13473,13473,13115,12414,12414,12055,11353,11353,10994,10290,10290,9931,9226,9226,8867,8161,8161,7802,7096,7096,6737,6030,6030,5671,4964,4964,4606,3898,3547,3547,2838,2838,2481,1773,1419,1419,710,710,355,-355,-710,-710,-1419,-1419,-1773,-2481,-2838,-2838,-3547,-3547,-3898,-4606,-4964,-4964,-5671,-6030,-6030,-6737,-7096,-7096,-7802,-8161,-8161,-8867,-9226,-9226,-9931,-10290, -10290,-10994,-11353,-11353,-12055,-12414,-12414,-13115,-13473,-13473,-14172, -14529,-14188,-14188,-14545,-14545,-15243,-15599,-15258,-15258,-15613,-15613, -16310,-16664,-16323,-16323,-16677,-16677,-17372,-17724,-17383,-17383,-17735, -17735,-18427,-18778,-18436,-18436,-18786,-18443,-18443,-18792,-18448,-18448, -18796,-18451,-18451,-18798,-18452,-18452,-18798,-18451,-18451,-18796,-18448, -18448,-18792,-18443,-18443,-18786,-18436,-18436,-18778,-18427,-17735,-17735, -17383,-17383,-17724,-17372,-16677,-16677,-16323,-16323,-16664,-16310,-15613, -15613,-15258,-15258,-15599,-15243,-14545,-14545,-14188,-14188,-14529,-14172, -13473,-13473,-13115,-12414,-12414,-12055,-11353,-11353,-10994,-10290,-10290, -9931,-9226,-9226,-8867,-8161,-8161,-7802,-7096,-7096,-6737,-6030,-6030,-5671, -4964,-4964,-4606,-3898,-3547,-3547,-2838,-2838,-2481,-1773,-1419,-1419,-710, -710,-355,355,710,710,1419,1419,1773,2481,2838,2838,3547,3547,3898,4606,4964,4964,5671,6030,6030,6737,7096,7096,7802,8161,8161,8867,9226,9226,9931,10290,10290,10994,11353,11353,12055,12414,12414,13115,13473,13473,14172,14529,14188,14188,14545,14545,15243,15599,15258,15258,15614,15614,16310,16665,16323,16323,16677,16677,17372,17724,17383,17383,17735,17735,18427,18778,18436,18436,18786,18443,18443,18792,18448,18448,18796,18451,18451,18798]

        xinner=[5026,5026,4667,5025,5025,4667,5025,5025,4666,4666,3949,3590,3590,2872,2513,2513,1795,1436,1436,718,359,359,-359,-359,-718,-1436,-1436,-1795,-2513,-2513,-2872,-3590,-3590,-3949,-4666,-4666,-5025,-5025,-4667,-5025,-5025,-4667,-5026,-5026,-4667,-5025,-5025,-4667,-5025,-5025,-4666,-4666,-3949,-3590,-3590,-2872,-2513,-2513,-1795,-1436,-1436,-718,-359,-359,359,359,718,1436,1436,1795,2513,2513,2872,3590,3590,3949,4666,4666,5025,5025,4667,5025,5025,4667]

        yinner=[0,0,622,1243,1243,1865,2487,2487,3108,3108,3109,3730,3730,3731,4352,4352,4353,4974,4974,4974,5595,5595,5595,5595,4974,4974,4974,4353,4352,4352,3731,3730,3730,3109,3108,3108,2487,2487,1865,1243,1243,622,-1,-1,-622,-1243,-1243,-1865,-2487,-2487,-3108,-3108,-3109,-3730,-3730,-3731,-4352,-4352,-4353,-4974,-4974,-4974,-5595,-5595,-5595,-5595,-4974,-4974,-4974,-4353,-4352,-4352,-3731,-3730,-3730,-3109,-3108,-3108,-2487,-2487,-1865,-1243,-1243,-622]
        coordsList=[numpy.array([x,y]).T/1000.,numpy.array([xinner,yinner]).T/1000.]
        self.setCoords(coordsList,telDiam)
        self.generate()
