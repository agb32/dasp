import numpy
import threading
from cmod.interp import gslCubSplineInterp,bicubicinterp,linearinterp,gslPeriodicCubSplineInterp
import util.tel,util.dist
import util.FITS
import scipy
import cmod.utils
import util.dot as quick
import time

"""Contains:
dmInfo class - info about a given DM
dmOverview class - overview of all DMs in a system
MirrorSurface class - class for interpolating a surface over actuators.
"""

class dmInfo:
    """A class holding info for a single DM (multiple source directions).
    """
    def __init__(self,label,idlist,height,nact,fov=None,coupling=0.1,minarea=0.25,actoffset=0.,closedLoop=1,actuatorsFrom="reconstructor",primaryTheta=0.,primaryPhi=0.,gainAdjustment=1.,zonalDM=1,actSpacing=None,reconLam=None,subpxlInterp=1,reconstructList="all",pokeSpacing=None,interpType="spline",maxActDist=None,slaving=None,actCoupling=0.,actFlattening=None,alignmentOffset=(0,0),infFunc=None,tiltAngle=0.,tiltTheta=0.,rotation=None,decayFactor=None):
        """idlist is a list of (dm ID,source ID) or just a list of source ID, where dm ID is the idstr for a particular DM object (ie at this height, for a particular direction), and source ID is the idstr for a given source direction.  If this list is just a list of source ID, the dm ID is made by concatenating label with source ID.
        fov is the field of view for this DM, or None in which case the minimum FOV will be computed.
        coupling is the coupling between actuators (fraction), or None.
        gainAdjustment is an adjustment factor which can be multiplied by the global gain to get the gain for this DM.
        if zonalDM==1, is a zonal DM.  Otherwise, is a modal DM, with nact modes.
        reconstructList is used for zonal DMs to specify which actuators should be controlled.  If equal to "all", a standard circular pupil will be assumed.  Othewise, it should be a list of source directions for which actuators are to be controlled, ie typically, a list of the wavefront sensor directions.  Using this can reduce the reconstruction time, and can make a smaller poke matrix.  For an example, compare the dmflag when using and not using this option.
        pokeSpacing is either None, or >0 and <nact, the spacing between actuators that is used when poking more than 1 at once.  It should be large enough so that centroids from a given subap aren't affected by more than one actuator, but small enough to allow the poking to be done quickly.  If equal to zero, a fairly good estimate will be used.

        maxActDist if not None is used instead of minarea.  This is the maximum distance (in units of actspacing) that an actuator is allowed to be from the edge of a dm for it to be used.  A sensible value may be 1.5 for interpolated DMs (sqrt(2)) in openloop, or the influence function width for others.

        slaving can be None (no slaving), "auto" to automatically slave to nearest, or a dictionary of indx:slavelist where indx is the index of the actuator in question, and slavelist is a list of (indx,val) where indx is the index of the slave driver, and val is the fraction of this actuator to use.
        reconLam is in nm.

        alignmentOffset (x,y) is the offset in pixels caused by misalignment.
        infFunc is the influence functions used if interpType=="influence".
        tiltAngle and tiltTheta are used to model a tilted DM.
        rotation is a value or function(that is called every iter and returns the rotation), the rotation angle of the DM in degrees.
        decayFactorAdjustment, if !=None, is used to multiply the output of the previous reconstruction, before adding the new one to it.  A traditional closed loop system would =1, openloop =0, but if wish to do integration with openloop, specify !=0.
        """
        self.label=label#the label for this DM.  This can be used as the same as vdmUser object idstr.
        self.height=height#dm conjugate height.  Zenith is calculated automatically.
        self.nact=nact
        self.fov=fov
        self.coupling=coupling
        self.minarea=minarea
        self.actoffset=actoffset
        self.actSpacing=actSpacing#used if nact == None.
        self.closedLoop=closedLoop
        self.gainAdjustment=gainAdjustment
        self.decayFactor=decayFactor
        self.zonalDM=zonalDM
        self.reconLam=reconLam#the reconstructor wavelength...
        self.subpxlInterp=subpxlInterp
        self.dmpup=None
        self.dmpupil=None
        self.dmflag=None
        self.subarea=None
        self.dmDiam=None
        self.mirrorModes=None
        self.mirrorScale=None
        self.mirrorSurface=None
        self.actCoupling=actCoupling
        self.actFlattening=actFlattening
        self.interpType=interpType
        self.tiltAngle=tiltAngle
        self.tiltTheta=tiltTheta
        self.rotation=rotation
        self.alignmentOffset=alignmentOffset
        self.infFunc=infFunc
        self.maxActDist=maxActDist
        self.slaving=slaving
        self.reconstructList=reconstructList#list of source directions to be reconstructed
        if self.zonalDM==1:# and pokeSpacing!=None and pokeSpacing>0 and pokeSpacing<self.nact:
            self.pokeSpacing=pokeSpacing
        else:
            self.pokeSpacing=None
        if type(actuatorsFrom)==type([]):
            self.actuatorsFrom=actuatorsFrom#list of where the actuator values come from.  EG a list of dm labels (in the case of this being a vdmUser which gets actuators from a number of virtual DMs), or "reconstructor" or an idstr for a reconstructor.
        else:
            self.actuatorsFrom=[actuatorsFrom]
        self.primaryTheta=primaryTheta#the direction the mirror axis points in (only used in MOAO).
        self.primaryPhi=primaryPhi
        if self.primaryTheta!=0 and self.tiltAngle!=0:
            raise Exception("Error - not yet implemented - primaryTheta!=0 and tiltAngle!=0 (are you sure this is neccessary?).")
        self.idlist=[]
        if type(idlist)!=type([]):
            raise Exception("util.dm: idlist must be a list (current %s)"%str(type(idlist)))
        for id in idlist:
            if type(id)==type(""):#this is just the source direction.
                self.idlist.append((self.label+id,id))#create the DM ID.
            else:#should be a tuple of (dm ID, source ID)
                self.idlist.append(id)
                
    def __repr__(self):
        txt="<util.dm.dmInfo instance %s>"%self.label
        return txt
    
    def computeCoords(self,telDiam,fname=None):
        """Creates an array containing coords of each actuator.
        0,0 means on axis.
        You can use gnuplot to plot the relative coordinates of the actuators and sub-aperture centres, e.g. using:
        set style line 2 lt 2#green
        plot "filename.csv" every :1::0::0 using 1:2,"" every :1::1::1 using 1:2 ls 2
        This assumes that there are 2 sets of coordinates in the file (e.g. one DM and one LGS).

        Offset specifies the offset as a fraction of actuator spacing.  0 will put actuators on the edge of subaps
        (if nact==nsubx+1), and 0.5 will put them in the centre (if nact==nsubx).  This helps to specify
        the geometry, eg hudgins, fried etc.
        Note - I'm not sure this has been modified completely to handle tiltAngle!=0.
        """
        if self.zonalDM==0:
            raise Exception("computeCoords failed for a modal DM")
        self.coords=numpy.zeros((self.nact,self.nact,2),numpy.float32)
        offset=self.actoffset
        dmDiam=(numpy.tan(self.fov/3600./180*numpy.pi)*self.height*2+telDiam)/numpy.cos(numpy.pi/180*self.tiltAngle)#width of the physical dm.
        actDiam=dmDiam/(self.nact-1.+2*offset)
        if fname!=None:
            f=open(fname,"a")
            f.write("#DM height %g, fov %g, nact %g, offset %g\n"%(self.height,self.fov,self.nact,offset))
        for i in range(self.nact):
            for j in range(self.nact):
                self.coords[i,j,0]=-dmDiam/2+(j+offset)*actDiam
                self.coords[i,j,1]=-dmDiam/2+(i+offset)*actDiam
                if fname!=None:
                    f.write("%g\t%g\n"%(self.coords[i,j,0],self.coords[i,j,1]))
        if fname!=None:
            f.write("\n")
            f.close()
    def computeEffectiveObscuration(self,npup,telDiam,r2):
        scale=npup/telDiam
        r2=(r2-self.height*numpy.tan(self.fov/3600./180.*numpy.pi)*scale)/numpy.cos(numpy.pi/180*self.tiltAngle)
        if r2<0.:
            r2=0.
        return r2

    def getSlaving(self):
        """returns dict of actuators slaved.  Keys are the actuator index, values are a list of (indx,val) the index of parent actuators and relative strength.
        Or can be None, for no slaving.
        """
        if self.slaving=="auto":
            #need to work out what should be slaved.
            self.slaving={}
            for y in range(self.nact):
                for x in range(self.nact):
                    if self.dmflag[y,x]==0:
                        l1=[]#work out which neighbouring acts are used
                        l2=[]#l1 is nearest, l2 is diagonal nearest.
                        if x>0 and self.dmflag[y,x-1]==1:
                            l1.append(y*self.nact+x-1)
                        if x<self.nact-1 and self.dmflag[y,x+1]==1:
                            l1.append(y*self.nact+x+1)
                        if y>0 and self.dmflag[y-1,x]==1:
                            l1.append((y-1)*self.nact+x)
                        if y<self.nact-1 and self.dmflag[y+1,x]==1:
                            l1.append((y+1)*self.nact+x)
                        if x>0 and y>0 and self.dmflag[y-1,x-1]==1:
                            l2.append((y-1)*self.nact+x-1)
                        if x>0 and y<self.nact-1 and self.dmflag[y+1,x-1]==1:
                            l2.append((y+1)*self.nact+x-1)
                        if x<self.nact-1 and y>0 and self.dmflag[y-1,x+1]==1:
                            l2.append((y-1)*self.nact+x+1)
                        if x<self.nact-1 and y<self.nact-1 and self.dmflag[y+1,x+1]==1:
                            l2.append((y+1)*self.nact+x+1)
                        if len(l1)==0:#no nearest actuators used...
                            l1=l2#use next nearest actuators...
                        if len(l1)>0:
                            l=1./len(l1)
                            self.slaving[y*self.nact+x]=[[xx,l] for xx in l1]
        return self.slaving
    def computeDMPupil(self,atmosGeom,centObscuration=0.,retPupil=1,reconstructList=None):
        """Computes the DM flag and pupil.
        centObscuration is the size of the central obscuration in pixels.  This is reduced here depending on the conjugate height.  Typically, it will equal pupil.r2
        reconstructList is either "all" or a list of sources for which actuators should be computed.  If None, the default is used.
        """
        if reconstructList==None:
            reconstructList=self.reconstructList
        if self.reconstructList==reconstructList:
            #may already have computed stuff...
            if self.dmflag!=None and self.subarea!=None:
                if retPupil and self.dmpupil!=None:
                    return self.dmflag,self.subarea,self.dmpupil
                else:
                    return self.dmflag,self.subarea
        self.reconstructList=reconstructList
        if reconstructList=="all":
            self.dmflag=None
            return self.computeDMPupilAll(atmosGeom,centObscuration,retPupil)
        dmpup=self.calcdmpup(atmosGeom)
        dmpupil=numpy.zeros((dmpup,dmpup),numpy.float32)
        if type(reconstructList)!=type([]):
            reconstructList=[reconstructList]
        xoff=self.height*numpy.tan(self.primaryTheta/60./60./180*numpy.pi)*numpy.cos(self.primaryPhi*numpy.pi/180.)
        yoff=self.height*numpy.tan(self.primaryTheta/60./60./180*numpy.pi)*numpy.sin(self.primaryPhi*numpy.pi/180.)
        pxlscale=atmosGeom.ntel/atmosGeom.telDiam
        if self.maxActDist!=None:
            self.dmflag=numpy.zeros((self.nact,self.nact),numpy.int32)
            #Compute the distance of each actuator from on-axis position.
            #actDist=numpy.zeros((self.nact,self.nact),numpy.float32)
            #actDist[:]=numpy.arange(self.nact)-self.nact/2.+0.5
            #actDist[:]=numpy.sqrt(actDist**2+actDist.transpose()**2)
            #actDist*=self.actSpacing
            #actDist+=numpy.sqrt(xoff*xoff+yoff*yoff)
            actDist=(numpy.arange(self.nact)-self.nact/2.+0.5)*self.actSpacing
        for source in reconstructList:
            #compute the central direction for this source, and then its effective diameter, and place this onto dmpupil at the correct point.
            s=atmosGeom.getSource(source,raiseerror=1)
            theta=s.theta
            phi=s.phi
            alt=s.alt
            telDiam=atmosGeom.telDiam
            height=self.height
            r=height*numpy.tan(theta/3600./180.*numpy.pi)
            sx=(r*numpy.cos(phi/180.*numpy.pi)-xoff)*((1/numpy.cos(numpy.pi/180*self.tiltAngle)-1)*numpy.cos(self.tiltTheta*numpy.pi/180)+1)
            sy=(r*numpy.sin(phi/180.*numpy.pi)-yoff)*((1/numpy.cos(numpy.pi/180*self.tiltAngle)-1)*numpy.sin(self.tiltTheta*numpy.pi/180)+1)
            x=dmpup/2.+sx*pxlscale
            y=dmpup/2.+sy*pxlscale
            #print x,y
            secDiam=centObscuration*telDiam*2/atmosGeom.ntel
            if alt>0:
                telDiam*=numpy.fabs(alt-height)/alt#rescale to lgs cone...
                secDiam*=numpy.fabs(alt-height)/alt
            w=pxlscale*telDiam/2.#width in pixels... (radius)
            diam2=(telDiam*pxlscale)**2/4.
            secDiam2=(secDiam*pxlscale)**2/4.
            if x-w<0 or y-w<0 or x+w>dmpup or y+w>dmpup:
                raise Exception("util.dm.computeDMPupil - bad coords: %g %g %g %g %g"%(x-w,y-w,x+w,y+w,dmpup))
            #now make a circle of this radius around x,y...
            for i in xrange(int(y-w),min(dmpup,int(numpy.ceil(y+w))+1)):
                yi2=int(i-y+0.5)**2
                for j in xrange(int(x-w),min(dmpup,int(numpy.ceil(x+w))+1)):
                    xi2=int(j-x+0.5)**2
                    #print numpy.sqrt([yi2,xi2,diam2,secDiam2])
                    if yi2+xi2<=diam2 and yi2+xi2>=secDiam2:
                        dmpupil[i,j]=1
            if self.maxActDist!=None:
                #for each actuator, find out how far it is from source edge.
                #If it lies within maxActDist*actSpacing of the edge of the mirror then it is allowed.
                #First move into the source coordinates...
                d=numpy.sqrt((actDist[numpy.newaxis]+sx)**2+(actDist[:,numpy.newaxis]+sy)**2)
                minrad=secDiam/2.-self.maxActDist*self.actSpacing
                maxrad=telDiam/2.+self.maxActDist*self.actSpacing
                self.dmflag|=((d<=maxrad).astype(numpy.int32) & (d>=minrad).astype(numpy.int32))
        #now compute the subarea...
        subarea=numpy.zeros((self.nact,self.nact),numpy.float32)
        dmsubapsize=dmpup/(self.nact+2*self.actoffset-1.)
        ypos=dmsubapsize*self.actoffset#position of first actuator...
        for i in xrange(self.nact):
            ys=max(0,int(numpy.floor(ypos-dmsubapsize/2.)))
            ye=min(dmpup,int(numpy.ceil(ypos+dmsubapsize/2.)))
            xpos=dmsubapsize*self.actoffset
            for j in xrange(self.nact):
                xs=max(0,int(numpy.floor(xpos-dmsubapsize/2.)))
                xe=min(dmpup,int(numpy.ceil(xpos+dmsubapsize/2.)))
                subarea[i,j]=numpy.sum(dmpupil[ys:ye,xs:xe])/float((xe-xs)*(ye-ys))
                xpos+=dmsubapsize
            ypos+=dmsubapsize
        if self.maxActDist==None:
            subflag=(subarea>=self.minarea).astype(numpy.int32)
            self.dmflag=subflag
        self.subarea=subarea
        self.dmpupil=dmpupil
        self.nacts=int(self.dmflag.sum())
        if retPupil:
            return self.dmflag,subarea,dmpupil
        else:
            return self.dmflag,subarea
    def computeDMPupilAll(self,atmosGeom,centObscuration=0.,retPupil=1):
        """Computes the DM flag and pupil, for a DM with nAct
        actuators.  actOffset is the offset of the first and last
        actuators into the pupil, ie a value of 0 means actuators are
        at the corners of subaps (if nact==nsubx+1, while a value of
        0.5 means they are in the centre of subaps (if nact==nsubx).

        dmminarea is the minimum area to count in the dmflag.
        centObscuration is the size of the central obscuration in pixels.  This is reduced here depending on the conjugate height.  Typically, it will equal pupil.r2
        """
        if self.zonalDM:
            if retPupil:
                if self.dmflag!=None and self.subarea!=None and self.dmpupil!=None:
                    return self.dmflag,self.subarea,self.dmpupil
            else:
                if self.dmflag!=None and self.subarea!=None:
                    return self.dmflag,self.subarea
            nAct=self.nact
            dmminarea=self.minarea
            actOffset=self.actoffset
            dmflag=numpy.zeros((nAct,nAct),numpy.int32)
            subarea=numpy.zeros((nAct,nAct),numpy.float32)
        else:
            if retPupil:
                if self.dmpupil!=None:
                    return None,None,self.dmpupil
            else:
                return None,None
            dmflag=None
            subarea=None
        dmpup=self.calcdmpup(atmosGeom)
        r1=dmpup/2.
        #reduce central obscuration depending on height - if we go high enough, all the dm will see stuff.
        scale=atmosGeom.npup/atmosGeom.telDiam
        r2=self.computeEffectiveObscuration(atmosGeom.npup,atmosGeom.telDiam,centObscuration)
        print "DM pupil central obscuration %g (radius %g)"%(r2,r1)
        if self.zonalDM:
            if retPupil:
                dmpupil=numpy.zeros((dmpup,dmpup),numpy.float32)
            dmsubapsize=dmpup/(nAct+2*actOffset-1.)
            for i in range(nAct):
                for j in range(nAct):
                    subarea[i,j]=self.integrateArea(nAct,actOffset,dmpup,r1,r2,i,j)
                    #subarea=numpy.sum(numpy.sum(self.fn[i*n:(i+1)*n,j*n:(j+1)*n]))

                    if self.maxActDist==None:
                        if subarea[i,j]>dmminarea*dmsubapsize**2:
                            dmflag[i,j]=1
                            if retPupil:
                                y1=int((i+actOffset-0.5)*dmsubapsize+0.5)
                                y2=int((i+0.5+actOffset)*dmsubapsize+0.5)
                                x1=int((j+actOffset-0.5)*dmsubapsize+0.5)
                                x2=int((j+0.5+actOffset)*dmsubapsize+0.5)
                                if x1<0:
                                    x1=0
                                if y1<0:
                                    y1=0
                                if x2>=dmpup:
                                    x2=dmpup
                                if y2>=dmpup:
                                    y2=dmpup
                                dmpupil[y1:y2,x1:x2]=1
            if self.maxActDist!=None:#have specified maxActDist.
                actDist=(numpy.arange(self.nact)-self.nact/2.+0.5)*self.actSpacing
                d=numpy.sqrt(actDist[numpy.newaxis]**2+actDist[:,numpy.newaxis]**2)
                secRad=r2*atmosGeom.telDiam/atmosGeom.ntel
                minrad=secRad-self.maxActDist*self.actSpacing
                maxrad=r1*atmosGeom.telDiam/atmosGeom.ntel+self.maxActDist*self.actSpacing
                dmflag|=((d<=maxrad).astype(numpy.int32) & (d>=minrad).astype(numpy.int32))
                if retPupil:
                    for i in range(nAct):
                        for j in range(nAct):
                            if dmflag[i,j]:
                                y1=int((i+actOffset-0.5)*dmsubapsize+0.5)
                                y2=int((i+0.5+actOffset)*dmsubapsize+0.5)
                                x1=int((j+actOffset-0.5)*dmsubapsize+0.5)
                                x2=int((j+0.5+actOffset)*dmsubapsize+0.5)
                                if x1<0:
                                    x1=0
                                if y1<0:
                                    y1=0
                                if x2>=dmpup:
                                    x2=dmpup
                                if y2>=dmpup:
                                    y2=dmpup
                                dmpupil[y1:y2,x1:x2]=1
                                
                        
        else:#a modal DM
            if retPupil:
                dmpupil=util.tel.Pupil(dmpup,r1,r2).fn
        self.dmflag=dmflag
        self.subarea=subarea
        if self.zonalDM:
            self.nacts=int(self.dmflag.sum())
        if retPupil:
            self.dmpupil=dmpupil
            return dmflag,subarea,dmpupil
        else:
            return dmflag,subarea        
    def calcdmpup(self,atmosGeom):
        """calculate the dmpup needed for a given dm (conjugate at height).
        If fov is specified (arcsec), the dm will have this fov.  Otherwise it will be just large enough to
        hold all sources.
        Also increased further if the DM is tilted.
        """
        if self.dmpup!=None:
            return self.dmpup
        height=self.height
        fov=self.fov
        #if fov==None:
        #    fov=max(self.atmosGeom.sourceThetas().values())
        npup=atmosGeom.npup
        telDiam=atmosGeom.telDiam
        
        scale=npup/telDiam#pxl per m.
        arcsecRad=numpy.pi/180/3600.

        dmpup=int(numpy.ceil((2*height*numpy.tan(fov*arcsecRad)*scale+npup)/numpy.cos(self.tiltAngle*numpy.pi/180.)))
        print "Calculating dmpup at height %g: %d (fov %g)"%(height,dmpup,fov)
        self.dmpup=dmpup
        return dmpup
    def integrateArea(self,nact,actOffset,dmpup,r1,r2,i,j):
        """Integrates the area of the actuator that will be visible.
        Result is returned in pixels^2."""
        #find the x/y coords of the actuator.
        dmsubapsize=dmpup/(nact+2*actOffset-1.)#size in pixels of inter-actuator distance.
        xc=dmsubapsize*(j+actOffset)
        yc=dmsubapsize*(i+actOffset)
        x1=xc-dmsubapsize/2
        y1=yc-dmsubapsize/2
        x2=xc+dmsubapsize/2
        y2=yc+dmsubapsize/2
        if x1<0:
            x1=0
        if y1<0:
            y1=0
        if x2>=dmpup:
            x2=dmpup
        if y2>=dmpup:
            y2=dmpup
        x1-=dmpup/2.
        x2-=dmpup/2.
        y1-=dmpup/2.
        y2-=dmpup/2.
        #Now, always work in the first quadrant... (simplify things).  If the coords span the axes, may need to do 2 or 4 area calcs to get the total area...
        xcoordList=[]
        ycoordList=[]
        coordList=[]
        if x1<0:
            if x2<=0:
                xcoordList.append((-x2,-x1))
            else:
                xcoordList.append((0,-x1))
                xcoordList.append((0,x2))
        else:
            xcoordList.append((x1,x2))
        if y1<0:
            if y2<=0:
                ycoordList.append((-y2,-y1))
            else:
                ycoordList.append((0,-y1))
                ycoordList.append((0,y2))
        else:
            ycoordList.append((y1,y2))
        for xs in xcoordList:
            for ys in ycoordList:
                if r1>0:
                    coordList.append((xs[0],xs[1],ys[0],ys[1],r1,+1))#add to area if in pupil
                if r2>0:
                    coordList.append((xs[0],xs[1],ys[0],ys[1],r2,-1))#but subtract if also in central obscuration
        #now integrate the quadrants and sum...
        totArea=0.
        for x1,x2,y1,y2,r,sign in coordList:#all x/y >= 0.
            if x1>=r:#x1 was outside the pupil... which means x2 will be as well... so nowt to integrate.
                continue
            if x2>r:
                x2=r
            cy2=self.getCircCoord(x1,r)#bigger than cy1.
            cy1=self.getCircCoord(x2,r)
            area=0
            #Now look to see where we are...
            if cy2<y2:
                if cy2<y1:#outside... (probably never get here)
                    pass
                else:
                    if cy1<y1:
                        #reduce x2 until its on the circle, then integrate.
                        x2=self.getCircCoord(y1,r)
                    #standard integration of circle between x1 and x2 and subtract (x2-x1)*y1 rectangle.
                    area=self.circIntegrate(x1,x2,r)-(x2-x1)*y1
            else:#cy2>=y2
                if cy1>=y2:#just integrate the square...
                    area=(x2-x1)*(y2-y1)
                else:#cy1<y2
                    # get x point where circle entres the square (top)
                    x11=self.getCircCoord(y2,r)
                    if cy1<y1:
                        #get x point where circle leaves the square (bottom)
                        x2=self.getCircCoord(y1,r)
                    # integrate circle from x11 to x2 and add/subtract rectangles.
                    area=(x11-x1)*(y2-y1)-(x2-x11)*y1+self.circIntegrate(x11,x2,r)
            totArea+=area*sign
        if totArea<0:
            print "Warning: Got area<0: %g - setting to zero."%totArea
            totArea=0
        return totArea

    def circIntegrate(self,x1,x2,r):
        """Integrate y=sqrt(r**2-x**2) along x axis.
        """
        t1=numpy.arccos(x1/float(r))
        t2=numpy.arccos(x2/float(r))
        
        area=((t2-numpy.sin(2*t2)/2)-(t1-numpy.sin(2*t1)/2))*r**2/2
        return -area
    def getCircCoord(self,x,r):
        a=r**2-x**2
        if a>=0:
            return numpy.sqrt(a)
        else:
            return None
##     def getMirrorSurface(self):
##         if self.mirrorSurface==None:
##             self.mirrorSurface=MirrorSurface(typ="spline",npup=self.dmpup,nact=self.nact,phsOut=None,actoffset=self.actoffset)
##         return self.mirrorSurface
    def makeMirrorModes(self,atmosGeom,r2,fitpup=1,mirrorSurface=None):
        """
        r2 is pixel radius of central obscuration
        """
        if self.mirrorModes!=None and self.mirrorModes.shape[1:]==(self.dmpup,self.dmpup):
            return self.mirrorModes

        import util.tel
        pupil=util.tel.Pupil(self.dmpup,self.dmpup/2.,self.computeEffectiveObscuration(atmosGeom.npup,atmosGeom.telDiam,r2)).fn
        if self.zonalDM:
            if mirrorSurface==None:
                mirrorSurface=self.getMirrorSurface()
            dmflag=self.computeDMPupil(atmosGeom,centObscuration=r2,retPupil=0)[0]
            nmode=numpy.sum(dmflag)
            dmindices=numpy.nonzero(dmflag.ravel())[0]
            actmap=numpy.zeros((self.nact,self.nact),numpy.float32)
            mirrorModes=numpy.zeros((nmode,self.dmpup,self.dmpup),numpy.float32)
            self.mirrorScale=numpy.zeros((nmode,),numpy.float32)
            #mirrorSurface=MirrorSurface(typ=interpType,npup=self.dmpup,nact=self.nact,phsOut=None,actoffset=self.actoffset,actCoupling=actCoupling,actFlattening=actFlattening)
            for i in xrange(nmode):
                y=dmindices[i]/self.nact
                x=dmindices[i]%self.nact
                actmap[y,x]=1
                mirrorModes[i]=mirrorSurface.fit(actmap)
                if fitpup:
                    mirrorModes[i]*=pupil
                self.mirrorScale[i]=numpy.sqrt(numpy.sum(mirrorModes[i]*mirrorModes[i]))
                mirrorModes[i]/=self.mirrorScale[i]
                actmap[y,x]=0
        else:
            nmode=self.nact
            mirrorModes=util.zernikeMod.Zernike(pupil,nmode,computeInv=0).zern
            print "todo - util.dm.computePhaseCovariance - do we need to scale the zernikes?"
        self.mirrorModes=mirrorModes
        return mirrorModes

    def makeLocalMirrorModes(self,atmosGeom,r2,fitpup=1,mirrorSurface=None,W=None):
        """Makes the mirror modes, but truncated, since they are zonal.
        """
        if not self.zonalDM:
            raise Exception("util.dm - makeLocalMirrorModes only applicable to zonal DMs")
        if W==None:#width as fraction of actuator spacing.  Default is 4, but actually, 8 might be better...
            W=int(4*self.actSpacing/self.dmDiam*self.dmpup+0.5)
            
        dmflag=self.computeDMPupil(atmosGeom,centObscuration=r2,retPupil=0)[0]
        import util.tel
        pupil=util.tel.Pupil(self.dmpup,self.dmpup/2.,self.computeEffectiveObscuration(atmosGeom.npup,atmosGeom.telDiam,r2)).fn
        nmode=numpy.sum(dmflag)
        if mirrorSurface==None:
            mirrorSurface=self.getMirrorSurface()

        if self.mirrorModes!=None and self.mirrorModes.shape==(nmode,W,W):
            return self.mirrorModes,self.mirrorModeCoords,self.vig

        dmindices=numpy.nonzero(dmflag.ravel())[0]
        actmap=numpy.zeros((self.nact,self.nact),numpy.float32)
        mirrorModes=numpy.zeros((nmode,W,W),numpy.float32)
        self.mirrorScale=numpy.zeros((nmode,),numpy.float32)
        mirrorModeCoords=numpy.zeros((nmode,2),numpy.int32)
        vig=numpy.zeros((nmode,),numpy.int32)
        for i in xrange(nmode):
            y=dmindices[i]/self.nact
            x=dmindices[i]%self.nact
            actmap[y,x]=1
            mirrorModeCoords[i,0]=min(self.dmpup-W,max(0,int((x+self.actoffset)/(self.nact-1.+2*self.actoffset)*self.dmpup-W/2.+0.5)))
            mirrorModeCoords[i,1]=min(self.dmpup-W,max(0,int((y+self.actoffset)/(self.nact-1.+2*self.actoffset)*self.dmpup-W/2.+0.5)))
            ymin=mirrorModeCoords[i,1]
            xmin=mirrorModeCoords[i,0]
            ymax=ymin+W
            xmax=xmin+W
            coords=(ymin,xmin,ymax,xmax)
            #print coords
            mirrorSurface.fit(actmap,mirrorModes[i],coords)
            if fitpup:
                mirrorModes[i]*=pupil[ymin:ymax,xmin:xmax]
                if pupil[ymin:ymax,xmin:xmax].sum()!=(ymax-ymin)*(xmax-xmin):
                    vig[i]=1#a vignetted actuator...
            self.mirrorScale[i]=numpy.sqrt(numpy.sum(mirrorModes[i]*mirrorModes[i]))
            mirrorModes[i]/=self.mirrorScale[i]
            actmap[y,x]=0
        self.mirrorModes=mirrorModes
        self.mirrorModeCoords=mirrorModeCoords
        self.vig=vig
        return mirrorModes,mirrorModeCoords,vig
            

    def getMirrorScale(self):
        return self.mirrorScale
        
    def computePhaseCovariance(self,atmosGeom,r2,r0=None,l0=None,lam=None,fitpup=1,nthreads=8,mirrorSurface=None,width=None,dovignetted=1,rescaleModes=0,rescalePhasecov=0):
        """r2 is the diameter of central obscuration of ground dm.
        r0 is frieds paramter
        l0 is outerscale
        lam is wavelength.
        If width is specified, it defines the size of the mirror modes (width x width rather than dmpup x dmpup).  This is good for zonal DMs.  Set to -1 to have defined as 4* actuator spacing (this is probably quite a good value).
        rescalePhasecov should probably be 1...
        """
        if r0!=None or l0!=None:
            print "WARNING util.dm.computePhaseCovariance - r0 and l0 no longer used - taken from atmosGeom."
        if lam==None:
            lam=self.reconLam
        import util.phaseCovariance,util.tel
        if mirrorSurface==None:
            mirrorSurface=self.getMirrorSurface()
        typ="vk"
        if atmosGeom.l0<0:
            typ="kol"
        if width==None:#do the full thing
            pupil=util.tel.Pupil(self.dmpup,self.dmpup/2.,self.computeEffectiveObscuration(atmosGeom.npup,atmosGeom.telDiam,r2))
            mirrorModes=self.makeMirrorModes(atmosGeom,r2,fitpup,mirrorSurface)#self.interpType,actCoupling,actFlattening)
            nmode=mirrorModes.shape[0]
            phasecov=util.phaseCovariance.make(typ=typ,npup=self.dmpup,nmode=nmode,pupil=pupil,modes=mirrorModes,r0=atmosGeom.r0,telDiam=self.dmDiam,l0=atmosGeom.l0,nthreads=nthreads,lam=lam)[3]
#         elif width=="test":
#             mirrorModes=self.makeMirrorModes(atmosGeom,r2,fitpup,mirrorSurface)#self.interpType,actCoupling,actFlattening)
#             phasecov=util.phaseCovariance.computeCov3(mirrorModes,self.nact,self.dmflag,atmosGeom.r0,atmosGeom.l0,atmosGeom.telDiam)
#         elif width=="testlocal":
#             width=None#use 4* actuator spacing by default.
#             mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=width)
#             phasecov=util.phaseCovariance.computeCov3local(mirrorModes,mirrorModeCoords,self.dmpup,self.nact,self.dmflag,atmosGeom.r0,atmosGeom.l0,atmosGeom.telDiam)
#         elif width=="testquick":
#             W=int(4*self.actSpacing/self.dmDiam*self.dmpup+0.5)

#             mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=W)
#             vsum=vig.sum()
#             print "\nVignetted modes: %d/%d\n"%(vsum,mirrorModes.shape[0])
#             if mirrorModes.shape[0]-vsum<self.nact*self.nact:
#                 print "Warning - so many vignetted modes means you may as well to this using the long method"
#             acts=numpy.zeros((self.nact,self.nact),numpy.float32)
#             acts[self.nact/2,self.nact/2]=1
#             phs=mirrorSurface.fit(acts)
#             f=int((self.actoffset+self.nact/2)*self.actSpacing/self.dmDiam*self.dmpup-W/2)
#             mode=phs[f:f+W,f:f+W].copy()
#             scale=numpy.sqrt(numpy.sum(mode*mode))
#             mode/=scale
#             xcoord=numpy.arange(self.nact)*self.actSpacing/self.dmDiam*self.dmpup
#             ycoord=xcoord
#             phasecov=util.phaseCovariance.computeCov4(mode,xcoord,ycoord,self.dmflag,mirrorModes,mirrorModeCoords,vig,self.dmpup,self.nact,atmosGeom.r0,atmosGeom.l0,atmosGeom.telDiam,dovignetted=dovignetted,nThreads=self.interpolationNthreads)
#         elif width=="testquickc":
#             W=int(4*self.actSpacing/self.dmDiam*self.dmpup+0.5)

#             mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=W)
#             vsum=vig.sum()
#             print "\nVignetted modes: %d/%d\n"%(vsum,mirrorModes.shape[0])
#             if mirrorModes.shape[0]-vsum<self.nact*self.nact:
#                 print "Warning - so many vignetted modes means you may as well to this using the long method"
#             acts=numpy.zeros((self.nact,self.nact),numpy.float32)
#             acts[self.nact/2,self.nact/2]=1
#             phs=mirrorSurface.fit(acts)
#             f=int((self.actoffset+self.nact/2)*self.actSpacing/self.dmDiam*self.dmpup-W/2)
#             mode=phs[f:f+W,f:f+W].copy()
#             scale=numpy.sqrt(numpy.sum(mode*mode))
#             mode/=scale
#             xcoord=numpy.arange(self.nact)*self.actSpacing/self.dmDiam*self.dmpup
#             ycoord=xcoord
#             phasecov=util.phaseCovariance.computeCov5(mode,xcoord,ycoord,self.dmflag,mirrorModes,mirrorModeCoords,vig,self.dmpup,self.nact,atmosGeom.r0,atmosGeom.l0,atmosGeom.telDiam,dovignetted=dovignetted,nthreads=nthreads)
#         elif width=="testfft":
#             width=None
#             mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=width)
#             phasecov=util.phaseCovariance.makeWithLocalModesFFT(self.dmpup,mirrorModes,mirrorModeCoords,r0=atmosGeom.r0,l0=atmosGeom.l0,telDiam=self.dmDiam,typ=typ,nthreads=nthreads,lam=lam)
#         elif width=="testfftThreaded":
#             width=None
#             mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=width)
#             phasecov=util.phaseCovariance.makeWithLocalModesFFTThreaded(self.dmpup,mirrorModes,mirrorModeCoords,r0=atmosGeom.r0,l0=atmosGeom.l0,telDiam=self.dmDiam,typ=typ,nthreads=nthreads,lam=lam)
#         elif width=="testfft3":
#             width=None
#             mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=width)
#             phasecov=util.phaseCovariance.makeWithLocalModesFFT3(self.dmpup,mirrorModes,mirrorModeCoords,r0=atmosGeom.r0,l0=atmosGeom.l0,telDiam=self.dmDiam,typ=typ,nthreads=nthreads,lam=lam)

        elif width=="testorig":#do a local mirror mode version.
            width=None
            mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=width)
            phasecov=util.phaseCovariance.makeWithLocalModes(self.dmpup,mirrorModes,mirrorModeCoords,r0=atmosGeom.r0,l0=atmosGeom.l0,telDiam=self.dmDiam,typ=typ,nthreads=nthreads,lam=lam)
        else:#do a local mirror mode version.
            if width<0:
                width=None
            print "Making mirror modes..."
            mirrorModes,mirrorModeCoords,vig=self.makeLocalMirrorModes(atmosGeom,r2,fitpup=fitpup,mirrorSurface=mirrorSurface,W=width)
            if rescaleModes:
                print "rescaling mirrorModes."
                for i in range(mirrorModes.shape[0]):
                    mirrorModes[i]*=self.mirrorScale[i]
            print "Making phase covariance... (mode shape=%s)"%str(mirrorModes.shape)
            phasecov=util.phaseCovariance.makeWithLocalModesFFTThreaded(self.dmpup,mirrorModes,mirrorModeCoords,r0=atmosGeom.r0,l0=atmosGeom.l0,telDiam=self.dmDiam,typ=typ,nthreads=nthreads,lam=lam)
        if rescalePhasecov:
            for i in range(mirrorModes.shape[0]):
                phasecov[i]*=self.mirrorScale[i]
                phasecov[:,i]*=self.mirrorScale[i]
        return phasecov

    def getMirrorSurface(self,interpType=None,actCoupling=None,actFlattening=None,couplingcoeff=0.1,gaussianIndex=2.,gaussianOverlapAccuracy=1e-6,phsOut=None,infFunc=None, interpolationNthreads = 0):
        """Create a MirrorSurface object for this DM"""
        if interpType==None:
            interpType=self.interpType
        if actCoupling==None:
            actCoupling=self.actCoupling
        if actFlattening==None:
            actFlattening=self.actFlattening
        if infFunc==None:
            infFunc=self.infFunc
        return MirrorSurface(typ=interpType,npup=self.dmpup,nact=self.nact,phsOut=phsOut,actoffset=self.actoffset,actCoupling=actCoupling,actFlattening=actFlattening,couplingcoeff=couplingcoeff,gaussianIndex=gaussianIndex,gaussianOverlapAccuracy=gaussianOverlapAccuracy,infFunc=infFunc, interpolationNthreads = interpolationNthreads)

class dmOverview:
    """DM object to hold info about DMs etc.
    Typically, this is used in the param file.
    """
    def __init__(self,dmInfoList,atmosGeom=None):
        """Create an object holding info about all DMs (virtual, physical, etc).  If atmosGeom is specified, this can be used for further
        stuff, eg calculation of number of pixels for each DM etc.
        
        dmInfoList is a list of dmInfo objects.
        atmosGeom is an instance of util.atmos.geom, or None, in which
        case, some functions (if called) will raise error.
        """
        self.dmInfoList=dmInfoList
        self.atmosGeom=atmosGeom


        for dm in dmInfoList:
            if atmosGeom!=None and atmosGeom.zenith!=0:
                print "Reconjugating dm to take zenith into account"
                dm.height/=numpy.cos(atmosGeom.zenith*numpy.pi/180.)
            if dm.fov==None:#atmosGeom must be specified in this case...
                dm.fov=0.
                #need to compute the FOV.  This is done such that all sources will just fit...
                xoff=dm.height*numpy.tan(dm.primaryTheta/60./60./180*numpy.pi)*numpy.cos(dm.primaryPhi*numpy.pi/180.)
                yoff=dm.height*numpy.tan(dm.primaryTheta/60./60./180*numpy.pi)*numpy.sin(dm.primaryPhi*numpy.pi/180.)
                for dmid,sourceid in dm.idlist:
                    x=dm.height*numpy.tan(atmosGeom.sourceTheta(sourceid)*numpy.pi/60./60./180)*numpy.cos(atmosGeom.sourcePhi(sourceid)*numpy.pi/180)-xoff
                    y=dm.height*numpy.tan(atmosGeom.sourceTheta(sourceid)*numpy.pi/60./60./180)*numpy.sin(atmosGeom.sourcePhi(sourceid)*numpy.pi/180)-yoff
                    fov=numpy.arctan2(numpy.sqrt(x**2+y**2),dm.height)/numpy.pi*180*60*60
                    #print sourceid,x,y,fov,dm.fov,dm.height,atmosGeom.sourceTheta(sourceid),atmosGeom.sourcePhi(sourceid),xoff,yoff
                    dm.fov=max(fov,dm.fov)

                print "FOV computed as %g for DM %s"%(dm.fov,dm.label)
            if dm.nact==None:
                #nact must be computed from FOV and actSpacing...
                if dm.zonalDM==0:
                    raise Exception("Modal DM must specify nact as the number of modes")
                if dm.actSpacing==None:
                    raise Exception("actSpacing or nact must be specified")
                dmDiam=(2*dm.height*numpy.tan(dm.fov/3600./180.*numpy.pi)+self.atmosGeom.telDiam)/numpy.cos(numpy.pi/180*dm.tiltAngle)
                dm.nact=int(dmDiam/dm.actSpacing-2*dm.actoffset+1)
                print "nact computed as %g for DM %s"%(dm.nact,dm.label)
                if dm.pokeSpacing!=None and (dm.pokeSpacing<0 or dm.pokeSpacing>=dm.nact):
                    dm.pokeSpacing=None
                #note, due to rounding, actspacing will now be inconsistent with nact.  Maybe it would be better to change actoffset in this case?
                old=dm.actoffset
                dm.actoffset=(dmDiam/dm.actSpacing-(dm.nact-1))/2.
                print "Recomputed actoffset to be %g (was %g)"%(dm.actoffset,old)
            else:
                if dm.actSpacing!=None:
                    print "WARNING: Overriding user defined actSpacing with something that works..."
                dmDiam=(2*dm.height*numpy.tan(dm.fov/3600./180.*numpy.pi)+self.atmosGeom.telDiam)/numpy.cos(numpy.pi/180*dm.tiltAngle)
                dm.actSpacing=dmDiam/(dm.nact+2*dm.actoffset-1.)
            if dm.reconLam==None:#use the wavelength at which the phase is at.
                dm.reconLam=atmosGeom.sourceLambda(dm.idlist[0][1])
                print "Assuming reconstructor wavelength of %s for DM %s"%(dm.reconLam,dm.label)
            dm.dmDiam=(2*dm.height*numpy.tan(dm.fov/3600./180.*numpy.pi)+self.atmosGeom.telDiam)/numpy.cos(numpy.pi/180*dm.tiltAngle)
            dm.calcdmpup(atmosGeom)
        #self.dmDict={}
        #self.dmOrder=[]
        #for dm in self.dmList:
        #    if self.dmDict.has_key(dm[0]):
        #        raise Exception("util.dm: DM idstr specified more than once")
        #    self.dmDict[dm[0]]=dm[1:,]
        #    self.dmOrder.append(dm[0])
        #v=0
        #for id in self.dmOrder:
        #    v+=self.getVirtualFlag(id)
        #if v:#there are virtual DMs.
        #    self.hasVirtual=1
        #else:
        #    self.hasVirtual=0
    def makeNactsList(self,reconID,centObs):
        dmlist=self.makeDMList(reconID)
        nactsList=[]
        nactsCumList=[0]
        for dm in dmlist:
            if dm.zonalDM:
                tmp=dm.computeDMPupil(self.atmosGeom,centObscuration=centObs,retPupil=0)
                nactsList.append(dm.nacts)#int(numpy.sum(tmp[0].ravel())))
            else:#modal.
                nactsList.append(dm.nact)
            nactsCumList.append(nactsList[-1]+nactsCumList[-1])
        return nactsCumList
        
    def getDM(self,idstr,raiseerror=0):
        """Get the DM which has object with DM ID idstr."""
        for dm in self.dmInfoList:
            for id in dm.idlist:
                if id[0]==idstr:
                    return dm
        #now check the dm.label to see if idstr is actually a dm label.
        for dm in self.dmInfoList:
            if dm.label==idstr:
                #print "Warning - found dm label %s to match idstr"%idstr
                return dm
        
        if raiseerror:
            print "DM idstr %s not found in getDM"%idstr
            dm=None
            raise Exception("DM idstr %s not found"%idstr)
        else:
            dm=self.dmInfoList[0]
            print "WARNING - util.dm.dmOverview - DM idstr %s not found in getDM, using DM %s instead"%(idstr,dm.idlist[0][0])
        return dm
    def getDMDiam(self,idstr):
        dm=self.getDM(idstr)
        dmDiam=(2*dm.height*numpy.tan(dm.fov/3600./180.*numpy.pi)+self.atmosGeom.telDiam)/numpy.cos(numpy.pi/180*dm.tiltAngle)
        return dmDiam

    def getDMFromLabel(self,label):
        for dm in self.dmInfoList:
            if dm.label==label:
                return dm
        print "Warning: util.dm - DM with label %s not found, returning None"%str(label)
        return None
    def getHeight(self,idstr):
        dm=self.getDM(idstr)
        if dm==None:
            print "dm idstr %s not found"%idstr
            idstr=self.dmInfoList[0].idlist[0][0]
            dm=self.getDM(idstr)
            print "using %s for idstr"%idstr
        return dm.height
    def getSourceID(self,idstr,raiseerror=0):
        dm=self.getDM(idstr)
        for id in dm.idlist:
            if id[0]==idstr:
                return id[1]
        if raiseerror:
            raise Exception("Source ID for DM id %s not found"%str(idstr))
        else:
            ret=dm.idlist[0][1]
            print "Warning - source ID for DM id %s not found.  Using %s"%(str(idstr),ret)
            return ret
    
    def getnAct(self,idstr):
        dm=self.getDM(idstr)
        if dm.zonalDM==0:
            raise Exception("getnAct called for a modal DM")
        return dm.nact

    def getnModes(self,idstr):
        dm=self.getDM(idstr)
        if dm.zonalDM:
            raise Exception("getnModes called for a zonal DM")
        return dm.nact
    def getfov(self,idstr):
        fov=self.getDM(idstr).fov
        if fov==None:
            print "Taking max fov for dm %s"%idstr
            fov=max(self.atmosGeom.sourceThetas().values())
        return fov
    def getcoupling(self,idstr):
        c=self.getDM(idstr).coupling
        if c==None:
            c=0.1
        return c
    def getminarea(self,idstr):
        if self.maxActDist!=None:
            raise Exception("Call to dm.getminarea when maxActDist not None")
        c=self.getDM(idstr).minarea
        if c==None:
            c=0.25
        return c
    def getactoffset(self,idstr):
        c=self.getDM(idstr).actoffset
        if c==None:
            c=0.
        return c
    def getClosedLoopFlag(self,idstr):
        return self.getDM(idstr).closedLoop
    #def getVirtualFlag(self,idstr):
    #    v=0
    #    if len(self.dmDict[idstr])>=9:
    #        v=self.dmDict[idstr][8]
    #    return v
    def getSourceTheta(self,idstr):
        return self.atmosGeom.sourceTheta(self.getSourceID(idstr))
    def getSourcePhi(self,idstr):
        return self.atmosGeom.sourcePhi(self.getSourceID(idstr))
    def calcdmpup(self,idstr):
        """calculate the dmpup needed for a given dm (conjugate at height).
        If fov is specified (arcsec), the dm will have this fov.  Otherwise it will be just large enough to
        hold all sources.
        """
        height=self.getHeight(idstr)
        fov=self.getfov(idstr)
        #if fov==None:
        #    fov=max(self.atmosGeom.sourceThetas().values())
        npup=self.atmosGeom.npup
        telDiam=self.atmosGeom.telDiam
        
        scale=npup/telDiam#pxl per m.
        arcsecRad=numpy.pi/180/3600.

        dmpup=int(numpy.ceil((2*height*numpy.tan(fov*arcsecRad)*scale+npup)/numpy.cos(numpy.pi/180*self.getDM(idstr).tiltAngle)))
        print "Calculating dmpup at height %g: %d (fov %g)"%(height,dmpup,fov)
        return dmpup
    def makeDMList(self,actsFrom=None):
        """Makes a list of the unique DMs (virtual or physical), for
        which actsFrom is in the dm.actuatorsFrom list or is the dm.label.
        So, actsFrom can be a list containing a combination of
        reconstructor idstr, and DM labels.

        This can then be put into
        createPokeMx stuff.  If actsFrom==None, will scan through all dmInfo
        objects, and use those who have actuatorsFrom==["reconstructor"].
        """
        if actsFrom==None:
            actsFrom=["reconstructor"]
        if type(actsFrom)!=type([]):
            actsFrom=[actsFrom]
        actsFrom=actsFrom[:]#make a copy...
        dml=[]
        if "reconstructor" in actsFrom:
            useRecon=1
            actsFrom.remove("reconstructor")
        else:
            useRecon=0
        for dm in self.dmInfoList:
            if useRecon:
                if "reconstructor" in dm.actuatorsFrom:
                        dml.append(dm)
            for af in actsFrom:
                if af in dm.actuatorsFrom:#af is the idstr of the reconstructor
                    if dm not in dml:
                        dml.append(dm)
                if af==dm.label:#af is the label of a virtual DM
                    if dm not in dml:
                        dml.append(dm)
                #Note, if you are using a vdm, and using actsFrom is the vdm.actuatorsFrom, then dml will include the vdm.  This should be rmeoved manually if necessary.  vdmUser does this automatically.
                
##             for af in dm.actuatorsFrom:
##                 if af in actsFrom:
##                     if dm not in dml:#if its not already in the list?
##                         dml.append(dm)
                    
#                if (dm[1],self.getVirtualFlag(dm[0])) not in height:
#                    idlist.append(dm[0])
#                    height.append((dm[1],self.getVirtualFlag(dm[0])))
#        for id in idlist:
#            if (dmType=="virtual" and self.getVirtualFlag(id)) or dmType=="all" or (dmType=="physical" and not self.getVirtualFlag(id)):
#                dml.append(physicalDM(self.getnAct(id),self.getHeight(id),self.getfov(id),self.getcoupling(id),self.getactoffset(id),self.getVirtualFlag(virtual)))
        return dml

    #def getAnID(self,height,virtual=0):
    #    """Gets one idstr at a conjugate height.  Used by tomoRecon."""
    #    poss=None
    #    for key in self.dmDict.keys():
    #        if self.getHeight(key)==height:
    #            if self.getVirtualFlag(key)==virtual:
    #                return key
    #            else:
    #                poss=key
    #    if poss!=None:
    #        print "WARNING: util.dm - DM at correct height found but not of correct type (virtual/physical)"
    #        return poss
    #    raise Exception("DM at height %g not found"%height)
            

    def computeDMPupil(self,idstr,centObscuration=0.,retPupil=1):
        """Computes the DM flag and pupil, for a DM with nAct
        actuators.  actOffset is the offset of the first and last
        actuators into the pupil, ie a value of 0 means actuators are
        at the corners of subaps (if nact==nsubx+1, while a value of
        0.5 means they are in the centre of subaps (if nact==nsubx).

        dmminarea is the minimum area to count in the dmflag.
        centObscuration is the size of the central obscuration in pixels.  This is reduced here depending on the conjugate height.  Typically, it will equal pupil.r2
        """
        return self.getDM(idstr).computeDMPupil(self.atmosGeom,centObscuration=centObscuration,retPupil=retPupil)
        
##         nAct=self.getnAct(idstr)
##         dmminarea=self.getminarea(idstr)
##         actOffset=self.getactoffset(idstr)
##         dmflag=numpy.zeros((nAct,nAct),numpy.int32)
##         subarea=numpy.zeros((nAct,nAct),numpy.float32)
##         dmpup=self.calcdmpup(idstr)
##         r1=dmpup/2.
##         #reduce central obscuration depending on height - if we go high enough, all the dm will see stuff.
##         scale=self.atmosGeom.npup/self.atmosGeom.telDiam
##         r2=centObscuration-self.getHeight(idstr)*numpy.tan(self.getfov(idstr)/3600./180.*numpy.pi)*scale
##         if r2<0.:
##             r2=0.
##         print "DM pupil central obscuration %g"%r2
##         if retPupil:
##             dmpupil=numpy.zeros((dmpup,dmpup),numpy.float32)
##         dmsubapsize=dmpup/(nAct+2*actOffset-1.)
##         for i in range(nAct):
##             for j in range(nAct):
##                 subarea[i,j]=self.integrateArea(nAct,actOffset,dmpup,r1,r2,i,j)
##                 #subarea=numpy.sum(numpy.sum(self.fn[i*n:(i+1)*n,j*n:(j+1)*n]))
                
##                 if subarea[i,j]>dmminarea*dmsubapsize**2:
##                     dmflag[i,j]=1
##                     if retPupil:
##                         y1=int((i+actOffset-0.5)*dmsubapsize+0.5)
##                         y2=int((i+0.5+actOffset)*dmsubapsize+0.5)
##                         x1=int((j+actOffset-0.5)*dmsubapsize+0.5)
##                         x2=int((j+0.5+actOffset)*dmsubapsize+0.5)
##                         if x1<0:
##                             x1=0
##                         if y1<0:
##                             y1=0
##                         if x2>=dmpup:
##                             x2=dmpup
##                         if y2>=dmpup:
##                             y2=dmpup
##                         dmpupil[y1:y2,x1:x2]=1
##         if retPupil:
##             return dmflag,subarea,dmpupil
##         else:
##             return dmflag,subarea

    def integrateArea(self,nact,actOffset,dmpup,r1,r2,i,j):
        """Integrates the area of the actuator that will be visible.
        Result is returned in pixels^2."""
        #find the x/y coords of the actuator.
        dmsubapsize=dmpup/(nact+2*actOffset-1.)#size in pixels of inter-actuator distance.
        xc=dmsubapsize*(j+actOffset)
        yc=dmsubapsize*(i+actOffset)
        x1=xc-dmsubapsize/2
        y1=yc-dmsubapsize/2
        x2=xc+dmsubapsize/2
        y2=yc+dmsubapsize/2
        if x1<0:
            x1=0
        if y1<0:
            y1=0
        if x2>=dmpup:
            x2=dmpup
        if y2>=dmpup:
            y2=dmpup
        x1-=dmpup/2.
        x2-=dmpup/2.
        y1-=dmpup/2.
        y2-=dmpup/2.
        #Now, always work in the first quadrant... (simplify things).  If the coords span the axes, may need to do 2 or 4 area calcs to get the total area...
        xcoordList=[]
        ycoordList=[]
        coordList=[]
        if x1<0:
            if x2<=0:
                xcoordList.append((-x2,-x1))
            else:
                xcoordList.append((0,-x1))
                xcoordList.append((0,x2))
        else:
            xcoordList.append((x1,x2))
        if y1<0:
            if y2<=0:
                ycoordList.append((-y2,-y1))
            else:
                ycoordList.append((0,-y1))
                ycoordList.append((0,y2))
        else:
            ycoordList.append((y1,y2))
        for xs in xcoordList:
            for ys in ycoordList:
                if r1>0:
                    coordList.append((xs[0],xs[1],ys[0],ys[1],r1,+1))#add to area if in pupil
                if r2>0:
                    coordList.append((xs[0],xs[1],ys[0],ys[1],r2,-1))#but subtract if also in central obscuration
        #now integrate the quadrants and sum...
        totArea=0.
        for x1,x2,y1,y2,r,sign in coordList:#all x/y >= 0.
            if x1>=r:#x1 was outside the pupil... which means x2 will be as well... so nowt to integrate.
                continue
            if x2>r:
                x2=r
            cy2=self.getCircCoord(x1,r)#bigger than cy1.
            cy1=self.getCircCoord(x2,r)
            area=0
            #Now look to see where we are...
            if cy2<y2:
                if cy2<y1:#outside... (probably never get here)
                    pass
                else:
                    if cy1<y1:
                        #reduce x2 until its on the circle, then integrate.
                        x2=self.getCircCoord(y1,r)
                    #standard integration of circle between x1 and x2 and subtract (x2-x1)*y1 rectangle.
                    area=self.circIntegrate(x1,x2,r)-(x2-x1)*y1
            else:#cy2>=y2
                if cy1>=y2:#just integrate the square...
                    area=(x2-x1)*(y2-y1)
                else:#cy1<y2
                    # get x point where circle entres the square (top)
                    x11=self.getCircCoord(y2,r)
                    if cy1<y1:
                        #get x point where circle leaves the square (bottom)
                        x2=self.getCircCoord(y1,r)
                    # integrate circle from x11 to x2 and add/subtract rectangles.
                    area=(x11-x1)*(y2-y1)-(x2-x11)*y1+self.circIntegrate(x11,x2,r)
            totArea+=area*sign
        if totArea<0:
            print "Warning: Got area<0: %g - setting to zero."%totArea
            totArea=0
        return totArea

    def circIntegrate(self,x1,x2,r):
        """Integrate y=sqrt(r**2-x**2) along x axis.
        """
        t1=numpy.arccos(x1/float(r))
        t2=numpy.arccos(x2/float(r))
        
        area=((t2-numpy.sin(2*t2)/2)-(t1-numpy.sin(2*t1)/2))*r**2/2
        return -area
    def getCircCoord(self,x,r):
        a=r**2-x**2
        if a>=0:
            return numpy.sqrt(a)
        else:
            return None
        


## class physicalDM:
##     """A class to hold info about the physical DMs in the system.  Typically used by dm.makeDMList().
##     """
##     def __init__(self,nact,height,fov,coupling=None,offset=0.,virtual=0):
##         """Offset is the offset of actuators relative to the dm surface/subaps.  If 0, actuators are around the edge of the DM, and corners of subaps (for a matched wfs, nact=nsubx+1).  If 0.5, actuators are spaced at the centre of the subaps (nact==nsubx).
##         virtual is a flag as to whether this is a physical or virtual DM.
##         """
##         self.fov=fov#field of view in arcsec (diameter, not radius).  Probably equal to max(GS.theta).
##         self.height=height#conjugate height in m.
##         self.nact=nact#number of actuators.
##         if coupling==None:
##             coupling=0.1
##         self.coupling=coupling#coupling between actuators (rudimentary)
##         self.offset=offset
##         self.virtual=virtual

class MirrorSurface:
    """A class for interpolating actuators onto a mirror surface.
    """
    def __init__(self,typ,npup,nact,phsOut=None,actoffset="fried",actCoupling=None,actFlattening=None,couplingcoeff=0.1,gaussianIndex=2.,gaussianOverlapAccuracy=1e-6,infFunc=None, interpolationNthreads=0):
        """typ can be spline, bicubic, gaussian, linear, influence or pspline.  Others can be added as necessary.
        actCoupling and actFlattening are used for bicubic only.
        actCoupling is also used for spline and pspline
        couplingcoeff, gaussianIndex and gaussianOverlapAccuracy are used for gaussian only.
        infFunc should be used if typ=="influence", and can be an array or a filename.  Shape should be (nact*nact,npup,npup)
        """
        self.typ=typ
        self.npup=npup
        self.nact=nact
        self.infFunc=infFunc
        if phsOut==None:
            self.phsOut=numpy.zeros((self.npup,self.npup),numpy.float32)
        else:
            self.phsOut=phsOut
        if actoffset=="fried":
            actoffset=0.
        elif actoffset=="hudgin":
            actoffset=0.5
        self.actoffset=actoffset
        self.calcCoords(self.npup,self.nact,self.actoffset)
        self.calcPsplineCoords()
        self.actCoupling=actCoupling
        self.actFlattening=actFlattening
        self.couplingcoeff=couplingcoeff
        self.gaussianIndex=gaussianIndex
        self.gaussianOverlapAccuracy=gaussianOverlapAccuracy
        self.interpolationNthreads = interpolationNthreads
        if self.typ=="spline":
            pass
        elif self.typ=="bicubic":
            pass
        elif self.typ=="gaussian":
            self.setupGaussian()
        elif typ=="linear":
            pass
        elif typ=="pspline":#periodic spline...
            pass
        elif typ=="influence":
            self.setupInfluence()

    def fit(self,actmap,phsOut=None,coords=None):
        """coords here can be a tuple of (ymin,xmin,ymax,xmax) over which the data is fitted.
        """
        if self.typ=="spline":
            return self.interpolateSpline(actmap,phsOut,coords=coords)
        elif self.typ=="bicubic":
            return self.interpolateBicubic(actmap,phsOut,coords=coords)
        elif self.typ=="gaussian":
            return self.fitGaussian(actmap,phsOut,coords=coords)
        elif self.typ=="linear":
            return self.interpolateLinear(actmap,phsOut,coords=coords)
        elif self.typ=="pspline":#periodic spline
            return self.interpolatePeriodicSpline(actmap,phsOut,coords=coords)
        elif self.typ=="influence":
            return self.fitInfluence(actmap,phsOut,coords=coords)
        else:
            print "WARNING: mirror surface unknown type - not fitting"

    def calcCoords(self,npup=None,nact=None,actoffset=None):
        if npup==None:
            npup=self.npup
        if nact==None:
            nact=self.nact
        if actoffset==None:
            actoffset=self.actoffset
        spacing=(npup-1.)/(nact-1+actoffset*2.)
        self.x=numpy.arange(npup).astype(numpy.float64)
        self.x2=(numpy.arange(nact)*spacing+spacing*actoffset).astype(numpy.float64)
        return self.x,self.x2
            
    def calcPsplineCoords(self):
        """Calculates coords for periodic spline - this is a special case because of the zero padding used to get the boundary conditions.
        """
        npup=self.npup
        nact=self.nact
        actoffset=self.actoffset

        zpad=10#the zero padding using in gslPeriodicCubSplineInterp.
        #spacing=(npup-1.)/(nact-1+actoffset*2.)
        spacing=1./(nact-1+actoffset*2.)
        self.x2ps=(numpy.arange(nact+2*zpad)*spacing+spacing*actoffset).astype(numpy.float64)
        self.xps=(numpy.arange(npup)+0.5)/npup+zpad/(nact-1+actoffset*2)#self.x2ps[zpad]
        return self.xps,self.x2ps
            
    def interpolateSpline(self,actmap,phsOut=None,coords=None):
        """Interpolation using bicubic spline.
        Geometry can be hudgin (actuators in centre of subaps) or fried (actuators at corners of subaps) or None, in which case actoffset is used.
        If actoffset==0, and nact==nsubx+1, same as fried (actuators in corners of subaps).  If actoffset=0.5 and nact=nsubx, actuators in centre of subaps - same as hudgin.
        """
        actmap=self.fudge(actmap,self.actCoupling)
        if phsOut==None:
            phsOut=self.phsOut
        x2=self.x2
        y=x=self.x
        if coords!=None:
            ymin,xmin,ymax,xmax=coords
            phsOut=phsOut[:ymax-ymin,:xmax-xmin]
            y=y[ymin:ymax]
            x=x[xmin:xmax]
        gslCubSplineInterp(actmap,x2,x2,y,x,phsOut,self.interpolationNthreads)
        return phsOut

    def interpolatePeriodicSpline(self,actmap,phsOut=None,coords=None):
        """Interpolation using bicubic spline.
        Geometry can be hudgin (actuators in centre of subaps) or fried (actuators at corners of subaps) or None, in which case actoffset is used.
        If actoffset==0, and nact==nsubx+1, same as fried (actuators in corners of subaps).  If actoffset=0.5 and nact=nsubx, actuators in centre of subaps - same as hudgin.
        """
        actmap=self.fudge(actmap,self.actCoupling)
        if phsOut==None:
            phsOut=self.phsOut
        x2=self.x2ps
        y=x=self.xps
        if coords!=None:
            ymin,xmin,ymax,xmax=coords
            phsOut=phsOut[:ymax-ymin,:xmax-xmin]
            y=y[ymin:ymax]
            x=x[xmin:xmax]
        gslPeriodicCubSplineInterp(actmap,x2,x2,y,x,phsOut)
        return phsOut

    def interpolateLinear(self,actmap,phsOut=None,coords=None):
        """Interpolation using linear.
        If actoffset==0, and nact==nsubx+1, same as fried (actuators in corners of subaps).  If actoffset=0.5 and nact=nsubx, actuators in centre of subaps - same as hudgin.
        """
        if phsOut==None:
            phsOut=self.phsOut
        x2=self.x2.astype("f")
        y=x=self.x.astype("f")
        if coords!=None:
            ymin,xmin,ymax,xmax=coords
            phsOut=phsOut[:ymax-ymin,:xmax-xmin]
            y=y[ymin:ymax]
            x=x[xmin:xmax]
        linearinterp(actmap.astype("f"),x2,x2,y,x,phsOut)
        return phsOut

        # Replace this with what ever interpolation routine you want to use
        # The input is an nact*nact array corresponding to the actuator values (in radians of phase?)
        # The output is a npup*npup array corresponding to the phase change after reflection from the DM
##         if type(phs_out)==type(None):
##             phs_out  = numpy.zeros((npup,npup),numpy.float32)
##         # May want to oversize interpolated map by a few pixels then truncate to npup*npup to match Fried geometry
##         phs_in = actmap

##         if geom=="fried":
##             step=(actmap.shape[0]-1)/float(npup-1)
##             x2=numpy.array(range(actmap.shape[0]),numpy.float64)
##             x=(numpy.arange(npup)*step).astype(numpy.float64)
##         elif geom=="hudgin":
##             step=actmap.shape[0]/float(npup-1)
##             x2=numpy.array(range(actmap.shape[0]),numpy.float64)+0.5
##             x=(numpy.arange(npup)*step).astype(numpy.float64)
##         else: #use actoffset.
##             x=numpy.arange(phs_out.shape[0])
##             spacing=(phs_out.shape[0]-1)/(actmap.shape[0]-1+actoffset*2)
##             x2=numpy.arange(actmap.shape[0])*spacing
##             x2+=spacing*actoffset#This will range from 0 to nact (not nact-1).



    def interpolateBicubic(self,actmap,phsOut=None,coords=None):
        """A second interpolation routine using bicubic
        interpolation.  This method also slaves adjacent
        actuatory to each other slightly, causing a rise or
        fall depending on what their neighbours are doing - it
        is hoped that this is a better model for the DM.
        """
        actCoupling=self.actCoupling
        actFlattening=self.actFlattening
        if type(phsOut)==type(None):
            phsOut  = self.phsOut#numpy.zeros((npup,npup),numpy.float32)
        actmap2=self.fudge(actmap,fiddle=actCoupling)#adjust actuators slightly.
        actmap2=self.fudge(actmap2,fiddle=actCoupling)
        dx,dy,dxy=self.gradients(actmap2,flatten=actFlattening)#compute gradients.
        if coords!=None:
            print "WARNING: util.dm - interpolateBicubic coords parameter not yet implemented"
        bicubicinterp(actmap2,dy,dx,dxy,phsOut)
        return phsOut

    def setupInfluence(self):
        if self.infFunc==None:
            raise Exception("Influence functions not specified")
        elif type(self.infFunc)==type(""):
            self.infFunc=util.FITS.Read(self.infFunc)[1]
        if self.infFunc.shape!=(self.nact*self.nact,self.npup,self.npup):
            raise Exception("Influence functions should be 3D (nact*nact,dmpup,dmpup")

    def fitInfluence(self,actmap,phsOut=None,coords=None):
        if type(phsOut)==type(None):
            phsOut=self.phsOut
        phsOut[:]=0
        am=actmap.ravel()
        for i in range(am.shape[0]):#for each actuator value...
            phsOut+=am[i]*self.infFunc[i]
        return phsOut

    def setupGaussian(self):
        self.influenceDict={}#keys are the offset in fraction of pixel of the gaussian peak from the nearest mirror pixel.
        self.influenceKey=numpy.zeros((self.nact,self.nact,2),numpy.float32)
        self.infCoords=numpy.zeros((self.nact,self.nact,4),numpy.int32)#coords of the influence function to use (it may be clipped if near the edge
        self.dmCoords=numpy.zeros((self.nact,self.nact,4),numpy.int32)#coords of the DM into which to place the influence function.
        npup=self.npup
        couplingcoeff=self.couplingcoeff
        gaussianIndex=self.gaussianIndex
        lnw=numpy.log(couplingcoeff)
        accuracy=self.gaussianOverlapAccuracy
        actspacing=self.x2[1]-self.x2[0]#in terms of phsout pixels.
        actrange=actspacing*(numpy.log(accuracy)/lnw)**(1./gaussianIndex)#the distance (in pixels) over which the actuator has an effect
        actrange=int(numpy.ceil(actrange))
        print "Taking actuator range to be %d pixels"%actrange
        actrange*=2#in each direction...
        
        for i in xrange(self.nact):
            y=self.x2[i]#get the y coord in terms of dm pixels.
            for j in xrange(self.nact):
                x=self.x2[j]#get the x coord
                self.influenceKey[i,j]=(y%1,x%1)
                key=(self.influenceKey[i,j,0],self.influenceKey[i,j,1])
                if not self.influenceDict.has_key(key):
                    #generate the influence function for this pixel offset
                    dists=util.dist.dist(actrange+1,dy=y%1,dx=x%1)/actspacing
                    self.influenceDict[key]=numpy.exp(lnw*dists**gaussianIndex)
                    actvaltmp=self.influenceDict[key]
                ys=int(y)-actrange/2
                ye=ys+actrange+1
                xs=int(x)-actrange/2
                xe=xs+actrange+1
                if ys<0:
                    self.infCoords[i,j,0]=-ys
                    self.dmCoords[i,j,0]=0
                else:
                    self.infCoords[i,j,0]=0
                    self.dmCoords[i,j,0]=ys
                if ye>npup:
                    self.infCoords[i,j,1]=actrange+1-(ye-npup)
                    self.dmCoords[i,j,1]=npup
                else:
                    self.infCoords[i,j,1]=actrange+1
                    self.dmCoords[i,j,1]=ye
                if xs<0:
                    self.infCoords[i,j,2]=-xs
                    self.dmCoords[i,j,2]=0
                else:
                    self.infCoords[i,j,2]=0
                    self.dmCoords[i,j,2]=xs
                if xe>npup:
                    self.infCoords[i,j,3]=actrange+1-(xe-npup)
                    self.dmCoords[i,j,3]=npup
                else:
                    self.infCoords[i,j,3]=actrange+1
                    self.dmCoords[i,j,3]=xe

    
    def fitGaussian(self,actmap,phsOut=None,coords=None):
        if type(phsOut)==type(None):
            phsOut=self.phsOut
        phsOut*=0
        dm=self.dmCoords
        inf=self.infCoords
        id=self.influenceDict
        ifk=self.influenceKey
        if coords==None:
            for i in xrange(self.nact):
                for j in xrange(self.nact):
                    key=(ifk[i,j,0],ifk[i,j,1])
                    phsOut[dm[i,j,0]:dm[i,j,1],dm[i,j,2]:dm[i,j,3]]+=id[key][inf[i,j,0]:inf[i,j,1],inf[i,j,2]:inf[i,j,3]]*actmap[i,j]
        else:
            print "WARNING util.dm - fitGaussian coords not yet implemeted"
        return phsOut


    def fudge(self,actmap,fiddle=0.1):
        """Alter actuator values depending on nearest neighbours."""
        if fiddle==0 or fiddle==None:
            return actmap
        nact=actmap.shape[0]
        diff=numpy.zeros(actmap.shape,actmap.dtype)
        diff2=diff.copy()
        diff2[1:]=actmap[:-1]
        diff2[:-1]+=actmap[1:]
        diff2[:,1:]+=actmap[:,:-1]
        diff2[:,:-1]+=actmap[:,1:]
        diff2[1:-1,1:-1]/=4
        diff2[1:-1,::diff2.shape[1]-1]/=3
        diff2[::diff2.shape[0]-1,1:-1]/=3
        diff2[::diff2.shape[0]-1,::diff2.shape[1]-1]/=2
#        for i in range(nact):
#            for j in range(nact):
#                nneighbour=0
#                sumneighbour=0.
#                if i>0:
#                    sumneighbour+=actmap[i-1,j]
#                    nneighbour+=1
#                if j>0:
#                    sumneighbour+=actmap[i,j-1]
#                    nneighbour+=1
#                if i<nact-1:
#                    sumneighbour+=actmap[i+1,j]
#                    nneighbour+=1
#                if j<nact-1:
#                    sumneighbour+=actmap[i,j+1]
#                    nneighbour+=1
#                #diff[i,j]=fiddle*(sumneighbour-4*actmap[i,j])/nneighbour
#                diff[i,j]=fiddle*(sumneighbour/nneighbour-actmap[i,j])
#        return actmap+diff
        res=actmap+(diff2-actmap)*fiddle
        if res.dtype!=actmap.dtype:
            res=res.astype(actmap.dtype)
        return res

    def gradients(self,actmap,flatten=1.):
        """Compute the x, y and xy gradients of actmap.
        """
        if flatten==None:
            flatten=1
        x=numpy.zeros(actmap.shape,actmap.dtype)
        y=numpy.zeros(actmap.shape,actmap.dtype)
        xy=numpy.zeros(actmap.shape,actmap.dtype)
        nact=actmap.shape[0]
        for i in range(nact):
            for j in range(nact):
                xm=j-1
                xp=j+1
                ym=i-1
                yp=i+1
                if xm<0:xm=0
                if ym<0:ym=0
                if xp>=nact:xp=nact-1
                if yp>=nact:yp=nact-1
                xy[i,j]=actmap[ym,xm]+actmap[yp,xp]-actmap[yp,xm]-actmap[ym,xp]/(xp-xm)*(yp-ym)
                y[i,j]=(actmap[yp,j]-actmap[ym,j])/(yp-ym)
                x[i,j]=(actmap[i,xp]-actmap[i,xm])/(xp-xm)
        flatten=numpy.array(flatten,actmap.dtype)
        x*=flatten
        y*=flatten
        xy*=flatten
        return x,y,xy

    def fitToWavefront(self,phase,subflag=None,pupmask=None,step=0.01):
        """Attempts to produce a least squares fit of the mirror to phase.
        If subflag is specified, only uses these actuators...
        pupmask can be specified for the area over which the RMS is computed.
        step is the initial step size.
        """
        #First calculate the actuator values from the value of phase at the
        #points of the actuators.
        actmap=numpy.zeros((self.nact,self.nact),numpy.float32)
        actsize=self.npup/(self.nact+2*self.actspacing-1)#pxl per actuator
        ypos=actsize*self.actspacing
        for i in xrange(self.nact):
            xpos=actsize*self.actspacing
            for j in xrange(self.nact):
                actmap[i,j]=phase[min(self.npup,max(0,int(numpy.round(ypos)))),min(self.npup,max(0,int(numpy.round(xpos))))]
                xpos+=actsize
            ypos+=actsize
        #Now do fitting step...
        if subflag!=None:
            actmap*=subflag
        converged=0
        while converged==0:
            for i in xrange(self.nact):
                for j in xrange(self.nact):
                    if subflag!=None and subflag[i,j]==0:
                        continue
                    dm=self.fit(actmap)
                    rms=self.calcRMS(dm-phase,pupmask)
                    rms2=0.
                    #first step the actuator in +ve direction.
                    s=step
                    while rms2<rms:#improvement
                        actmap[i,j]+=s
                        dm=self.fit(actmap)
                        rms2=self.calcRMS(dm-phase,pupmask)
                        if rms2<rms:#improvement, going in right direction
                            rms=rms2
                            s*=1.1
                        else:#not improvement - reduce step, and try again.
                            s*=0.5
                            actmap[i,j]-=s

    def calcRMS(self,phase,mask=None):
        """calcs rms in radians (or whatever phase is in)"""
        if type(mask)==type(None):
            rms=numpy.sqrt(numpy.average(phase.flat**2)-numpy.average(phase.flat)**2)
        else:
            p=(phase*mask).ravel()
            p2=numpy.sum(p*p)
            s=numpy.sum(mask)
            m=numpy.sum(p)/s
            rms=numpy.sqrt(p2/s-m*m)
        return rms

    def rotate(self,angle,phase=None):
        """Angle in degrees"""
        print "Rotating by %g (util.dm)"%angle
        overwrite=0
        res=None
        if phase==None:
            phase=self.phsOut
            overwrite=1
        else:
            res=numpy.zeros(phase.shape,numpy.float32)
        if phase.dtype!=numpy.float32:
            print "Copying phase to float32"
            phase=phase.astype(numpy.float32)
        if res==None:#rotate inplace
            cmod.utils.rotateArray(phase,angle)
            res=phase
        else:
            cmod.utils.rotateArray(phase,angle,res)
        return res

    def rotateOld(self,angle,phase=None):
        """angle in degrees..."""
        #the scipy version leaves blank lines when rotating about 90 degrees...
        #Also returns a uint8 result.  Doesn't even get 0 degrees correct.
        #return scipy.misc.pilutil.imrotate(phase,angle)
        #So, write my own version... which works well...
        overwrite=0
        if phase==None:
            phase=self.phsOut
            overwrite=1
        yr=numpy.arange(phase.shape[0])#-phase.shape[0]/2.+0.5
        xr=numpy.arange(phase.shape[1])#-phase.shape[0]/2.+0.5
        c=numpy.cos(angle/180.*numpy.pi)
        s=numpy.sin(angle/180.*numpy.pi)
        res=numpy.zeros(phase.shape,phase.dtype)
        for yy in yr:
            y=yy-phase.shape[0]/2.+0.5
            for xx in xr:
                x=xx-phase.shape[1]/2.+0.5
                yold=-s*x+c*y + phase.shape[0]/2.-0.5
                xold=c*x+s*y + phase.shape[1]/2.-0.5
                x1=int(numpy.floor(xold))
                x2=x1+1
                y1=int(numpy.floor(yold))
                y2=y1+1
                xm=xold-x1
                ym=yold-y1
                if y2==phase.shape[0]:
                    y2=y1
                    ym=0
                if x2==phase.shape[1]:
                    x2=x1
                    xm=0
                if x1==-1:
                    x1=0
                    xm=0
                if y1==-1:
                    y1=0
                    ym=0
                #print y1,x1,ym,xm
                val=0
                if y1>=0 and y1<phase.shape[0]:
                    if x1>=0 and x1<phase.shape[1]:
                        val+=phase[y1,x1]*(1-xm)*(1-ym)
                    if x2<phase.shape[1] and x2>=0:
                        val+=phase[y1,x2]*xm*(1-ym)
                if y2<phase.shape[0] and y2>=0:
                    if x2<phase.shape[1] and x2>=0:
                        val+=phase[y2,x2]*xm*ym
                    if x1>=0 and x1<phase.shape[1]:
                        val+=phase[y2,x1]*(1-xm)*ym
                #print yy,xx,val
                res[yy,xx]=val
        if overwrite:
            self.phsOut[:]=res
        return res

# # # # # # # # # # #  END OF CLASS MirrorSurface  # # # # # # # # # # # # # # # # # # # 

# # The functions below are not used anywhere in the aosim:  # # # # # # # # # # # # # # 

def dmProjectionQuick(config=None,batchno=0,vdmidstr="vdm",rmx=None,rmxOutName=None,reconIdStr=None):
    """Uses vdmUser to do a geometrical projection"""
    if config==None:
        config="params.xml"
    if type(config) in [type(""),type([]),type(())]:
        import base.readConfig
        config=base.readConfig.AOXml(config,batchno=batchno)
    if config.getVal("projMxFilename",raiseerror=0)==None:
        config.this.globals.projMxFilename="projmx%d.fits"%batchno
    if reconIdStr!=None:#should be specified in config file, but if not, can specify here...
                        #(though may still be overwritten by config file one)
        config.this.globals.reconIdStr=reconIdStr
    import science.vdmUser
    v=science.vdmUser.vdmUser(None,config,idstr=vdmidstr)
    #the projection matrix will have been saved now.
    if type(rmx)==type(""):
        import util.FITS
        rmx=util.FITS.Read(rmx)[1]
    if rmx!=None:#apply projection matrix to create specific reconstructor.
        #res=numpy.empty((rmx.shape[0],v.projmx.shape[0]),numpy.float32,order='F')
        try:
            import cmod.mkl
        except:
            print "Unable to import cmod.mkl - using quick.dot instead"
            res=quick.dot(rmx,v.projmx.T)
        else:
            res=numpy.empty((rmx.shape[0],v.projmx.shape[0]),numpy.float32,order='F')
            cmod.mkl.gemm(rmx,v.projmx.T,res)
        if rmxOutName==None:
            rmxOutName="rmxProjected%d.fits"%batchno
        util.FITS.Write(res,rmxOutName)
    return v.projmx

def calcActuators(hlist,fov,telDiam,r0list=None,strList=None):
    """Computes the ideal ratio of number of actuators for multi-DMs.
    r0list is computed from globalRo*strLayer**(-3./5) I think.  Or maybe is just the strengths of the layers.  Or 1/strength of layers.  Oops - I forget!!!  Use strlayer**-0.6.  If r0List==None, then strList is used and assumes this.
    This is taken from a 2 page paper by Don Gavel (google deformable mirror fitting error)
    returns the ratios of number of actuators in 1 dimension (nact) required.
    """
    if r0list==None:#strList cannot be None in this case.
        r0list=map(lambda x:x**-0.6,strList)
    mu=1.
    sigma2=1.
    M=len(hlist)
    alpha=mu*(numpy.pi/4)**(5./6)
    dlist=numpy.zeros((M,),numpy.float32)
    diam=numpy.zeros((M,),numpy.float32)
    nact=numpy.zeros((M,),numpy.float32)
    for i in range(M):
        diam[i]=telDiam+2*hlist[i]*numpy.tan(fov/3600./180.*numpy.pi)
        dlist[i]=diam[i]/r0list[i]/diam[0]*r0list[0]
    nact[0]=(1./sigma2*alpha*dlist.sum()**(10./11))**(6./5)
    for i in range(1,M):
        nact[i]=(dlist[0]**(5./3)*nact[0]**(-11./6)/dlist[i]**(5./3))**(-6./11)
    nact/=nact[0]
    nact=numpy.sqrt(nact)
    return nact,diam#diam/nact gives actSpacing (normalised)...

def dmProjection(dm,vdmList,npup,telDiam,interpolate=1,rcond=1e-15,usemkl=0,usemmap=0,stage=0,basefilename="tmpinfluence.fits",nthreads=1,pupfn=None):
    """Here is TB code...
    This method makes this runnable for ELT scale simulations.
    User can choose to start at different stages - ie to recover from a crash or something.
    Stage=0 - generate influence
    stage=1 - svd to get invInf

   inf_funcs=[]
   for idm in range(N_DM):
       inf_funcs.append(numpy.zeros((DM_TOTNACT[idm],NPUP**2),numpy.float))

   for idm in range(N_DM):
       r=DM_H[idm]*numpy.tan(deg2rad(WFS_OFFAXIS[iwfs]/3600.))
       pr_dx=r*numpy.cos(deg2rad(WFS_PHI[iwfs]))/APEL_SIZE
       pr_dy=r*numpy.sin(deg2rad(WFS_PHI[iwfs]))/APEL_SIZE

       x=mynint(DM_NPUP[idm]/2+pr_dx)
       y=mynint(DM_NPUP[idm]/2+pr_dy)

       map=numpy.zeros((DM_NACT[idm],DM_NACT[idm]),numpy.float)
       dmphs=numpy.zeros((DM_NPUP[idm],DM_NPUP[idm]),numpy.float)
       iact=0
        for i in range(DM_NACT[idm]):
           for j in range(DM_NACT[idm]):
               map*=0.
                map[i,j]=1.
                DM_OBJ[idm].fit(map,dmphs)
                inf_funcs[idm][iact]=dmphs[x-NPUP/2:x+NPUP/2,y-NPUP/2:y+NPUP/2].flatten()
                iact+=1

   projmx=numpy.zeros((DM_ALLTOTNACT,DM_TOTNACT[0]),numpy.float)
   iact=0
   tmp=numpy.linalg.pinv(inf_funcs[0])
   for idm in range(N_DM):
        projmx[iact:iact+DM_TOTNACT[idm]]=quick.dot(inf_funcs[idm],tmp)
        iact+=DM_TOTNACT[idm]

   return projmx
    """
    #First, compute the inverse part... then use this to do the rest.
    #This is done for the true MOAO DM.
    #import util.computeRecon
    #import cmod.mkl
    import time
    import os
    projmx=None
        
    if stage<1:
        t1=time.time()
        if usemmap==0:
            influence=numpy.zeros((dm.nacts,npup**2),numpy.float32)#this may need to be mmaped.
        else:
            mmap=numpy.memmap(basefilename+".mmap",dtype=numpy.float32,mode="w+",shape=(dm.nacts*npup**2+2880/4))
            influence=mmap[2880/4:]
            influence.shape=dm.nacts,npup**2
            hdr=mmap[:2880/4].view("c")
            hdr[:]=numpy.array(list(util.FITS.MakeHeader(influence.shape,"f",doByteSwap=0)))
        actmap=numpy.zeros((dm.nact,dm.nact),numpy.float32)
        #note, for this real DM, dm.dmpup==npup.
        phs=numpy.zeros((dm.dmpup,dm.dmpup),numpy.float32)
        pos=0
        ms=dm.getMirrorSurface()
        #This bit could be multi-processed quite easily.  But, note, only python2.5 on gig4X, so can't use multiprocessing module.  os.fork?
        if nthreads==1 or usemmap==0:
            for i in range(dm.nact):
                print "Doing dm: %d/%d      \r"%(i+1,dm.nact),
                sys.stdout.flush()
                for j in range(dm.nact):
                    if dm.dmflag[i,j]:
                        actmap[i,j]=1.
                        ms.fit(actmap,phs)
                        if pupfn!=None:
                            phs*=pupfn
                        #Here, we just select all of the phase since we assume that this true DM is of size npup anyway.  If this is too restrictive, will have to think a bit and recode slightly.
                        influence[pos]=phs.ravel()
                        pos+=1
                        actmap[i,j]=0.
        else:
            pidlist=[]
            for n in range(nthreads):
                pid=os.fork()
                if pid==0:#nth child
                    s=n*((dm.nact+nthreads-1)//nthreads)
                    e=(n+1)*((dm.nact+nthreads-1)//nthreads)
                    if e>dm.nact:
                        e=dm.nact
                    pos=dm.dmflag[:s].sum()
                    for i in range(s,e):
                        print "Thread %d doing dm: %d/%d     \r"%(n,i+1,dm.nact),
                        sys.stdout.flush()
                        for j in range(dm.nact):
                            if dm.dmflag[i,j]:
                                actmap[i,j]=1.
                                ms.fit(actmap,phs)
                                if pupfn!=None:
                                    phs*=pupfn
                                #Here, we just select all of the phase since we assume that this true DM is of size npup anyway.  If this is too restrictive, will have to think a bit and recode slightly.
                                influence[pos]=phs.ravel()
                                pos+=1
                                actmap[i,j]=0.
                    #Now, we've done our work, so exit.
                    sys.exit(0)
                else:
                    pidlist.append(pid)
            for pid in pidlist:
                os.waitpid(pid,0)
        #End of stage 0...
        if usemmap==0:#otherwise, already written...
            import util.FITS
            util.FITS.Write(influence,basefilename,doByteSwap=0)
        else:
            influence.flush()
        print "Stage 0 took %gs"%(time.time()-t1)
    if stage<2:
        #Now do the inverse of this... using mkl...
        t1=time.time()
        if stage==1:#need to load influence...
            if usemmap:
                influence=numpy.memmap(basefilename+".mmap",dtype=numpy.float32,mode="r",offset=2880)[:dm.nacts*npup**2]
                influence.shape=(dm.nacts,npup**2)
            else:
                influence=util.FITS.Read(basefilename)[1]
            #print min(influence.ravel()),max(influence.ravel())
        if usemkl:
            import util.FITS
            import util.computeRecon
            #del(influence)
            mr=util.computeRecon.makeRecon(basefilename,rcond)
            res=mr.dotdot(influence)
            u,e,vt=mr.dosvd(issparse=0,a=res)
            del(res)
            del(u)
            del(e)
            del(vt)
            inv=mr.makeInv()
            #Note, invInf might be too large to fit in memory - may need to be mmaped.
            if usemmap:
                #Need to rearrange influence so that it is fortran contiguous.  I think.  inv is already F contig (so invT is c contig)  OR maybe it doesn't matter
                print "TODO: rearrange influence so is F contig (current: %s)"%str(influence.flags).replace("\n",",")

                mmap=numpy.memmap(basefilename[:-5]+"_rmxden%g.fits"%rcond,dtype=numpy.float32,mode="w+",shape=(npup**2*dm.nacts+2880/4,))
                invInf=mmap[2880/4:]
                invInf.shape=npup**2,dm.nacts
                hdr=mmap[:2880/4].view("c")
                hdr[:]=numpy.array(list(util.FITS.MakeHeader(invInf.shape,"f",doByteSwap=0)))
                import cmod.mkl
                cmod.mkl.gemm(inv.transpose(),influence,invInf.transpose())
                del(inv)
                del(influence)
            else:
                del(inv)
                invInf=mr.denseDotDense(transA=1)
        else:
            #invInf=numpy.linalg.pinv(influence,rcond)
            #do it this way, so that gives same results as mkl...
            invInf=quick.dot(numpy.linalg.pinv(quick.dot(influence,influence.transpose()),rcond),influence).transpose()
            util.FITS.Write(invInf,basefilename[:-5]+"_rmxden%g.fits"%rcond,doByteSwap=0)
        print "Stage 1 took %gs"%(time.time()-t1)
        #print max(invInf.ravel())

    if stage<3:
        t1=time.time()
        #And now compute the projection matrix - this is done for each of the virtual DMs.
        #This part, we may consider forking to do multi-process.  So, arrays need to be shm/memmapped.
        #The big arrays that will be needed are:
        #invInf
        #projmx
        if stage==2:
            #have to load stuff computed previously.
            if usemmap:
                invInf=numpy.memmap(basefilename[:-5]+"_rmxden%g.fits"%rcond,dtype=numpy.float32,mode="r",offset=2880)[:npup**2*dm.nacts]
                invInf.shape=npup**2,dm.nacts
            else:
                invInf=util.FITS.Read(basefilename[:-5]+"_rmxden%g.fits"%rcond)[1]
        #it would be a good idea to make invInf F contiguous.  So, do the swap here...
        
                
        nactTotAll=0
        for vdm in vdmList:
            nactTotAll+=vdm.nacts
        if nthreads>1:
            projmx=numpy.memmap("/dev/shm/projmx.mmap",dtype=numpy.float32,mode="w+",shape=(nactTotAll,dm.nacts),order='F')
            os.unlink("/dev/shm/projmx.mmap")
            projmx[:]=0
        else:
            projmx=numpy.zeros((nactTotAll,dm.nacts),numpy.float32,order='F')
        print "total number of actuators:",nactTotAll
        from util.computeRecon import getMem
        mem=getMem()
        memleft=mem-projmx.size*projmx.itemsize-max(vdm.dmpup for vdm in vdmList)**2*4-max(vdm.nact for vdm in vdmList)**2*4
        if memleft>0:
            nblock=int(numpy.ceil((invInf.size*invInf.itemsize)/float(memleft)))
            if nblock<1:
                nblock=1
        else:
            print "Warning - used all memory - will be slow"
            nblock=1
        if nblock>1:
            nblock*=2
        print mem,memleft,projmx.shape,invInf.shape,nblock,invInf.flags
        bstart=0
        bend=0
        invInfPart=numpy.empty((npup**2,invInf.shape[1]%nblock+invInf.shape[1]/nblock),numpy.float32,order='F')
        for ii in range(nblock):
            bstart=bend
            bend+=invInf.shape[1]/nblock
            if ii==nblock-1:
                bend=invInf.shape[1]
            invInf2=invInfPart[:,:bend-bstart]
            #Now copy invInf (which will convert it from C to F contig
            print "Copying invInf to F config"
            invInf2[:]=invInf[:,bstart:bend]
            print "Doing block %d/%d"%(ii+1,nblock)
            tlist=[]
            for i in range(nthreads):
                #projectionWorker(vdmList,projmx,invInf,nblock,i,nthreads,interpolate,pupfn)
                tlist.append(threading.Thread(target=projectionWorker,args=(vdmList,projmx[:,bstart:bend],invInf2,nblock,i,nthreads,interpolate,pupfn,usemkl)))
                tlist[-1].start()
                #projectionWorker(vdmList,projmx[:,bstart:bend],invInf2,nblock,i,nthreads,interpolate,pupfn)
            for t in tlist:
                t.join()
        print "Stage 2 took %gs"%(time.time()-t1)
    return projmx

def projectionWorker(vdmList,projmx,invInf,nblock,threadno,nthreads,interpolate,pupfn,usemkl):
    pos=0

    #Now, to avoid excessing swapping, calculate the number of blocks invInf should be broken into...
    #We will then have to calculate the surfaces this many times, but still is probably faster.
    print "Running process %d"%threadno
    # bstart=0
    # bend=0#invInf.shape[1]/nblock
    # for ii in range(nblock):
    #     bstart=bend
    #     bend+=invInf.shape[1]/nblock
    #     if ii==nblock-1:
    #         bend=invInf.shape[1]
    ii=999
    if 1:
        pos=0

        for i in range(len(vdmList)):
            vdm=vdmList[i]
            ms=vdm.getMirrorSurface()
            actmap=numpy.zeros((vdm.nact,vdm.nact),numpy.float32)
            phs=numpy.zeros((vdm.dmpup,vdm.dmpup),numpy.float32)
            r=vdm.height*numpy.tan(dm.primaryTheta*numpy.pi/180/3600.)
            xpos=r*numpy.cos(dm.primaryPhi*numpy.pi/180.)*npup/telDiam+vdm.dmpup/2.-npup/2.
            ypos=r*numpy.sin(dm.primaryPhi*numpy.pi/180.)*npup/telDiam+vdm.dmpup/2.-npup/2.
            if interpolate==0:
                xf=int(round(xpos))
                yf=int(round(ypos))
                xt=xf+npup
                yt=yf+npup
            else:
                import cmod.interp
                xf=int(numpy.floor(xpos))
                yf=int(numpy.floor(ypos))
                xfrac=xpos-xf
                yfrac=ypos-yf
                xt=xf+npup+1
                yt=yf+npup+1
                interpolated=numpy.zeros((npup,npup),numpy.float32)
                xin=numpy.arange(npup+1).astype(numpy.float64)
                xout=numpy.arange(npup).astype(numpy.float64)
                yout=numpy.arange(npup).astype(numpy.float64)
                xout+=xfrac
                yout+=yfrac

            nstart=(vdm.nact/nthreads)*threadno
            nend=(vdm.nact/nthreads)*(threadno+1)
            if nend>vdm.nact:
                nend=vdm.nact
            if threadno==nthreads-1:#last thread
                nend=vdm.nact
            #and move to the right offset...
            pos+=vdm.dmflag[:nstart].sum()
            print pos,nstart,nend
            for j in xrange(nstart,nend):
                if threadno==nthreads/2:
                    print "Thread %d doing block %d/%d vdm %d/%d: %d/%d       \n"%(threadno,ii+1,nblock,i+1,len(vdmList),j+1,nend),
                    sys.stdout.flush()
                for k in xrange(vdm.nact):
                    if vdm.dmflag[j,k]:
                        if threadno==nthreads/2:
                            print "Doing %d/%d"%(k+1,vdm.nact)
                        actmap[j,k]=1.
                        ms.fit(actmap,phs)
                        #TODO: Optionally, here, we might do a subpxl interpolation
                        #print xf,xt,yf,yt,phs[yf:yt,xf:xt].shape,invInf.shape,phs.shape
                        #print max(phs.ravel())
                        #If invInf is larger than memory, may make sense to do in sections invInf[:,f:t], generating each phs several times to avoid swapping.
                        if interpolate==0:#just select the correct bit of phase
                            interpolated=phs[yf:yt,xf:xt]
                        else:
                            #   Comment, UB: ProjectionWorker is not used in aosim and is called only
                            #   interactively, therefore the number of threads is hard-coded and set to 1:
                            cmod.interp.gslCubSplineInterp(phs[yf:yt,xf:xt],xin,xin,yout,xout,interpolated,1)
                        if pupfn!=None:
                            interpolated*=pupfn
                        #tmp=quick.dot(interpolated.ravel(),invInf[:,bstart:bend])
                        #print max(tmp)
                        #print tmp.shape
                        #projmx[pos,bstart:bend]=tmp

                        #tmp=numpy.zeros((1,bend-bstart),numpy.float32,order='F')
                        tmp=numpy.zeros((1,projmx.shape[1]),numpy.float32,order='F')
                        intflat=interpolated.ravel()
                        intflat.shape=1,intflat.size
                        #cmod.mkl.gemm(intflat,invInf[:,bstart:bend],tmp)
                        if usemkl:
                            cmod.mkl.gemm(intflat,invInf,tmp)
                        else:
                            tmp=quick.dot(intflat,invInf)
                        #projmx[pos,bstart:bend]=tmp
                        projmx[pos]=tmp
                        pos+=1
                        actmap[j,k]=0.
            pos+=vdm.dmflag[nend:].sum()




if __name__=="__main__":
    import sys
    import time
    if sys.argv[1]=="projection":
        #Note on execution time - doubling npu quadruples the runtime of this calc.  Therefore, if you can get away with undersampling npup for this part do.  The difference in projection matrix is small (0.1-1%), but you should see for your own case whether this makes a difference.
        t1=time.time()
        print sys.argv
        #calculate DM projection
        pfile=sys.argv[2].split(",")#params.xml
        batchno=int(sys.argv[3])#15
        dmid=sys.argv[4]#pdmm
        vdmidList=eval(sys.argv[5])#['0H','10000H']
        outfile=sys.argv[6]#/var/ali/spmx2layer42mb30dense_rmxden1e-05projm.fits
        rcond=float(sys.argv[7])#1e-15
        usemkl=int(sys.argv[8])#set to 1 if the problem is large.
        usemmap=int(sys.argv[9])#1
        stage=int(sys.argv[10])#0
        nthreads=int(sys.argv[11])#1
        fitpup=int(sys.argv[12])
        if len(sys.argv)>13:
            rmxfile=sys.argv[13]
        else:
            rmxfile=None
        import base.readConfig
        import util.FITS
        config=base.readConfig.AOXml(pfile,batchno=batchno)
        npup=config.getVal("npup")
        telDiam=config.getVal("telDiam")
        dmObj=config.getVal("dmObj")
        pupil=config.getVal("pupil")
        config.setSearchOrder(["tomoRecon","globals",])
        if rmxfile==None:
            rmxfile=config.getVal("reconmxFilename",raiseerror=0)
        if rmxfile=="":
            rmxfile=None
        basefilename=outfile[:-5]+"_base.fits"
        dmObj.computeDMPupil(dmid,pupil.r2,retPupil=1)
        pupfn=None
        if fitpup==1:
            pupfn=config.getVal("pupil").fn
        dm=dmObj.getDM(dmid)
        vdmList=[]
        for vdmid in vdmidList:
            dmObj.computeDMPupil(vdmid,pupil.r2,retPupil=0)
            vdmList.append(dmObj.getDM(vdmid))
        interpolate=1
        if nthreads!=1:
            print "Warning - dones't work with >1 thread"
        projmx=dmProjection(dm,vdmList,npup,telDiam,interpolate,rcond,usemkl,usemmap,stage,basefilename,nthreads,pupfn)
        if rmxfile!=None:
            print "Producing reconstructor..."
            if usemkl:
                import util.computeRecon
                mr=util.computeRecon.makeRecon(basefilename,rcond)
                pfile=basefilename[:-5]+"_projmx.fits"
                util.FITS.Write(projmx,pfile)
                projmx=mr.denseDotDense(a=rmxfile,b=pfile,resname=outfile)
                print projmx.shape
            else:
                rmx=util.FITS.Read(rmxfile)[1]
                if rmx.shape[1]==projmx.shape[0]:
                    projmx=quick.dot(rmx,projmx)
                else:
                    print "Failed to dot rmx and projmx: %s %s"%(str(rmx.shape),str(projmx.shape))
        else:
            print "Failed to find rmx filename - unable to dot with projmx"
        print "Writing to %s shape %s"%(outfile,str(projmx.shape))
        util.FITS.Write(projmx,outfile)
        t2=time.time()
        print "Projection matrix took %gs"%(t2-t1)
