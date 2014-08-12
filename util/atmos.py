from cmod.interp import gslCubSplineInterp,linearshift
import time
import numpy
import util.dm,util.guideStar
import util.tel
import util.zernikeMod
class layer:
    """a holder for info about a single layer"""
    def __init__(self,height,direction,speed,strength,seed,strtype="fraction"):
        """strtype can be fraction of cn2, depending on what the values are in strength"""

        self.height=height#m
        self.direction=direction#degrees
        self.speed=speed#m/s
        if strtype=="cn2":
            strength=strength**(1./2)
            print "Converting layer strength Cn^2 value to fraction (TODO: check in util.atmos.py that the power is correct - currently I sqrt the Cn2 value - may depend on what the units of your CN2 value are... some units allow you to use fraction for the strtype)."
        self.strength=strength#fraction
        self.seed=seed#random number seed
        self.zenith=None
        self.zenithAz=None
    def adjustForZenith(self,zenith,zenithAz):
        if self.zenith!=None:
            raise Exception("Zenith already adjusted for")
        self.zenith=zenith
        self.zenithAz=zenithAz
        if zenith==0.:
            return
        degrad=numpy.pi/180.
        self.height/=numpy.cos(zenith*degrad)
        #Now adjust the velocity and direction...
        vx=self.speed*numpy.cos((zenithAz-self.direction)*degrad)*numpy.cos(zenith*degrad)
        vy=self.speed*numpy.sin((zenithAz-self.direction)*degrad)
        self.speed=numpy.sqrt(vx*vx+vy*vy)
        self.direction=zenithAz-numpy.arctan(numpy.tan((zenithAz-self.direction)*degrad)/numpy.cos(zenith*degrad))/degrad
        
class source:
    """a holder for source direction objects
    sourcelam is the wavelength for this direction.
    phslam is the wavelength at which the received phase is at.  Defaults to the same as sourcelam, but can be different.
    lam in nm.
    """
    def __init__(self,idstr,theta,phi,alt,nsubx,sourcelam=None,reconList=None,phslam=None):
        self.idstr=idstr#id string
        self.theta=theta#arcsec
        self.phi=phi#degrees
        self.alt=alt#-1 or lgs alt
        self.nsubx=nsubx#None or the number of subaps.
        self.sourcelam=sourcelam
        if phslam==None:
            self.phslam=sourcelam
        else:
            self.phslam=phslam
        if reconList!=None:
            if type(reconList)!=type([]):
                reconList=[reconList]
        self.reconList=reconList#a list of reconstructor idstr which use this wfs.
class geom:
    """A class to hold information about atmospheric/source geometry.  Typically used by the parameter file."""
    def __init__(self,layerDict,sourceList,ntel,npup,telDiam,r0,l0,zenith=0.,zenithAz=0.):
        """layerDict is a dictionary with keys equal to the layer name (the idstr used by infScrn objects)
        and values equal to either a layer() object or (depreciated) to a tuple of:
        (height, wind direction, speed, strength,initSeed) (metres, degrees, metres per second, fraction, int).

        sourceList is a list of source objects or a tuple of:
        (idstr,theta, phi, alt, nsubx) (idstr, arcsec, degrees, metres or -1 for star at infinity, number of subaps or None if not a GS).  Note, if the source is science only, nsubx should be None.  If the source is GS and science, nsubx should not be None.
        Here, idstr is the idstr used by infAtmos objects.  It is specified as a list since ordering may be important (eg for the reconstructor).
        zenith is in degrees.
        zenithAz the azimuthal angle, is in degrees too.
        """
        self.zenith=zenith
        self.zenithAz=zenithAz
        self.rotateDirections=0#Used in getScrnXPxls etc, whether to rotate the phase screens when considering max width to use.
        self.r0=r0*numpy.cos(self.zenith*numpy.pi/180.)**0.6#adjust r0 for zenith.
        self.l0=l0
        self.zenithOld=0.#this was used when making eliptical pupils.  However, that idea is no longer in vogue.  So, this is always set at zero.
        totstr=0
        self.layerDict=layerDict
        for key in self.layerDict.keys():
            if type(self.layerDict[key])==type(()):
                print "DEPRECIATION WARNING - atmosGeom object layerDict should have values equal to util.atmos.layer() objects"
                self.layerDict[key]=layer(*self.layerDict[key])
            totstr+=self.layerDict[key].strength
            self.layerDict[key].adjustForZenith(self.zenith,self.zenithAz)
        if totstr!=1.:
            print "WARNING: atmospheric layer strenghts do not add up to 1 - correcting..."
            for key in self.layerDict.keys():
                self.layerDict[key].strength/=totstr
        self.sourceList=[]#sourceList
        for obj in sourceList:
            if type(obj)==type(()):
                print "DEPRECIATION WARNING - atmosGeom object sourceList should be a list of util.atmos.source() objects"
                self.sourceList.append(source(*obj))
            else:
                self.sourceList.append(obj)
        
        self.sourceDict={}
        self.sourceOrder=[]
        for s in self.sourceList:
            if self.zenith!=0 and s.alt>0:
                print "WARNING Changing source height to take zenith into account"
                s.alt/=numpy.cos(self.zenith*numpy.pi/180.)
            if self.sourceDict.has_key(s.idstr):
                print s.idstr,self.sourceDict.keys()
                raise Exception("util.atmos.geom - source dictionary defines an idstr more than once")
            self.sourceDict[s.idstr]=s
            self.sourceOrder.append(s.idstr)
        self.npup=npup
        self.ntel=ntel
        self.telDiam=telDiam
        self.scrnSize=None
        self.layerOffset=None
    def getSource(self,id,raiseerror=0):
        if self.sourceDict.has_key(id):
            return self.sourceDict[id]
        if raiseerror:
            raise Exception("source with id %s not found in sourceDict"%id)
        else:
            key=self.sourceDict.keys()[0]
            
            print "WARNING: util.atmos.geom - source with id %s not found - using %s instead."%(id,key)
            return self.sourceDict[key]
    def layerHeight(self,id):
        return self.layerDict[id].height#/numpy.cos(self.zenith*numpy.pi/180.)
    def layerWind(self,id):
        return self.layerDict[id].direction
    def layerSpeed(self,id):
        return self.layerDict[id].speed
    def windDirections(self):
        w={}
        for key in self.layerDict.keys():
            w[key]=self.layerWind(key)
        return w
    def windSpeeds(self):
        s={}
        for key in self.layerDict.keys():
            s[key]=self.layerSpeed(key)
        return s
    def layerAltitudes(self):
        a={}
        for key in self.layerDict.keys():
            a[key]=self.layerHeight(key)
        return a
    def layerStrength(self,id):
        return self.layerDict[id].strength
    def layerInitSeed(self,id):
        return self.layerDict[id].seed

    def sourceTheta(self,id):
        return self.getSource(id).theta
    def sourcePhi(self,id):
        return self.getSource(id).phi
    def sourceAlt(self,id):
        return self.getSource(id).alt
    def sourcensubx(self,id):
        return self.getSource(id).nsubx
    def sourceLambda(self,id):
        return self.getSource(id).sourcelam
    def phaseLambda(self,id):
        return self.getSource(id).phslam
    def sourceThetas(self):
        t={}
        for key in self.sourceDict.keys():
            t[key]=self.sourceTheta(key)
        return t
    def getLayerWidth(self,height):
        """computes layer width at a given height...
        NOT USED ANYWHERE - AND NOT RECOMMENDED IF ZENITH!=0
        """
        if self.zenith!=0:
            print "WARNING util.atmos.getLayerWidth - do not use if zenith!=0"
        fov=max(self.sourceThetas().values())/60./60./180*numpy.pi
        return numpy.tan(fov)*height*2/numpy.cos(self.zenith*numpy.pi/180.)+self.telDiam
    def getScrnSize(self,rotateDirections=0):
        """Compute the screen sizes"""
        self.rotateDirections=rotateDirections
        thetas={}#source direction
        phis={}#source direction
        ntel=self.ntel
        npup=self.npup
        telDiam=self.telDiam
        arcsecrad=2*numpy.pi/360./3600.
        degrad=2*numpy.pi/360.
        scrnSize={}
        rotangle=0.
        for altkey in self.layerDict.keys():
            xposlist=[]
            yposlist=[]
            layerAlt=self.layerHeight(altkey)
            if rotateDirections:#this comes from the new phase screen method (generated along 1 axis and then rotated in iatmos).  To work out the max size of this screen, we consider the rotation here.  
                wd=self.layerWind(altkey)
                rotangle=(wd+90)*numpy.pi/180.#Add 90 because they are generated going upwards, but in simulation, we define 0 to be wind travelling from right to left.
            for key in self.sourceDict.keys():
                #xposlist.append(layerAlt*Numeric.fabs(Numeric.tan(self.sourceTheta(key)*arcsecrad)*Numeric.cos(self.sourcePhi(key)*degrad)))
                #yposlist.append(layerAlt*Numeric.fabs(Numeric.tan(self.sourceTheta(key)*arcsecrad)*Numeric.sin(self.sourcePhi(key)*degrad)))
                xposlist.append(layerAlt/numpy.cos(self.zenithOld*degrad)*numpy.tan(self.sourceTheta(key)*arcsecrad)*numpy.cos(self.sourcePhi(key)*degrad+rotangle))
                yposlist.append(layerAlt*numpy.tan(self.sourceTheta(key)*arcsecrad)*numpy.sin(self.sourcePhi(key)*degrad+rotangle))
            maxx=max(xposlist)
            minx=min(xposlist)
            maxy=max(yposlist)
            miny=min(yposlist)
            #scrnXPxls=int(Numeric.ceil(maxx*ntel/telDiam+npup+Numeric.ceil(minx*ntel/telDiam))+1)
            #scrnYPxls=int(Numeric.ceil(maxy*ntel/telDiam+npup+Numeric.ceil(miny*ntel/telDiam))+1)
            print altkey,maxx,minx,(maxx-minx),maxy,miny,(maxy-miny)
            scrnXPxls=int(numpy.ceil(npup/numpy.cos(self.zenithOld*degrad)))+1+int(numpy.ceil((maxx-minx)*ntel/telDiam))#agb 090313 - changed from ceil(maxx-minx) to ceil of whole thing. 090518 added zenith part.
            scrnYPxls=npup+1+int(numpy.ceil((maxy-miny)*ntel/telDiam))
            scrnSize[altkey]=(scrnXPxls,scrnYPxls)
        print "scrnSize: %s"%str(scrnSize)
        self.scrnSize=scrnSize
        return scrnSize

    def getScrnXPxls(self,id,rotateDirections=0):
        """return screen x size"""
        if self.scrnSize==None or self.rotateDirections!=rotateDirections:
            self.getScrnSize(rotateDirections)
        return self.scrnSize[id][0]
    def getScrnYPxls(self,id,rotateDirections=0):
        """return screen y size"""
        if self.scrnSize==None or self.rotateDirections!=rotateDirections:
            self.getScrnSize(rotateDirections)
        return self.scrnSize[id][1]
    def getLayerOffset(self):
        """compute the layer offsets"""
        arcsecrad=2*numpy.pi/360./3600.
        degrad=2*numpy.pi/360.
        npup=self.npup
        ntel=self.ntel
        telDiam=self.telDiam
        layerXOffset={}
        layerYOffset={}
        for altKey in self.layerDict.keys():
            xpos=[]
            ypos=[]
            for sourceKey in self.sourceDict.keys():
                xpos.append(numpy.tan(self.sourceTheta(sourceKey)*arcsecrad)*numpy.cos(self.sourcePhi(sourceKey)*degrad)*self.layerHeight(altKey)/numpy.cos(self.zenithOld*degrad)/telDiam*ntel-npup/numpy.cos(self.zenithOld*degrad)/2.)#scrnSize[altKey][0]/2.)
                ypos.append(numpy.tan(self.sourceTheta(sourceKey)*arcsecrad)*numpy.sin(self.sourcePhi(sourceKey)*degrad)*self.layerHeight(altKey)/telDiam*ntel-npup/2.)#scrnSize[altKey][1]/2.)
            minx=min(xpos)
            miny=min(ypos)
            #print "minx,miny %g %g"%(minx,miny)
            if minx<0.:#less than zero...
                layerXOffset[altKey]=int(numpy.ceil(numpy.fabs(minx)))
            else:
                layerXOffset[altKey]=-int(numpy.floor(minx))
            if miny<0.:#less than zero...
                layerYOffset[altKey]=int(numpy.ceil(numpy.fabs(miny)))
            else:
                layerYOffset[altKey]=-int(numpy.floor(miny))
        self.layerOffset=(layerXOffset,layerYOffset)
        #print self.layerOffset
        return self.layerOffset

    def getLayerXOffset(self):
        """return the layer x offsets"""
        if self.layerOffset==None:
            self.getLayerOffset()
        return self.layerOffset[0]
    def getLayerYOffset(self):
        """return the layer y offsets"""
        if self.layerOffset==None:
            self.getLayerOffset()
        return self.layerOffset[1]
#makeDMList depreciated - use util.dm instead...
##     def makeDMList(self,nact,height,fov=None,coupling=0.1):
##         """Makes a list of util.dm.DM objects ready for theoretical poke matricees etc.
##         nact can be an int or a list of same length as height.  Similarly for fov.
##         """
##         if type(nact)!=type([]):
##             nact=[nact]*len(height)
##         if type(fov)!=type([]):
##             fov=[fov]*len(height)
##         dmlist=[]
##         for i in range(len(height)):
##             if fov[i]==None:
##                 fov[i]=max(self.sourceThetas().values())*2+0.1
##             dmlist.append(util.createPokeMx.DM(nact[i],height[i],fov[i],coupling))
##         return dmlist
    def makeNGSList(self,reconID=None,minarea=0.5):#,nsubx,keylist):
        """make a list of util.guideStar.NGS objects.  nsubx can be
        int, a dict or a list of length keylist length.
        keylist is the list of keys to use for the NGS objects.  Any
        with a +ve altitude won't be included.
        If reconID is specified, only GS with a matching reconList entry
        will be included.
        """
##         if keylist==None:
##             if type(nsubx)==type({}):
##                 keylist=nsubx.keys()
##             else:
##                 keylist=self.sourceDict.keys()
##         if type(nsubx)!=type({}):
##             if type(nsubx)!=type([]):
##                 nsubx=[nsubx]*len(keylist)
##             n={}
##             for i in range(len(keylist)):
##                 n[keylist[i]]=nsubx[i]
##             nsubx=n
        if reconID=="":
            reconID=None
        ngslist=[]
        for key in self.sourceOrder:
            nsubx=self.sourcensubx(key)#If None, this is a science object only!
            if self.sourceAlt(key)<0 and nsubx!=None:#a NGS.
                if reconID==None:#any will do
                    ngslist.append(util.guideStar.NGS(nsubx,self.sourceTheta(key),self.sourcePhi(key),minarea=minarea))
                else:
                    rl=self.getSource(key).reconList
                    if rl!=None and (rl==[] or reconID in rl):
                        ngslist.append(util.guideStar.NGS(nsubx,self.sourceTheta(key),self.sourcePhi(key),minarea=minarea))

        return ngslist
    def makeLGSList(self,reconID=None,minarea=0.5):#,nsubx,keylist):
##         if keylist==None:
##             if type(nsubx)==type({}):
##                 keylist=nsubx.keys()
##             else:
##                 keylist=self.sourceDict.keys()
##         if type(nsubx)!=type({}):
##             if type(nsubx)!=type([]):
##                 nsubx=[nsubx]*len(keylist)
##             n={}
##             for i in range(len(keylist)):
##                 n[keylist[i]]=nsubx[i]
##             nsubx=n
        if reconID=="":
            reconID=None
        lgslist=[]
        for key in self.sourceOrder:
            nsubx=self.sourcensubx(key)
            if self.sourceAlt(key)>=0.:
                if nsubx==None:
                    raise Exception("util.atmos.geom - LGS specified without nsubx")
                if reconID==None:
                    lgslist.append(util.guideStar.LGS(nsubx,self.sourceTheta(key),self.sourcePhi(key),self.sourceAlt(key),minarea=minarea))
                else:
                    rl=self.getSource(key).reconList
                    if rl!=None and (rl==[] or reconID in rl):
                        lgslist.append(util.guideStar.LGS(nsubx,self.sourceTheta(key),self.sourcePhi(key),self.sourceAlt(key),minarea=minarea))
                        
        return lgslist
    
    def getWFSOrder(self,reconID=None):
        """Gets the order in which WFS data should be used.  This is the same order as the
        objects makeNGSList()+makeLGSList() would be provided.
        """
        idlist=[]
        if reconID=="":
            reconID=None
        for key in self.sourceOrder:
            nsubx=self.sourcensubx(key)
            if self.sourceAlt(key)<0 and nsubx!=None:
                if reconID==None:
                    idlist.append(key)
                else:
                    rl=self.getSource(key).reconList
                    if rl!=None and (rl==[] or reconID in rl):
                        idlist.append(key)
        for key in self.sourceOrder:
            nsubx=self.sourcensubx(key)
            if self.sourceAlt(key)>=0.:
                if nsubx==None:
                    raise Exception("util.atmos.geom - LGS specified w/o nsubx")
                if reconID==None:
                    idlist.append(key)
                else:
                    rl=self.getSource(key).reconList
                    if rl!=None and (rl==[] or reconID in rl):
                        idlist.append(key)
        return idlist
    
class atmos:
    """A class to carry out computation of a pupil phase screen
    for a given source direction.  This can (possibly) be used stand-alone
    and is used as part of the AO simulation (infAtmos module)."""
    def __init__(self,parent,sourceAlt,sourceLam,sourceTheta,sourcePhi,npup,pupil,colAdd,rowAdd,layerAltitude,phaseScreens,scrnScale,layerXOffset,layerYOffset,layerList=None,zenith=0.,intrinsicPhase=None,storePupilLayers=0,computeUplinkTT=0,launchApDiam=0.35,ntel=None,telDiam=0.,interpolationNthreads=0,):
        """parent is a dictionary of parents (keys are the layers).  Source altitude, wavelength and position are given, and info about all the phase screen layers.
        storePupilLayers - if 1, a dictionary of layers along this line of sight will be created - unexpanded in the case of LGS - i.e. full pupil.
        computeUplinkTT - if 1 then uplink tip/tilt will be computed.
        """
        
        self.npup=npup
        if ntel==None:
            ntel=npup
        self.ntel=ntel
        self.telDiam=telDiam
        self.interpolationNthreads=interpolationNthreads
        self.pupil=pupil
        self.parent=parent
        self.colAdd=colAdd
        self.rowAdd=rowAdd
        self.sourceAlt=sourceAlt
        self.sourceLam=sourceLam
        self.sourceTheta=sourceTheta
        self.sourcePhi=sourcePhi
        self.intrinsicPhase=intrinsicPhase#added to the output each iteration
        degRad=2*numpy.pi/360
        arcsecRad=2*numpy.pi/360/3600
        self.layerList=layerList
        if self.layerList==None:
            self.layerList=self.parent.keys()
        self.zenith=zenith
        self.zenithOld=0.#used to be used when stretching screens - now no longer.
        self.layerAltitude=layerAltitude#these are pre-scaled by zenith.
        self.sortedLayerList=[]
        hh=sorted(self.layerAltitude.values())
        for h in hh:
            for key in self.layerList:
                if self.layerAltitude[key]==h:
                    self.sortedLayerList.append(key)
                    break
        print "atmos: Layers %s, sorted %s"%(self.layerList,self.sortedLayerList)
        if len(self.layerList)!=len(self.sortedLayerList):
            raise Exception("Error sorting layerList in atmos.py")

        self.phaseScreens=phaseScreens
        self.scrnScale=scrnScale#length per phase pixel (teldiam/ntel)
        self.layerXOffset=layerXOffset
        self.layerYOffset=layerYOffset
        # if colAdd<zero, we add new columns on the right of the array.
        # If rowAdd<zero, we add new rows at the bottom of the array (note, that if plotting in Gist, this is the top of the array).
        self.storePupilLayers=storePupilLayers
        self.computeUplinkTT=computeUplinkTT
        if self.computeUplinkTT and self.sourceAlt==-1:
            raise Exception("Cannot compute uplink for source at infinity")
        self.launchApDiam=launchApDiam
        if self.storePupilLayers or self.computeUplinkTT:
            self.uplinkPositionDict={}
            if self.storePupilLayers:
                self.storedLayer={}
                for key in self.layerList:
                    self.storedLayer[key]=numpy.zeros((self.npup,self.npup),numpy.float64)
            else:
                self.tmpPup=numpy.zeros((self.npup,self.npup),numpy.float64)
            if self.computeUplinkTT:
                self.uplinkTTDict={}
                self.calibrateUplinkTilt()
                self.uplinkPhs=numpy.zeros((self.npup,self.npup),numpy.float64)
                #self.distToNextLayer={}
                self.distToFocus={}
                for i in range(len(self.sortedLayerList)):
                    key=self.sortedLayerList[i]
                    alt=self.layerAltitude[key]
                    if alt<self.sourceAlt:
                        self.distToFocus[key]=self.sourceAlt-alt
                #for i in range(len(self.sortedLayerList)):
                #    key=self.sortedLayerList[i]
                #    alt=self.layerAltitude[key]
                #    if i<len(self.sortedLayerList)-1:
                #        keynext=self.sortedLayerList[i+1]
                #        altnext=self.layerAltitude[keynext]
                #    else:
                #        altnext=self.sourceAlt
                #    if altnext>self.sourceAlt:
                #        altnext=self.sourceAlt
                #    if alt<self.sourceAlt:
                #        self.distToNextLayer[key]=altnext-alt
                #        if self.distToNextLayer[key]<0:
                #            raise Exception("Distance to next layer negative: %s (key %s)"%(str(self.distToNextLayer),key))

        #Now, for each layer, compute where abouts in this layer we should take phase from.

        self.positionDict={}
        xpos=numpy.tan(self.sourceTheta*arcsecRad)*numpy.cos(self.sourcePhi*degRad)
        ypos=numpy.tan(self.sourceTheta*arcsecRad)*numpy.sin(self.sourcePhi*degRad)
        for key in self.layerList:
            if self.sourceAlt<0 or self.sourceAlt>=self.layerAltitude[key]:
                #use the layer...
                shape=self.phaseScreens[key].shape
                x=xpos*self.layerAltitude[key]/numpy.cos(self.zenithOld*numpy.pi/180.)/self.scrnScale-self.npup/numpy.cos(self.zenithOld*degRad)/2+self.layerXOffset[key]#+shape[0]/2#zenith added 090518
                y=ypos*self.layerAltitude[key]/self.scrnScale-self.npup/2+self.layerYOffset[key]#+shape[1]/2
                #if x<0 or y<0 or numpy.ceil(x+self.npup)>shape[1] or numpy.ceil(y+self.npup)>shape[0]:#agb removed +1 in x+npup and y+npup, and replaced with numpy.ceil, date 070831.
                if x<0 or y<0 or x+self.npup/numpy.cos(self.zenithOld*degRad)+1>shape[1] or y+self.npup+1>shape[0]:#zenith added 090518
                    print "ERROR: util.atmos - phasescreen %s is not large enough to contain this source %g %g %g %g %g %g"%(str(key),x,y,x+self.npup/numpy.cos(self.zenithOld*degRad)+1,y+self.npup+1,shape[1],shape[0])#agbc changed shape[0->1] and vice versa zenith added 090518
                    raise Exception("ERROR: util.atmos - phasescreen is not large enough to contain this source")
                axis1=axis2=axis3=None
                if (self.storePupilLayers or self.computeUplinkTT) and self.sourceAlt>0:
                    #Need to make a note of the position of full pupil.
                    width=self.npup
                    widthx=int(width/numpy.cos(self.zenithOld*numpy.pi/180.)+0.5)
                    if self.zenithOld!=0:
                        axis2=numpy.arange(width).astype(numpy.float64)
                        axis3=numpy.arange(widthx).astype(numpy.float64)*((width-1.)/(widthx-1.))
                    self.uplinkPositionDict[key]=(int(x),int(y),widthx,width,x-numpy.floor(x),y-numpy.floor(y),axis2,axis1,axis3)

                if self.sourceAlt>0:
                    #this is a LGS... need to project (interpolate)
                    #so we only select the central portion of the pupil
                    #since beam diverges from LGS spot.
                    npup2=self.npup*(1.-self.layerAltitude[key]/self.sourceAlt)
                    x=x+(self.npup-npup2)/numpy.cos(self.zenithOld*degRad)/2#select central portion zenith added 090518
                    y=y+(self.npup-npup2)/2
                    width=int(npup2+0.5)#round correctly...
                    widthx=int(npup2/numpy.cos(self.zenithOld*degRad)+0.5)
                    axis2=numpy.arange(width).astype(numpy.float64)
                    step=(width-1)/float(self.npup-1)
                    axis1=numpy.arange(self.npup).astype(numpy.float64)*step
                    #and the axis used for binning down to width when zenith!=0
                    if self.zenithOld!=0:
                        axis3=numpy.arange(widthx).astype(numpy.float64)*((width-1.)/(widthx-1.))
                else:#select a whole pupil area (collimated from NGS).
                    width=self.npup
                    widthx=int(width/numpy.cos(self.zenithOld*numpy.pi/180.)+0.5)
                    if self.zenithOld!=0:
                        axis2=numpy.arange(width).astype(numpy.float64)
                        axis3=numpy.arange(widthx).astype(numpy.float64)*((width-1.)/(widthx-1.))
                    
                self.positionDict[key]=(int(x),int(y),widthx,width,x-numpy.floor(x),y-numpy.floor(y),axis2,axis1,axis3)
                #print "positionDict %s: %s %s"%(str(key),str(self.positionDict[key]),str((self.layerXOffset[key],self.layerYOffset[key])))
            else:#layer above source, so don't use.
                self.positionDict[key]=()
        
    def initMem(self,outputData,interpPhs):
        """Initialise memory to shared memory regions... in this case, they
        are always the same size, so access directly without using
        arrayFromArray.
        """
        self.outputData=outputData    #this will be shared by all resource 
        self.interpPhs=interpPhs#sharing infAtmos objects.  

    def calibrateUplinkTilt(self):
        pup=util.tel.Pupil(self.npup,self.ntel/self.telDiam*self.launchApDiam,0)
        ztt=util.zernikeMod.Zernike(pup,3,computeInv=0).zern[1:]
        pupsum=pup.fn.sum()
        
        a=numpy.zeros((self.npup,self.npup),numpy.float32)
        a[:]=numpy.arange(self.npup)
        est=(ztt[0]*a).sum()/pupsum
        #est should equal theta=arctan((npup-1)/npup*(ntel/telDiam*launchApDiam)*500e-9/(2*numpy.pi)/launchApDiam)
        #Can use small angle approx - so est== (npup-1)*ntel*500e-9/(npup*telDiam*2*numpy.pi)
        #So, scaling factor is this/est.  Also include the pupil sum here, to avoid 2 further divisions each iteration later on.
        scalingFactor=((self.npup-1)*self.ntel*500e-9)/(self.npup*self.telDiam*2*numpy.pi*est*pupsum)
        #This means that given some phase, (phs*ztt[0]).sum()*scale should give the tilt in radians.
        #So, multiply scale with ztt, so then (phs*ztt[0]).sum() gives tilt in radians.
        ztt*=scalingFactor
        self.uplinkZTT=ztt

        #Now compute the zernikes that go from tip/tilt shift in m on-sky to radians across the pupil again.
        #Angle is arctan(dist/height) ~ dist/height.
        #Across the pupil this corresponds to a tilt of diam*tan(theta) = diam*dist/height, giving the P-V tilt in m.  diam*dist/height/500e-9*2*pi gives in radians.  ie this is the range of tt values that need to be applied.  But what is the range of the zernike generated?
        dtt=util.zernikeMod.Zernike(self.pupil.fn,3,computeInv=0).zern[1:]
        #Find range - same for tip and tilt.
        r=numpy.ptp(dtt[0])#numpy.max(dtt[0])-numpy.min(dtt[0])
        dtt*=self.telDiam*2*numpy.pi/(self.sourceAlt*500e-9*r)
        self.uplinkDownTT=-dtt
        print "Computed uplink related calibrations"

    def createPupilPhs(self,phaseScreens,interpPosColDict,interpPosRowDict,control):
        """Here, interpolate the screens, add together, and voila,
        have the pupil phase.
        Assumes the screens have already been created (eg in infAtmos).
        """
        #t1=time.time()
        if self.intrinsicPhase!=None:
            self.outputData[:]=self.intrinsicPhase
        else:
            self.outputData[:,]=0.#clear the array.
        if self.computeUplinkTT:
            self.uplinkPhs[:]=0#reset the uplink phase.
            tottip=0.
            tottilt=0.
        for key in self.sortedLayerList:#for each atmosphere layer... (increasing in height)
            posDict=self.positionDict[key]
            if len(posDict)>0:#this layer is used (ie below star height)
                #print "atmos time2 %g"%(time.time()-t1)
                interpPosCol=posDict[4]#static offset due to source direction.
                interpPosRow=posDict[5]#static offset due to source direction.
                if self.colAdd[key]>=0:#add on to start...
                    if interpPosColDict[key]>0:
                        interpPosCol+=1-interpPosColDict[key]
                else:
                    if interpPosColDict[key]>0:
                        interpPosCol+=interpPosColDict[key]
                    else:
                        interpPosCol+=1
                phsShiftY=numpy.floor(interpPosCol)
                tmpcol=interpPosCol
                interpPosCol-=phsShiftY

                if self.rowAdd[key]>=0:#add on to start
                    if interpPosRowDict[key]>0:
                        interpPosRow+=1-interpPosRowDict[key]
                else:
                    if interpPosRowDict[key]>0:
                        interpPosRow+=interpPosRowDict[key]
                    else:
                        interpPosRow+=1
                phsShiftX=numpy.floor(interpPosRow)
                tmprow=interpPosRow
                interpPosRow-=phsShiftX

                #now select the area we're interested in for this target...
                phs=phaseScreens[key][posDict[1]+phsShiftX:posDict[1]+posDict[3]+1+phsShiftX,
                                      posDict[0]+phsShiftY:posDict[0]+posDict[2]+1+phsShiftY]
                #now do interpolation... (sub-pxl)
                #print "atmos time3 %g"%(time.time()-t1)
                #This next section takes a long time (3/4)
                #Note, posDict[2,3] will be npup unless an LGS... in which case
                #will be the projected width, which we then interpolate out
                #to the full pupil array.
                #c function call about 10x faster.
                #print "Doing linshift"
                linearshift(phs,interpPosCol,interpPosRow,self.interpPhs[:posDict[3],:posDict[2]])
                #print "Linshift done"
                if self.storePupilLayers or self.computeUplinkTT:
                    #Want to store the individual layers at pupil size, before LGS compression.
                    #  e.g. for uplink tiptilts...
                    #Note - we do this a bit stupidly - ie brute force... the launch ap is only a small
                    #  fraction of the pupil, but we still use the full pupil, mostly multiplied by zeros... 
                    #  so if need speedup in future, this is one place...
                    if self.sourceAlt>0:#lgs
                        # want full pupil, so will have to shift
                        pd=self.uplinkPositionDict[key]
                        phs2=phaseScreens[key][pd[1]+phsShiftX:pd[1]+pd[3]+1+phsShiftX,
                                               pd[0]+phsShiftY:pd[0]+pd[2]+1+phsShiftY]
                        if self.storePupilLayers:
                            phs=self.storedLayer[key]
                        else:
                            phs=self.tmpPup
                        #print "Doing linearshift"
                        linearshift(phs2,interpPosCol,interpPosRow,phs)
                        #print "Linearshift done"
                    else:#already shifted
                        phs=self.interpPhs[:posDict[3],:posDict[2]]
                        if self.storePupilLayers:
                            self.storedLayer[key][:]=phs

                    if self.computeUplinkTT:
                        #add phase of this layer to the total uplink phase
                        self.uplinkPhs+=phs
                        #Now compute the combined tip/tilt of the uplink
                        #uplinkZTT contain a zernike tip and tilt mode on the diameter of the 
                        # launch telescope, scaled so that sum(phs*uplinkZTT[0]) gives tip in radians.
                        #tip,tilt=(numpy.inner(self.uplinkZTT[0].ravel(),self.uplinkPhs.ravel()),
                        # (numpy.inner(self.uplinkZTT[1].ravel(),self.uplinkPhs.ravel())))
                        #self.uplinkTTDict[key]=tip,tilt
                        #Now compute how this will shift the spot.
                        #hdiff=self.distToNextLayer[key]
                        #Instead, just compute tip/tilt of this layer.
                        tip,tilt=(numpy.inner(self.uplinkZTT[0].ravel(),phs.ravel()),
                                  (numpy.inner(self.uplinkZTT[1].ravel(),phs.ravel())))
                        self.uplinkTTDict[key]=tip,tilt
                        #Now compute how this will shift the spot.
                        hdiff=self.distToFocus[key]
                        #Assuming tip/tilt gives the mean tip and tilt in radians, then we 
                        # can compute the angle according to:  tan-1(tip*500e-9/(2*numpy.pi*launchApDiam))
                        tottip+=hdiff*numpy.tan(tip)
                        tottilt+=hdiff*numpy.tan(tilt)

                if self.zenithOld!=0:
                    #need to bin down to the correct size...
                    x=posDict[6]
                    x2=posDict[8]
                    #make it square (was rectangular):
                    gslCubSplineInterp(self.interpPhs[:posDict[3],:posDict[2]],x,x2,x,x,
                                       self.interpPhs[:posDict[3],:posDict[3]],self.interpolationNthreads)
                # print "atmos time4 %g"%(time.time()-t1)

                if self.sourceAlt>0:
                    #this is a LGS... need to project (interpolate/expand up to the full npup sized array.)
                    x2=posDict[6]
                    x=posDict[7]
                    #Bicubic interpolation (in C) for LGS projection
                    #print "gslCubS",self.interpolationNthreads,posDict[3]
                    #util.FITS.Write(self.interpPhs[:posDict[3],:posDict[3]],"tmp.fits")
                    #util.FITS.Write(x2,"tmp.fits",writeMode="a")
                    #util.FITS.Write(x,"tmp.fits",writeMode="a")
                    #util.FITS.Write(self.interpPhs,"tmp.fits",writeMode="a")
                    gslCubSplineInterp(self.interpPhs[:posDict[3],:posDict[3]],x2,x2,x,x,
                                       self.interpPhs,self.interpolationNthreads)
                    #print "gslCubs done"
                #print "atmos time5 %g"%(time.time()-t1)
                self.outputData+=self.interpPhs#and finally add the phase.
        if self.computeUplinkTT:
            if self.sourceLam!=500:
                self.uplinkPhs*=(500./self.sourceLam)
            self.uplinkTT=numpy.array([tottip,tottilt])
            #This value should be added to a lgs wfs.  Note, tip is in x direction, tilt is in y direction.
            #print "Uplink tt (m shift on-sky): %s"%str(self.uplinkTT)
            #So, how much zernike does this correspond to?
            self.outputData+=self.uplinkDownTT[0]*tottip
            self.outputData+=self.uplinkDownTT[1]*tottilt
        #print "atmos time6 %g"%(time.time()-t1)
        #from here to end takes a long time (1/4), equally spaced between each.
        # Remove overall piston
        if control["fullPupil"]:
            pfn=1
            pfnArea=self.npup*self.npup
        else:
            pfn=self.pupil.fn
            pfnArea=self.pupil.area
        if control["removePiston"]:
            pist=numpy.sum(self.outputData*pfn)/pfnArea
            self.outputData-=pist
        else:
            pist=0
        #print "atmos time7 %g"%(time.time()-t1)
        # Multiply by pupil and scale to output wavelength
        if not control["fullPupil"]:
            self.outputData*=pfn
        if self.sourceLam!=500.:
            self.outputData[:,]*=(500./self.sourceLam)
        #print "atmos time8 %g"%(time.time()-t1)
                   
    def createSingleLayerPhs(self,phaseScreens,interpPosColDict,interpPosRowDict,key,control):
        """Just do the interpolation for a single layer."""
        posDict=self.positionDict[key]
##         interpPosCol=interpPosColDict[key]+posDict[4]
##         interpPosRow=interpPosRowDict[key]+posDict[5]
        
##         #print "%g %g"%(interpPosCol,interpPosRow)
##         phsShiftY=int(interpPosCol)
##         phsShiftX=int(interpPosRow)
##         interpPosCol-=phsShiftY
##         interpPosRow-=phsShiftX
##         if self.colAdd[key]>=0:
##             interpPosCol=1-interpPosCol
##             #phsShiftY*=-1
##         if self.rowAdd[key]>=0:
##             interpPosRow=1-interpPosRow
##             #phsShiftX*=-1

        interpPosCol=posDict[4]#static offset due to source direction.
        interpPosRow=posDict[5]#static offset due to source direction.
        if self.colAdd[key]>=0:#add on to start...
            if interpPosColDict[key]>0:
                interpPosCol+=1-interpPosColDict[key]
        else:
            if interpPosColDict[key]>0:
                interpPosCol+=interpPosColDict[key]
            else:
                interpPosCol+=1

        phsShiftY=numpy.floor(interpPosCol)
        tmpcol=interpPosCol
        interpPosCol-=phsShiftY
        #if interpPosCol==0 and self.colAdd[key]<0:
        #    interpPosCol=1.


        if self.rowAdd[key]>=0:#add on to start
            if interpPosRowDict[key]>0:
                interpPosRow+=1-interpPosRowDict[key]
        else:
            if interpPosRowDict[key]>0:
                interpPosRow+=interpPosRowDict[key]
            else:
                interpPosRow+=1
        phsShiftX=numpy.floor(interpPosRow)
        tmprow=interpPosRow
        interpPosRow-=phsShiftX
        #if interpPosRow==0 and self.rowAdd[key]<0:
        #    interpPosRow=1.

        #phsShiftX+=1
        #phsShiftY+=1
        #now select the area we're interested in for this target...
        phs=phaseScreens[key][posDict[1]+phsShiftX:posDict[1]+posDict[3]+1+phsShiftX,
                              posDict[0]+phsShiftY:posDict[0]+posDict[2]+1+phsShiftY]
        #now do interpolation... (sub-pxl)
        linearshift(phs,interpPosCol,interpPosRow,self.interpPhs[:posDict[3],:posDict[2]])
        if self.zenithOld!=0:
            x=posDict[6]
            x2=posDict[8]
            #make it square (was rectangular):
            gslCubSplineInterp(self.interpPhs[:posDict[3],:posDict[2]],x,x2,x,x,
                               self.interpPhs[:posDict[3],:posDict[3]],self.interpolationNthreads)
            
        if self.sourceAlt>0:
            #this is a LGS... need to project (interpolate)
            x2=posDict[6]
            x=posDict[7]
            #Bicubic interpolation (in C) for LGS projection
            gslCubSplineInterp(self.interpPhs[:posDict[3],:posDict[3]],x2,x2,x,x,
                               self.interpPhs,self.interpolationNthreads)
        if control["fullPupil"]:
            return self.interpPhs
        else:
            return self.interpPhs*self.pupil.fn
