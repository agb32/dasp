"""A module to compute DM actuators from reconstructed tomographic
wavefronts, which are for virtual DMs.  IE this is used as an extra
reconstructor, taking reconstructor outputs, and compressing along a
given line of sight, to pass them into a DM model.

This basically takes the virtual DM actuators, sums (compresses) them
along a single line of sight and then outputs actuators for a true
(eg xinterp_dm) DM.

This true DM can have a non-zero field of view - in which case a
larger area of higher vDMs are used, and then squashed down onto the
true DM here.  The conjugate height of this true DM still doesn't
matter, since we are just squashing other layers onto it (or
expanding, if they are below its height).  This then gives a slightly
over-sized real DM, which can then be used as normal for different
source directions.

However, the height does matter once the shape of the DM has been
computed... this is then used to select which parts of the DM are
visible in which source directions (if this dm is not ground
conjugated).

the MOAO buttons have a non-zero field of view.  So, what we have to
do here is use the reconstructor data to create the virtual DM
actuators.  Then for each DM direction (note, there could be many
source directions on a single DM - we're concerned here with the
centre of the DM) we interpolate and compress the virtual DMs.  This
then gives actuator values for the true MOAO DM.  These are used with
the interpolation routine, to create the DM surface.  Then, for every
direction viewed by this DM, the appropriate part of this surface is
selected for the outputData.

Field of view: If the moao mirror is ground conjugate, the field of
view doesn't matter since all paths use all the mirror.  However, if
conjugated above ground, different parts of the mirror should be
selected depending on source direction.  The field of view defines how
large the DM should be.

Compressing along the line of sight - we can do two things here.
Either just select the npup portion of the virtual DMs.  Or use a
larger portion of virtual DMs with increasing height, such that we
average over the field of view.  This may give better off-axis
performance.

So, we need two fov factors.  One specifies the actual FOV and hence
size of the DM (which if ground conjugate is irrelevant).  The other
specifies the diameter of vDM segments to be used.  The two FOV values
are independant of each other.

This chould be a resource sharing object, sharing for one single real
DM object (multiple source directions on this moao dm).

Here, we assume the source is at infinity.

"""

#from science.xinterp_dm import interpolateSpline,interpolateBicubic
import cmod.interp
import numpy
import time
import base.aobase
import util.FITS
class vdmUser(base.aobase.aobase):
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.dmObj=self.config.getVal("dmOverview",raiseerror=0)
        if self.dmObj==None or type(self.dmObj)!=type(self):
            print "DEPRECATION: warning: dmObj should now be dmOverview"
            self.dmObj=self.config.getVal("dmObj")
        self.dmlabel=self.config.getVal("dmLabel",default=self.idstr[0],raiseerror=0)
        self.interpolationNthreads=self.config.getVal("interpolationNthreads",default=0)
        if self.dmlabel==None:
            print "Warning: science.vdmUser - dmlabel not found, using vdm.  idstr=%s"%str(self.idstr)
            self.dmlabel="vdm"
        self.dmInfo=self.dmObj.getDMFromLabel(self.dmlabel)
        self.pupil=self.config.getVal("pupil")
        self.dmflag=self.dmInfo.computeDMPupil(self.dmObj.atmosGeom,centObscuration=self.pupil.r2,retPupil=0)[0]
        self.nacts=numpy.sum(numpy.sum(self.dmflag))
        if forGUISetup==1:
            self.outputData=[(self.nacts,),numpy.float32]
        else:
            self.projMxFilename=self.config.getVal("projMxFilename",raiseerror=0)
            self.nact=self.dmInfo.nact
            self.actoffset=self.dmInfo.actoffset
            self.telDiam=self.config.getVal("telDiam")
            self.control={"dm_update":1}
            self.outputData=numpy.zeros((self.nacts),numpy.float32)
            self.actmap=numpy.zeros((self.nact,self.nact),numpy.float32)
            self.atmosGeom=self.config.getVal("atmosGeom")
            #if the reconstructor (tomoRecon) has an idstr, this should be used here... 
            # otherwise it won't agree... however, by default this is None, so should be ok:
            self.reconIdStr=self.config.getVal("reconIdStr",raiseerror=0)
            #the DMs reconstructed by the reconstructor.This is the order that the reconstructor will have them too.
            self.reconstructorDmList=self.dmObj.makeDMList(self.reconIdStr)
            #the DMs used here as virtual ones:
            virtDmListUnordered=self.dmObj.makeDMList(actsFrom=self.dmInfo.actuatorsFrom)
            if self.dmInfo in virtDmListUnordered:
                virtDmListUnordered.remove(self.dmInfo)
            if self.projMxFilename:
                nactstot=0
                for dm in self.reconstructorDmList:
                    if dm.zonalDM:
                        #dmindices=numpy.nonzero(dm.computeDMPupil(
                        #self.dmObj.atmosGeom,centObscuration=self.pupil.r2,retPupil=0)[0].ravel())[0]
                        #nacts=int(dmindices.shape[0])
                        tmp=dm.computeDMPupil(self.dmObj.atmosGeom,centObscuration=self.pupil.r2,retPupil=0)
                        nacts=int(numpy.sum(tmp[0].ravel()))

                    else:
                        nacts=dm.nact
                    print "nacts %s: %d"%(str(dm),nacts)
                    nactstot+=nacts
                self.projmx=numpy.zeros((self.nacts,nactstot),numpy.float32)
            else:
                self.projmx=None
            self.npup=self.config.getVal("npup")
            #the actuator offset from the edge of the mirror:
            self.offset=self.dmInfo.actoffset#config.getVal("offsetmdm")
            #get the direction of the centre of the DM, and the DM FOV.
            #the direction in which the moao mirror is facing:
            self.mdmTheta=self.dmInfo.primaryTheta#self.config.getVal("mdmTheta")
            #the direction in which the moao mirror is facing:
            self.mdmPhi=self.dmInfo.primaryPhi#config.getVal("mdmPhi")
            self.fov=self.dmInfo.fov#config.getVal("fov")#mirror field of view (irrelevant if ground conjugate).
            print "MOAO vdmUser is in direction theta=%g, phi=%g, with fov=%g"%(self.mdmTheta,self.mdmPhi,self.fov)
            self.compressFov=self.config.getVal("compressFov",default=0.)#specifies the amount of each vdm that is 
                                                # used.  If this is zero, only a telDiam sized section of the 
                                                # vdm is used.  If non-zero, means that a slightly larger 
                                                # section will be compressed onto the mdm.
            self.nactsList=[]
            self.nactsCumList=[]
            self.actValList=[]
            self.vdmCoordList=[]
            self.mdmActPosList=[]
            self.vdmActPosList=[]
            self.interpolated=numpy.zeros((self.nact,self.nact),numpy.float32)
            self.dmindices=numpy.nonzero(self.dmflag.ravel())[0]
            self.dmindicesList=[]
            nactstot=0
            pos=0
            self.virtDmList=[]
            if len(self.reconstructorDmList)==0:
                raise Exception("Unable to compute reconstructorDmList in vdmUser.py")
            print "vdm reconstructor dm list: %s"%str(self.reconstructorDmList)
            for dm in self.reconstructorDmList:
                if dm.zonalDM:
                    dmindices=numpy.nonzero(dm.computeDMPupil(self.dmObj.atmosGeom,
                                    centObscuration=self.pupil.r2,retPupil=0)[0].ravel())[0]
                    nacts=int(dmindices.shape[0])
                else:
                    nacts=dm.nact
                nactstot+=nacts
                if dm in virtDmListUnordered:
                    self.virtDmList.append(dm)
                    self.dmindicesList.append(dmindices)
                    self.nactsList.append(nacts)
                    self.nactsCumList.append((nactstot-nacts,nactstot))
                    #self.closedLoopList.append(self.dmInfo.closedLoop)
                    self.actValList.append(numpy.zeros((dm.nact,dm.nact),numpy.float32))
                    #find what section of the vdm is required for this MOAO DM.
                    widthvdm=2*numpy.tan(dm.fov/60./60./180.*numpy.pi)*dm.height+self.telDiam
                    widthmdm=2*numpy.tan(self.compressFov/60./60./180.*numpy.pi)*dm.height+self.telDiam
                    xc=numpy.tan(self.mdmTheta/60./60./180.*numpy.pi)*dm.height*numpy.cos(self.mdmPhi/180.*numpy.pi)
                    yc=numpy.tan(self.mdmTheta/60./60./180.*numpy.pi)*dm.height*numpy.sin(self.mdmPhi/180.*numpy.pi)
                    xmin=xc-widthmdm/2
                    ymin=yc-widthmdm/2
                    xmax=xc+widthmdm/2
                    ymax=yc+widthmdm/2
                    tol=1e-10#add a timy bit of tolerance to avoid floating point rounding/precision errors.
                    if xmin+tol<-widthvdm/2 or ymin+tol<-widthvdm/2 or xmax-tol>widthvdm/2 or ymax-tol>widthvdm/2:
                        print xmin,ymin,xmax,ymax,widthvdm,dm.height
                        raise Exception("Error: vdmUser.py - mdm will not fit in vdm")
                    actSpacingMdm=widthmdm/(self.nact-1+2*self.actoffset)
                    actSpacingVdm=widthvdm/(dm.nact-1+2*dm.actoffset)
                    #onsky positions of the actuators of the MOAO DM at this layer height.
                    actPosMdmx=(numpy.arange(self.nact)-(self.nact-1)/2.)*actSpacingMdm+xc
                    actPosMdmy=(numpy.arange(self.nact)-(self.nact-1)/2.)*actSpacingMdm+yc
                    self.mdmActPosList.append((actPosMdmy,actPosMdmx))
                    #onsky positions of the actuators of the virtual DM (same in x and y since symmetrical).
                    actPosVdm=(numpy.arange(dm.nact)-(dm.nact-1)/2.)*actSpacingVdm
                    #Now select and store the portions of this which contain the mdm actuators.
                    oe=((dm.nact-1)%2)/2.#odd or even?
                    #sign=int(actPosMdmx[0]>=0)*2-1#get the sign.
                    sx=numpy.floor(actPosMdmx[0]/actSpacingVdm+oe+tol)-oe
                    ex=numpy.ceil(actPosMdmx[-1]/actSpacingVdm+oe-tol)-oe
                    sy=numpy.floor(actPosMdmy[0]/actSpacingVdm+oe+tol)-oe
                    ey=numpy.ceil(actPosMdmy[-1]/actSpacingVdm+oe-tol)-oe

                    #now get the positions of the relevant vdm actuators (in m).These are then used for interpolation.
                    actPos=(numpy.arange(sy,ey+1)*actSpacingVdm,numpy.arange(sx,ex+1)*actSpacingVdm)
                    self.vdmActPosList.append(actPos)
                    #and get the array coordinates for these actuator values.
                    off=(dm.nact-1)/2.
                    coords=(int(sy+off),int(ey+off)+1,int(sx+off),int(ex+off)+1)
                    self.vdmCoordList.append(coords)
                    if self.projmx!=None:#shape=nacts,totnacts
                        actVal=self.actValList[-1].ravel()
                        for i in dmindices:
                            print "Computing projection matrix column %d/%d  \r"%(pos,self.projmx.shape[1]),
                            try:
                                actVal[i]=1
                                res=numpy.take(self.compress(len(self.actValList)-1).ravel(),self.dmindices)
                                self.projmx[:,pos]=res
                                pos+=1
                                actVal[i]=0
                            except:
                                print "(%d %d proj shape %s, interpolated shape %s, resshape %s)"%(
                                    pos,self.projmx.shape[1],str(self.projmx.shape),
                                    str(self.interpolated.shape),str(res.shape))
                                raise
                        print "Projection matrix done next mirror... %s %s %s                     "%(str(self.interpolated.shape),str(self.projmx.shape),str(res.shape))
                    #     tmp=numpy.zeros((dm.nact,dm.nact),numpy.float32)
                    #     tmp2=numpy.zeros((dm.nact,dm.nact),numpy.float32)
                    #     tmp=tmp[coords[0]:coords[1],coords[2]:coords[3]]#select the 
                    #         portion that impactsthis direction.
                    #     #now work out which vdm actuators affect dm.
                    #     yin=actPos[0]
                    #     xin=actPos[1]
                    #     st=nactsCumList[-1][0]
                    #     en=nactsCumList[-1][1]
                    #     #interpolate
                    #     cmod.interp.gslCubSplineInterp(tmp,yin,xin,actPosMdmy,actPosMdmx,tmp2,4)
                    #     #and then for each dm actuator, compute how it is affected by neighbours.
                    #     #compute how actuators should be selected in x and y directions.
                    #     xinterp=numpy.interp(numpy.arange(dm.nact),numpy.arange(actPosMdmx.size,actPosMdmx))
                    #     yinterp=xxx
                    #     ipos=0
                    #     for i in range(dm.nact):
                    #         for j in range(dm.nact):
                    #             if dm.dmflag[i,j]:
                    #                 self.projmx[pos,c1]=xinterp[j]*yinterp[i]
                    #                 self.projmx[pos,c2]=(1-xinterp[j])*yinterp[i]
                    #                 self.projmx[pos,c3]=xinterp[j]*(1-yinterp[i])
                    #                 self.projmx[pos,c4]=(1-xinterp[j])*(1-yinterp[i])
                    #                 pos+=1
                
                        
            self.nactsInput=sum(self.nactsList)
            if self.projmx!=None:
                util.FITS.Write(self.projmx,self.projMxFilename)
                util.FITS.Write(self.dmflag,self.projMxFilename,writeMode="a")
                print "Saved projection matrix  as %s (shape %s)"%(self.projMxFilename,str(self.projmx.shape))
    def generateNext(self):
        """Compression main loop - compress along a line of sight to get DM actuator values."""
        if self.generate==1:
            if self.newDataWaiting:
                if self.parent.dataValid==1:
                    self.dataValid=1
                    if self.control["dm_update"]==1 and self.currentIdObjCnt==0:
                        self.reconData=self.parent.outputData
                        self.update()
                else:
                    print "vdmUser waiting for data from parent, but not valid"
                    self.dataValid=0
        else:
            self.dataValid=0

    def update(self):
        """put the reconstructor data into arrays, interpolate these and compress to get the actual DM shape.
        """
        #reset the actual actuator map.
        #self.actmap[:,]=0.
        for i in range(len(self.virtDmList)):
            self.getInputActuators(i)
        #Now compress all the DM actuators along the line of sight.
        self.actmap[:,]=0#reset the mirror...
        for i in range(len(self.virtDmList)):
            #find the coords of this vDM that we need to work with.
            #Coords of -1 and 1 specify the very edges of the vDM.
            #Depending on offset, there may be an actuator at -1 and 1 (offset=0) or maybe not.
            self.actmap+=self.compress(i)#self.interpolated#add the actuator values to the mirror.
        #now take the required actuators to the output...
        #self.doInterpolation(self.actmap,self.outputData)
        self.outputData[:,]=numpy.take(self.actmap.ravel(),self.dmindices)

    def compress(self,indx):
        """Compress actuators... of vdm index indx."""
        c=self.vdmCoordList[indx]
        vdm=self.actValList[indx][c[0]:c[1],c[2]:c[3]]#the virtual DM actuator values that are in the fov.
        yin=self.vdmActPosList[indx][0]
        xin=self.vdmActPosList[indx][1]
        #this is just a sub-pxl interpolation I think...(or a compression for LGS)
        # - not interpolation of the DM surface...
        # Number of threads set to 1, since the function is not used in aosim
        cmod.interp.gslCubSplineInterp(vdm,yin,xin,self.mdmActPosList[indx][0],self.mdmActPosList[indx][1],
                                       self.interpolated,self.interpolationNthreads)
        return self.interpolated

    def getInputActuators(self,dmno):
        """Put the 1D array of actuators for the given dm no into a 2D array"""
        i=dmno
        dm=self.virtDmList[i]
        numpy.put(self.actValList[i].ravel(),self.dmindicesList[i],
                  self.reconData[self.nactsCumList[i][0]:self.nactsCumList[i][1]])
        return self.actValList[i]

        
    def plottable(self,objname="$OBJ"):
        """Return a XML string which contains the commands to be sent
        over a socket to obtain certain data for plotting.  The $OBJ symbol
        will be replaced by the instance name of the object - e.g.
        if scrn=mkscrns.Mkscrns(...) then $OBJ would be replaced by scrn."""
        if self.idstr[0]==None or self.idstr[0]=="":
            id=""
        else:
            id=" (%s)"%self.idstr[0]
        txt="""<plot title="vdmUser Output data%s" cmd="data=%s.outputData" ret="data" type="pylab" when="rpt" palette="gray"/>\n"""%(id,objname)
        for i in range(len(self.virtDmList)):
            txt+="""<plot title="vdmUser input actuators %d%s" cmd="data=-%s.getInputActuators(%d)" ret="data" when="rpt" type="pylab"/>\n"""%(i,id,objname,i)
        txt+="""<plot title="vdmUser output actuators%s" cmd="data=-%s.actmap" ret="data" when="rpt" type="pylab"/>\n"""%(id,objname)
        for i in range(len(self.virtDmList)):
            txt+="""<plot title="vdmUser compressed actators %d%s" cmd="data=-%s.compress(%d)" ret="data" when="rpt" type="pylab"/>\n"""%(i,id,objname,i)
        return txt
	
