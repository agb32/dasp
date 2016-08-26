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

#import string,math
import numpy
#import util.dot as quick
import time,os
from util.tel import Pupil
import cmod.sor
import util.sor
from cmod.zernike import zern
#from LinearAlgebra import singular_value_decomposition
import numpy.linalg
import base.aobase
import util.zernikeMod,cmod.binimg,util.FITS


class recon(base.aobase.aobase):
    """Reconstructor for segmented DM.
    Calculate segment piston and tilt updates from WFS centroids using
    SOR or Zernike reconstruction.
    If more than one centroider is present as parent, will average the centroids together before doing the reconstruction.
    """
    
    def __init__(self,parent,config,args={},forGUISetup=0,debug=None,idstr=None):
        if type(parent)!=type({}):
            parent={1:parent}
        base.aobase.aobase.__init__(self,parent,config,args,forGUISetup=forGUISetup,debug=debug,idstr=idstr)
        self.fpDatatype=self.config.getVal("fpDatatype",default=numpy.float64)
        if forGUISetup==1:
            nsubx=self.config.getVal("wfs_nsubx")
            self.outputData=[(3,nsubx,nsubx),self.fpDatatype]
        else:
            self.control={"zero_dm":0,"close_dm":1}
            #self.generate=1
            self.dataValid=1#data is valid even before the first iteration because we assume the mirror starts zerod
            #self.newDataWaiting=1
            self.nParents=len(self.parent.keys())
            npup=self.npup=self.config.getVal("npup")
            self.timing=self.config.getVal("timing")
            wfs_nsubx=nsubx=self.wfs_nsubx=self.config.getVal("wfs_nsubx")
            wfs_minarea=self.config.getVal("wfs_minarea")
            n=wfs_n=self.wfs_n=self.config.getVal("wfs_n")
            self.zernike=self.config.getVal("zernike",default=0)
            self.maxiters=self.config.getVal("SORmaxiters",default=1000)
            self.soriters=0
            self.conv=self.config.getVal("SORconv",default=0.01)
            #self.makez=self.config.getVal("makez",default=0)
            self.sorpist=numpy.zeros((nsubx,nsubx),self.fpDatatype)                         # Arrays for segment pistons + tilts
            self.outputData=numpy.zeros((3,nsubx,nsubx),self.fpDatatype)
            self.pist=self.outputData[0]  # Numeric.zeros((nsubx,nsubx),Numeric.Float64)
            self.xtilt=self.outputData[1] # Numeric.zeros((nsubx,nsubx),Numeric.Float64)
            self.ytilt=self.outputData[2] # Numeric.zeros((nsubx,nsubx),Numeric.Float64)
            self.SRCDIR=self.config.getVal("srcdir",default=os.getcwd())#location to write temporary files.
            if self.SRCDIR[-1]!="/":
                self.SRCDIR+="/"
            #self.dmlist=self.config.getVal("dmlist")

            self.gamma=self.config.getVal("gamma")                  # The loop gain
            self.tiltSens=self.config.getVal("tiltSens",default=0)
            self.pupil=self.config.getVal("pupil")
            self.sormask=util.sor.createSorMask(self.pupil,wfs_n,wfs_minarea).ravel()
            self.avPistDivisor=numpy.sum((self.sormask>0).flat)

            self.pupsub=self.pupil.pupsub
            self.subarea=self.pupil.subarea
            self.subflag=self.pupil.subflag
            self.ndata=self.pupil.ndata#number of centroids that will be computed.
            self.subflaggamma=self.subflag*self.gamma
            self.subflaggamma2pi=self.subflaggamma*2*numpy.pi
            self.nsubUsed=self.ndata/2
            if (type(self.subflag.itemsize)==type(1) and self.subflag.itemsize==8) or (type(self.subflag.itemsize)==type(numpy.zeros) and self.subflag.itemsize()==8):
                print "WARNING: el_recon not tested with 8 byte longs (subflag)"

            self.centroids=None
            self.inputData=numpy.zeros((nsubx,nsubx,2),self.fpDatatype)
            #self.centx=None
            #self.centy=None
            self.centx=self.inputData[:,:,0]#Numeric.zeros((nsubx,nsubx),Numeric.Float64)                           # Centroid arrays
            self.centy=self.inputData[:,:,1]#Numeric.zeros((nsubx,nsubx),Numeric.Float64)
            self.tmparr=numpy.zeros((2,nsubx,nsubx),numpy.float32)
            self.raiseReconError=1#will raise an error if sor convergence not met.  This can be set to zero by the GUI during testing... ie to see what happens if you limit the number of iterations.
            self.sorStruct=cmod.sor.init(self.sormask,self.avPistDivisor,self.conv,self.maxiters,self.tmparr)

            if(self.zernike): # Setup Zernike fitting if reqd
                self.data=numpy.zeros(self.ndata,numpy.float64)
                self.recon_jmin=self.config.getVal("recon_jmin")
                self.control["nz"]=self.config.getVal("nz")#max number of zernikes to use... can be changed during simulation. should never be bigger than nzmax.
                self.nzrad=self.config.getVal("nzrad")#max radial order
                nzmax=self.nzmax=(self.nzrad+1)*(self.nzrad+2)/2#self.config.getVal("nzmax")#max number of zernikes during initialisation.
                self.control["zern_tilts"]=self.config.getVal("zern_tilts")#whether to use tilts during fitting zernikes.
                self.coeff=numpy.zeros((nzmax),numpy.float64)
                files=os.listdir(self.SRCDIR)
                files.sort()
                filefound=[]
                filestart="shz%d_%d_"%(npup,nsubx)
                for file in files:
                    if file[:len(filestart)]==filestart and file[-5:]==".fits":
                        try:
                            tmpnzmax=int(file[len(filestart):-5])
                            if tmpnzmax>=nzmax:
                                filefound.append(tmpnzmax)
                        except:
                            pass
                filefound.sort()
                
                if len(filefound)==0:#not os.path.exists(self.SRCDIR+"shz%d.fits"%nsubx):
                    # (Re)calc and save subap averaged Zern gradients
                    print "Computing zernike gradients"
                    self.subzpist,self.subzxgrad,self.subzygrad=self.zgrads(npup,nsubx,nzmax,numpy.where(self.subarea==0,0.,1./self.subarea))
                else:# Read zernike pist and grads from file
                    subgrad=util.FITS.Read(self.SRCDIR+"shz%d_%d_%d.fits"%(npup,nsubx,filefound[0]))[1]
                    self.subzpist=subgrad[0]#Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)
                    self.subzxgrad=subgrad[1]#Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)
                    self.subzygrad=subgrad[2]#Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)
                    #f=open(self.SRCDIR+'shz%d.dat'%nsubx,'r')
                    #for i in range(nsubx):
                    #    for j in range(nsubx):
                    #        for m in range(nzmax):
                    #            stuff=split(f.readline())
                    #            self.subzpist[i][j][m]=float(stuff[0])
                    #            self.subzxgrad[i][j][m]=float(stuff[1])
                    #            self.subzygrad[i][j][m]=float(stuff[2])
                    #f.close()
                k=-2
                subzgrad=numpy.zeros((self.ndata,nzmax),numpy.float64)
                for i in range(nsubx): # Zernike SVD fit setup
                    for j in range(nsubx):
                        if(self.subflag[i][j]==1):
                            k+=2
                            for m in range(nzmax):
                                subzgrad[k,m]=   self.subzxgrad[m,i,j]
                                subzgrad[k+1,m]= self.subzygrad[m,i,j]
                print "Performing zernike setup"
                #print "Doing svd recon"
                self.zernikeReconMx,self.eigenVals=self.calcZernikeRecon(subzgrad,1e-6)
                #print "done"
                # save some multiplication later on in a loop...
                self.subzpist*=self.wfs_nsubx*numpy.pi
                self.subzxgrad*=2.*numpy.pi
                self.subzygrad*=2.*numpy.pi
            self.niter=0

    def newSorStruct(self):
        """Frees existing sorStruct, and creates new one.  You should
        change the values that are to be changed before this is called
        (using the GUI).  This will be called using the GUI."""
        cmod.sor.free(self.sorStruct)
        self.sorStruct=cmod.sor.init(self.sormask,self.avPistDivisor,self.conv,self.maxiters,self.tmparr)

    def calcZernikeRecon(self,pokemx,corr_thresh=None):
        """Make control matrix from poke matrix
        corr_thresh is the fractional value of the threshold,
        ie the threshold will be placed at corr_thresh*max value."""
        svec=numpy.linalg.svd(pokemx)#singular_value_decomposition(pokemx)
        u,a,vt=svec

        #util.FITS.Write(a,'singular_values.fits')
        #a[self.n_modes_corrected:]=0.

        n_removed=0
        thresh=max(a)*corr_thresh
        for i in range(len(a)):
            if a[i] < thresh:
                a[i]=0.
                n_removed += 1
        print 'Removed',n_removed,'modes from control matrix'
        #print 'Eigenvalues:',a

        ut=numpy.transpose(u)
        v=numpy.transpose(vt)
        id=numpy.identity(len(a))
        ai=numpy.multiply(a,id)
        ai=numpy.where(ai != 0, 1/ai, 0)

        reconmx = quick.dot(v, quick.dot(ai, ut))
        return reconmx,a


    def generateNext(self,msg=None):    
        """Reconstructor main loop"""
        t1=time.time()
        if self.debug!=None:
            print "el_recon: Doing generateNext (debug=%s)"%str(self.debug)

        if self.generate==1:
            self.niter+=1
            if self.newDataWaiting:
                nin=0
                for key in self.parent.keys():
                    if self.parent[key].dataValid==1:
                        if nin==0:
                            self.inputData[:,]=self.parent[key].outputData#.copy() # Make a copy because we alter it
                        else:
                            self.inputData+=self.parent[key].outputData
                        nin+=1
                    else:
                        #print "el_recon: Waiting for data from cent, but not valid (iter %d)"%self.niter
                        pass
                if nin>0:
                    self.dataValid=1
                    self.inputData/=numpy.array(nin,self.inputData.dtype)
                    if nin!=self.nParents:
                        print "el_recon: Received unexpected number of centroid object results (ie not all parents)"
                else:
                    self.dataValid=0
            if self.dataValid:
                #self.centx=self.inputData[:,:,0]
                #self.centy=self.inputData[:,:,1]
                self.calc()                                                       # Calculate DM updates
        else:
            self.dataValid=0
        if self.debug!=None:
            print "el_recon: Done generateNext (debug=%s)"%str(self.debug)
        self.generateNextTime=time.time()-t1



    def calc(self):
        """Reconstruct piston errors for current tilt data"""
    
        nsubx=self.wfs_nsubx
        ndata=self.ndata  
        subflag=self.subflag
        if self.control["close_dm"]:
            if(self.zernike):  # Zernike reconstruction
                #print 'Zernike...'
                #first reorganise the centroids...
                idata=-1
                for i in range(nsubx):
                    for j in range(nsubx):
                        if(subflag[i,j]==1):
                            idata+=1
                            self.data[idata]=self.centx[i,j]
                            idata+=1
                            self.data[idata]=self.centy[i,j]
                #and then fit them.
                #t1=time.time()
                self.coeff[:,]=quick.dot(self.zernikeReconMx,self.data)
                #t2=time.time()
                #self.pist*=0.#now done by zero_dm...
                #self.xtilt*=0.
                #self.ytilt*=0.
                for i in range(nsubx):                                            # Pistons + tilt updates for DM:
                    for j in range(nsubx):
                        if(subflag[i,j]==1):
                            #self.pist[i][j]=0.
                            #self.xtilt[i][j]=0.
                            #self.ytilt[i][j]=0.
                            for m in range(self.recon_jmin,min(self.control["nz"]+1,self.nzmax)):
                                self.pist[i,j]+=self.subzpist[m,i,j]*self.coeff[m]#*(float(self.wfs_nsubx)*math.pi)-now done at start...
                                if(self.control["zern_tilts"]):    # Piston + tilts, or
                                    self.xtilt[i,j]+=self.subzxgrad[m,i,j]*self.coeff[m]#*2.*math.pi-now done at start
                                    self.ytilt[i,j]+=self.subzygrad[m,i,j]*self.coeff[m]#*2.*math.pi-now done at start
                            if not self.control["zern_tilts"]:  # Piston only, Tilts direct
                                self.xtilt[i,j]+= self.centx[i,j]*2.*numpy.pi
                                self.ytilt[i,j]+= self.centy[i,j]*2.*numpy.pi


            else:                    # SOR reconstruction
                if(self.tiltSens):   # Remove mean tilt - for Zern recon
                    #avcentx=0.       # this is done by excluding j=2,3
                    #avcenty=0.
                    #nflag=0.
                    #for i in range(nsubx):
                    #    for j in range(nsubx):
                    #        if(subflag[i][j]==1):
                    #            avcentx=avcentx+self.centx[i,j]               
                    #            avcenty=avcenty+self.centy[i,j]
                    #            nflag=nflag+1.
                    #self.centx=self.centx-avcentx/nflag
                    #self.centy=self.centy-avcenty/nflag
                    #this does the same as the last section (with avcentx etc)
                    #remove mean tilt: for zern this is done by excluding j=2,3
                    self.centx-=numpy.sum(self.centx*self.subflag)/self.nsubUsed
                    self.centy-=numpy.sum(self.centy*self.subflag)/self.nsubUsed

                self.sorpist[:,]=0.0
                self.soriters,err=cmod.sor.fit(self.centx,self.centy,self.sorpist,self.sorStruct)#mask,self.avPistDivisor,self.conv,self.maxiters,self.tmparr)
                if self.soriters==self.maxiters and self.raiseReconError:
                    print "WARNING - SOR didn't converge after %d iterations (err=%g, conv=%g).  Trying again with a larger convergence factor."%(self.soriters,err,self.conv)
                    cmod.sor.setconv(self.sorStruct,self.conv*10)
                    self.soriters,err=cmod.sor.fit(self.centx,self.centy,self.sorpist,self.sorStruct)#mask,self.avPistDivisor,self.conv*10,self.maxiters,self.tmparr)
                    cmod.sor.setconv(self.sorStruct,self.conv)
                    if self.soriters==self.maxiters:
                        print "Didn't converge this time either.  Giving up (err=%g)."%err
                        util.FITS.Write(self.centx,"sordebugx.fits")
                        util.FITS.Write(self.centy,"sordebugy.fits")
                        print max(self.centx.copy().flat)
                        print max(self.centy.copy().flat)
                        raise Exception("sor module didn't converge. Centroids written to sordebugx.fits and sordebugy.fits.  If this happens and you are running the wfscent FPGA module on node1, try running on a different node - possibly there is a problem with the FPGA...")
                if self.debug!=None:
                    print "%s: done sor fit in %d iters"%(str(self.debug),self.soriters)
                #for i in range(nsubx):                                            # Pistons + tilt updates for DM:
                #    for j in range(nsubx):
                #        if(subflag[i][j]==1):
                #            self.pist[i][j]=  -self.sorpist[i][j]
                #            self.xtilt[i][j]= self.centx[i][j]*2.*math.pi
                #            self.ytilt[i][j]= self.centy[i][j]*2.*math.pi
                #this does the same as the previous for loop... (I think).
                self.pist[:,]-=self.sorpist*self.subflaggamma
                self.xtilt[:,]+=self.centx*self.subflaggamma2pi#2*math.pi*subflag*self.gamma
                self.ytilt[:,]+=self.centy*self.subflaggamma2pi#2*math.pi*subflag*self.gamma
    ##             self.pist[:,]=-self.sorpist*subflag
    ##             self.xtilt[:,]=self.centx*2*math.pi*subflag
    ##             self.ytilt[:,]=self.centy*2*math.pi*subflag
        if self.control["zero_dm"]:
            print "el_recon: Zeroing DM"
            self.control["zero_dm"]=0
            self.pist*=0
            self.xtilt*=0
            self.ytilt*=0
                        
    def zgrads(self,npup,nsubx,nzmax,invarea):
        """Get mean zernike pist and grads over subaps + write to file
        invarea is 1/subarea, except ==0 where subarea==0.
        """
        n=npup/nsubx
        j=0                                                               # Make the Zernike fns 
        pp=[0]
        qq=[0]
        tt=[1]
        trig=0
        for p in range(1,self.nzrad+1):
            for q in range(p+1):
                if(numpy.fmod(p-q,2)==0):
                    if(q>0):
                        pp.append(p)
                        qq.append(q)
                        trig=not(trig)
                        tt.append(trig)
                        pp.append(p)
                        qq.append(q)
                        trig=not(trig)
                        tt.append(trig)
                    else:
                        pp.append(p)
                        qq.append(q)
                        tt.append(1)
                        trig=not(trig)
        #nzin=55
        #f1=open(self.SRCDIR+'gamx.dat','r')# Read in the gamma matrices - can be found in aosim/data/
        #stuff1=split(f1.read())
        #f1.close()    
        #f2=open(self.SRCDIR+'gamy.dat','r')    
        #stuff2=split(f2.read())
        #f2.close()
        #gamx=Numeric.zeros((nzin,nzin),Numeric.Float64)
        #gamy=Numeric.zeros((nzin,nzin),Numeric.Float64)
        #k=-1
        #for i in range(nzin):
        #    for j in range(nzin):
        #        k=k+1
        #        gamx[i][j] = float(stuff1[k])
        #        gamy[i][j] = float(stuff2[k])
        gamx=None
        gamy=None
        #look for suitable gamma file...
        files=os.listdir(self.SRCDIR)
        filex=[]
        filey=[]
        for file in files:
            if file[:4]=="gamx" and file[-5:]==".fits":
                try:
                    nzradtmp=int(file[4:-5])
                    if nzradtmp>=nzrad:
                        filex.append(nzradtmp)
                except:
                    pass
            if file[:4]=="gamy" and file[-5:]==".fits":
                try:
                    nzradtmp=int(file[4:-5])
                    if nzradtmp>=nzrad:
                        filey.append(nzradtmp)
                except:
                    pass
        #if os.path.exists(self.SRCDIR+"gamx.fits"):#found in aosim/data/
        #    gamx=util.FITS.Read(self.SRCDIR+"gamx.fits")[1]
        #if os.path.exists(self.SRCDIR+"gamy.fits"):
        #    gamy=util.FITS.Read(self.SRCDIR+"gamy.fits")[1]
        #if type(gamx)==type(None) or type(gamy)==type(None) or gamx.shape[0]<self.nzrad or gamy.shape[0]<self.nzrad:
        #    gamx,gamy=util.zernikeMod.makegammas(self.nzrad,returntype="Numeric")
        if len(filex)==0 or len(filey)==0:
            gamx,gamy=util.zernikeMod.makegammas(self.nzrad,returntype="numpy")
        else:
            filex.sort()
            filey.sort()
            gamx=util.FITS.Read(self.SRCDIR+"gamx%d.fits"%filex[0])[1]
            gamy=util.FITS.REad(self.SRCDIR+"gamy%d.fits"%filey[0])[1]
        #zxgrad=Numeric.zeros((npup,npup,nzmax),Numeric.Float64)           # Make Zernike derivative fns
        #zygrad=Numeric.zeros((npup,npup,nzmax),Numeric.Float64)
        #zxgrad=Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)           # Make Zernike derivative fns
        #zygrad=Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)

        subgrad=numpy.zeros((3,nzmax,nsubx,nsubx),numpy.float64)
        subzpist=subgrad[0]#Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)       # write to file
        subzxgrad=subgrad[1]#Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)      # Get mean zernike pist
        subzygrad=subgrad[2]#Numeric.zeros((nzmax,nsubx,nsubx),Numeric.Float64)      # +grads over subaps and

        #f=open(self.SRCDIR+'shz%d.dat'%nsubx,'w')

        #z=Numeric.zeros((npup,npup,nzmax),Numeric.Float64)
        ztmp=numpy.zeros((npup,npup),numpy.float64)
        #for j in range(nzmax):
        #    zern(ztmp,pp[j],qq[j],tt[j])
        #    z[:,:,j]=ztmp
        tmp=numpy.zeros((nsubx,nsubx),numpy.Float64)
            
        for m in range(nzmax):
            print "Creating zern %d/%d"%(m,nzmax-1)
            zern(ztmp,pp[m],qq[m],tt[m])
            cmod.binimg.binimg(ztmp,tmp)#bin the zernike...
            subzpist[m]=tmp*invarea#divide by the area...
            #for i in range(nsubx):
            #    for j in range(nsubx):
            #        if(area[i][j]>0.):
            #            subzpist[m,i,j]=Numeric.sum(Numeric.sum(ztmp[i*n:(i+1)*n,j*n:(j+1)*n]))/area[i][j]
            for i in range(nzmax):
                cmod.binimg.binimg(gamx[i,m]*ztmp,tmp)
                subzxgrad[i]+=tmp
                cmod.binimg.binimg(gamy[i,m]*ztmp,tmp)
                subzygrad[i]+=tmp
        for m in range(nzmax):
            subzxgrad[m]*=invarea
            subzygrad[m]*=invarea
            #for i in range(nsubx):
            #    for j in range(nsubx):
            #        if(area[i,j]>0.):
            #            subzxgrad[m,i,j]/=area[i,j]
            #            subzygrad[m,i,j]/=area[i,j]

                    #f.write(`subzpist[i,j,m]`+' '+`subzxgrad[i][j][m]`+' '+`subzygrad[i][j][m]`+'\n')
        #f.close()
        util.FITS.Write(subgrad,self.SRCDIR+'shz%d_%d_%d.fits'%(npup,nsubx,nzmax))
        return subzpist,subzxgrad,subzygrad




if __name__=="__main__":
    recon=Recon()

