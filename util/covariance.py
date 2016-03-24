
#
# Modified fitting error program
#  - uses csim instead of old phase screen/DM code
#  

"""Computes simplexes of DM to phase, and writes the actuators to a file.  Continues doing this for many realisations of the phase until terminated.
"""

import sys
import time

#import Numeric
import numpy
#import pylab

import amoeba
#import AO
import util.FITS
import util.dm
import science.infScrn
#import turbulence
#import util.dm
import base.readConfig


# DM wrapper function to cope with unused actuators
def DMsurface(DM,flags,actvec):
	assert flags.sum()==actvec.shape[0]
	n=flags.shape[0]
	actmap=numpy.zeros((n,n),numpy.float)
	ii=0
	for i in range(n):
		for j in range(n):
			if flags[i,j]>0:
				actmap[i,j]=actvec[ii]
				ii+=1
	return DM.fit(actmap)


# Pupil function
def pupil(npup,r1=-1.,r2=-1.):
    if r1 < 0.:
        r1 = float(npup/2)
    pup=numpy.zeros((npup,npup),numpy.int)
    for i in range(npup):
        for j in range(npup):
            r=numpy.sqrt((float(i-npup/2)+0.5)**2+(float(j-npup/2)+0.5)**2)
            #r=sqrt((float(i-npup/2))**2+(float(j-npup/2))**2)
            if ((r<=r1)&(r>=r2)):
                pup[i,j]=1
    return pup



class FitDM(amoeba.amoeba):
	def initialise(self,nact,npup,pupil,DM,phs,dmflags,subtractPiston):
		self.nact=nact
		self.npup=npup
		self.pupil=pupil
		self.phs=phs
		self.DM=DM
		self.dmflags=dmflags
		self.subtractPiston=subtractPiston

	def fun(self,soln,final=0): #overwrite Nirmal's merit function
		self.funEvals+=1

		dmphs=DMsurface(self.DM,self.dmflags,soln)

		area=float(self.pupil.sum())
		#print "shapes",dmphs.shape,self.phs.shape,self.pupil.shape
		diff=(dmphs-self.phs)*self.pupil
		if self.subtractPiston:
			pist=diff.sum()/area
			diff-=pist*self.pupil
		rms=numpy.sqrt((diff**2).sum()/area)
		print "RMS/nm %g"%(500.*rms/numpy.pi/2.)
		return rms

	def showme(self):
		dmphs=DMsurface(self.DM,self.dmflags,self.bestSoln)
		print phs.min(),phs.max(),dmphs.min(),dmphs.max()
		

diii=int(sys.argv[1])#step
iii=eval(sys.argv[2])#start
paramfile=sys.argv[3]
try:
	paramfile=eval(paramfile)
except:
	pass
paramfile=paramfile.split(",")
batchno=int(sys.argv[4])
dmid=sys.argv[5]
outdir=sys.argv[6]#'cov3'
subtractPiston=int(sys.argv[7])#probably should be set.
writeFailed=0
if len(sys.argv)>8:
	writeFailed=int(sys.argv[8])#write the failed ones?
config=base.readConfig.AOXml(paramfile,batchno=batchno)
atmosGeom=config.getVal("atmosGeom")
dmObj=config.getVal("dmOverview",raiseerror=0)
if dmObj==None or type(dmObj)!=type(atmosGeom):
        print "DEPRECATION: warning: dmObj should now be dmOverview"
        dmObj=config.getVal("dmObj")
pupil=config.getVal("pupil")
dmflag,subarea,dmpupil=dmObj.computeDMPupil(dmid,pupil.r2,1)
#wfs_n=config.getVal("wfs_n")

#dmflag=util.FITS.Read('dmflag.fits')[1].astype(numpy.int)
nact=dmflag.shape[0]
#nsubx=nact-1
#wfs_n=12
npup=dmpupil.shape[0]
dmDiam=dmObj.getDM(dmid).dmDiam
#npup=(nact-1)*wfs_n
#if npup!=config.getVal("npup"):
#	raise Exception("npup != npup %d %d"%(npup,config.getVal("npup")))
totnact=dmflag.sum()

#telDiam=config.getVal("telDiam")#10
r0=atmosGeom.r0
#r0=0.106
#npup=160 #128 #64
#nact=21 #17 #9
#r0pix=wfs_n*r0/(telDiam/nsubx)  #8. #16. #0.12
#lam=500.e-9
#nscrn=1024 #512

print 'nact =',nact
print 'npup =',npup
#print 'nscrn/npup =',float(nscrn)/float(npup)
print 'totnact =',totnact
print "dmpupil.shape=",dmpupil.shape
print "dmDiam",dmDiam
#D=4.2
#D2=1.2

SIMPLEX_TOLERANCE=0.001#1e-4#1e-6
SIMPLEX_MAX_ITER=100000
#pup=pupil(npup, float(npup/2), float(npup/2)/7.)
pup=dmpupil#pupil.fn
#pup=pupil(npup,float(npup/2)) #D2/D*float(npup/2))
#pup=numpy.ones((npup,npup),numpy.int)

#DM=AO.dm.CubicSpline(nact,npup)
DM=dmObj.getDM(dmid).getMirrorSurface()
#DM=util.dm.MirrorSurface("pspline",npup,nact)

#turb=AO.atmos.KolScrn(nscrn,r0pix) #float(npup)*r0/D)
strLayer=None
dmheight=dmObj.getDM(dmid).height
for key in atmosGeom.layerDict.keys():
	l=atmosGeom.layerDict[key]
	if l.height==dmheight:#DM conjugate to this layer.
		strLayer=l.strength
		break
warnStr=0
if strLayer==None:
	strLayer=1.
	warnStr=1

print "Layer strength:",strLayer 
while 1:
	#turb.refresh()
	t0=time.time()
	scrn=science.infScrn.makeInitialScreen(dpix=npup,Dtel=dmDiam,L0=atmosGeom.l0,scrnXPxls=None,scrnYPxls=None,seed=int(time.time()*100),tstep=0.05,globR0=r0,strLayer=strLayer,windDirection=190.,vWind=10.)
	if warnStr:
		print "WARNING - Layer not found conjugate to this DM - using layer strength of 1"
	#note - windDirection and vWind aren't used at all, except to put the phase in the right part of scrn (leaving blank rows/cols where they won't be used).
	t1=time.time()
	phs=scrn[:npup,:npup]*pup
	if subtractPiston:
		pist=phs.sum()/float(pup.sum())#npup*npup)
		phs-=pist
		phs*=pup

	var=(phs**2).sum()/float(pup.sum())
	print 'Phase variance/nm = ',var*500/numpy.pi/2.
	
	#simplex=FitDM(nact*nact,SIMPLEX_TOLERANCE,SIMPLEX_MAX_ITER,SIMPLEX_MAX_ITER)
	simplex=FitDM(totnact,SIMPLEX_TOLERANCE,SIMPLEX_MAX_ITER,SIMPLEX_MAX_ITER)
	simplex.initialise(nact,npup,pup,DM,phs,dmflag,subtractPiston)
	#simplex.firstGuess(numpy.zeros(nact*nact,numpy.float),30.)
	indx=(numpy.arange(nact)*(npup-1.)/(nact-1.)).astype("i")
	firstGuess=scrn[indx][:,indx].ravel()
	firstGuess=numpy.take(firstGuess,numpy.nonzero(dmflag.ravel())[0])
	t2=time.time()
	simplex.firstGuess(firstGuess,1.)
	t3=time.time()
	success=simplex.run()
	if success==0:
		print "Failure to converge."
		if writeFailed==0:
			print "Will not write"
	
	#print simplex.bestSoln
	print "Took %d iters"%simplex.iterNo
	t4=time.time()
	RMSerror_nm=500.*simplex.fun(simplex.bestSoln)/(2.*numpy.pi)
	print 'RMS/nm =',RMSerror_nm,'nm'
	simplex.showme()
	#raw_input()

	fname=outdir+'/act'+str(1000000+iii)[1:]+'.fits'
	if success!=0 or writeFailed==1:
		print fname
		util.FITS.Write(simplex.bestSoln,fname,extraHeader=["CONVERGD= %d"%success,"RMSERRNM= %g"%RMSerror_nm])
	#util.FITS.Write(simplex.tolarr,fname,writeMode="a")
		iii+=diii
	t5=time.time()
	print "Timings:",t1-t0,t2-t1,t3-t2,t4-t3,t5-t4,"Total:",t5-t0
	#iii+=1

	#f=open('fitting_error'+sys.argv[1]+'.txt','a')
	#f.write(str(RMSerror_nm)+'\n')
	#f.close()



