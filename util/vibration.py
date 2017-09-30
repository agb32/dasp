import numpy
import util.zernikeMod
class Vibration:
    """Class to introduce tip and tilt vibration - modified from Bharmal code."""
    def __init__(self,pupil,telDiam,tstep,freqList=[10.],angleList=[0.05],phaseList=[0.],dirnList=[0.]):
        """telDiam - in metres
        wfslam - in nm.
        npup - in pixels
        tstep - in seconds
        freqList - in Hz
        angleList - in arcsec
        phaseList - in radians
        dirnList - in degrees.
        """
        self.itercnt=0
        npup=pupil.shape[0]
        self.pupil=pupil
        self.telDiam=telDiam
        #self.wfslam=wfslam
        self.npup=npup
        self.tstep=tstep

        self.freqList=freqList
        self.angleList=angleList
        self.phaseList=phaseList
        self.dirnList=dirnList

        z=util.zernikeMod.Zernike(pupil,3,computeInv=0)
        z.zern[1]*=(npup-1)/2./z.zern[1].max()#scale to Bharmal values.
        z.zern[2]*=(npup-1)/2./z.zern[2].max()#scale to Bharmal values.
        self.ttModes=z.zern[1:]
        
    def angScaling(self,ang, phs, step, freq,lam):
        scale=(ang*4.848e-6*self.telDiam)*(2*numpy.pi/lam*1e9)\
            *(2*self.npup**-1.0)*numpy.sin( phs + 2*numpy.pi*step*self.tstep*freq )
        return scale


    def vibrationMode(self,ang, phs, freq, rot,lam):
        # \/ args:- angle, phase, frequency, and rotation (direction)
        #    as [arcsec],[rad],[Hz],[rad]

        scale=self.angScaling(ang,phs,self.itercnt,freq,lam)
        t=numpy.cos(rot*numpy.pi/180.)*self.ttModes[0]+numpy.sin(rot*numpy.pi/180.)*self.ttModes[1]
        t*=scale
        return t

    def addVibration(self,out=None,lam=None):
        """Out is the phase to be added too, and lam is the wavelength in nm."""
        if out is None:
            out=numpy.zeros((self.npup,self.npup),numpy.float32)
        for i in range(len(self.freqList)):
            out+=self.vibrationMode(self.angleList[i],self.phaseList[i],self.freqList[i],self.dirnList[i],lam)
        self.itercnt+=1
        return out
            
