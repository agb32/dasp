# IMPORTs on any system
import math
import numpy
import numpy as na
import numpy.random as ra
import pylab
import sys
import os
import time
from scipy.special import gamma
sys.path.append('/home/dmw')
import numpyFITS as FITS
from cmod.interp import gslCubSplineInterp_UB
from gp2 import g_surface
import spline


### A couple of rather useful functions ###


def findrms(arr1,arr2):
    mask = circ((numpy.shape(arr2)[0]-1)/2,numpy.shape(arr2)[0])
    phase = (arr1-arr2)*mask
    rms=numpy.sqrt(((numpy.sum(phase**2))/numpy.sum(mask))-(
                                numpy.sum(phase)/numpy.sum(mask)**2))
    return rms

def circ(r, gsize, outside=0., inside=1.):
    z = numpy.zeros((gsize,gsize))+outside      # Initialise array for aperture.
    for i in range(0,gsize):                    # Loop over all elements
        for j in range(0,gsize):
            if ((i-gsize/2)**2+(j-gsize/2)**2) < r**2:
                z[i,j] = inside     # Elements within 'r' centre set to inside
    return z                        #                               (default 1.)


# makewavefront is replaced by 'makeInitialScreen', below
def makewavefront(gsize, r0, pixelscale):
    L_0 = 30.
    rndgrid = numpy.random.randn(gsize,gsize)
    f = numpy.fft.fftshift(numpy.fft.fftfreq(gsize,d=pixelscale))
    f2D = (f**2+(f.reshape(gsize,1))**2)**0.5
    powerspectrum = 6.883877*((f2D/r0)**(-5./3.))
#    powerspectrum = 0.022883*(r0**(-5./3.))*(L_0**(11./3.))*(1.+
#                                           (L_0*numpy.abs(f2D))**2)**(-11./6.)
    powerspectrum[powerspectrum.shape[0]/2,powerspectrum.shape[0]/2] = 0.
    wavefront = numpy.abs(numpy.fft.ifft2(powerspectrum*rndgrid))
    return wavefront

# Python function for linear interpolation... would be quicker wrapped into C
def linint(xvals,yvals,xinterp):
    yinterp = numpy.zeros((numpy.shape(xinterp)))
    y1vals = numpy.zeros(numpy.shape(xvals)[0])
    for i in range(numpy.shape(y1vals)[0]-1):
        y1vals[i]=(yvals[i+1]-yvals[i])/(xvals[i+1]-xvals[i])
    spacing = (numpy.shape(xinterp)[0]-1)/(numpy.shape(xvals)[0]-1)
    for n in range(numpy.shape(xinterp)[0]):
        x = xinterp[n]
        lo = math.floor(x)
        xindex = int( n/spacing )
        yinterp[n]= yvals[xindex]+ (x-lo)*y1vals[xindex]
    return yinterp


#   borrowed code for phasescreen generation...
def makeInitialScreen(dpix=1024,Dtel=42.,L0=30.,scrnXPxls=None,scrnYPxls=None,
                      seed=0,tstep=0.05,globR0=0.2,strLayer=1.,
                      natype=na.float64,windDirection=0.,vWind=10.):
    """dpix is the telescope aperture diameter, and dtel is tel diameter.
    The actual number of pixels used is scrnXPxls x scrnYPxls.
    """
    if scrnXPxls==None:
        scrnXPxls=dpix
    if scrnYPxls==None:
        scrnYPxls=dpix
    scrnPxls=max(scrnXPxls,scrnYPxls)
    pixScale=Dtel/float(dpix)
    ro=globR0*(strLayer**(-3./5.))##we compute the ro in the considered layer
    colAdd=-vWind*na.cos(windDirection*na.pi/180)/pixScale*tstep
    #number of pixels to step each iteration (as float).
    rowAdd=-vWind*na.sin(windDirection*na.pi/180)/pixScale*tstep
    #number of pixels to step each iteration (as float).


    if seed==None:
        print "ERROR (possibly): computeInitialScreen - seed is None, so \
        timer will be used, meaning that the initial screen cannot be \
        replicated, so if both infScrn and infAtmos try to create, you \
        will get a bad phasescreen.  If you wish to use a random seed, \
        use int(time.time()) in the parameter file - though this will \
        only work if all running in the same process."
    if L0>=Dtel:
        scrnSize=2*L0
    elif Dtel>=2*L0:
        scrnSize=Dtel
    else:
        scrnSize=2*L0
    ##size in pixels of the required phase screen
    nfft=int(na.around(dpix*1.*scrnSize/Dtel))#agb: scrnPxls was dpix
    nfft=max(nfft,scrnPxls)#now choose the maximum size...
    print "Initial phase screen size %d"%nfft
    ##gaussian standard 2D random variable
    ra.seed(seed)
    rn=ra.randn(nfft,nfft)

    ##creation of Von Karman spectrum with power spectrum in
    ##                                  (f**2+fo**2)**(-11/6.)
    ##circular grid
    axe_x=na.arange(nfft)-nfft/2.
    f=na.sqrt((axe_x[:,na.newaxis])**2.+(axe_x[na.newaxis,:])**2.)

    #for FFT computations
    ff = na.fft.fftshift(f) # to have 0 frequency in lower left corner
    fo=1./L0; #f0 = 1/L0

    ##factor in the power spectrum (0.023)
    fact=(gamma(11./6.))**2.;
    fact/=2*na.pi**(11./3.);
    fact*=((24./5)*gamma(6./5))**(5./6);

    ##f0 in pixels for power spectrum (cf FA's thesis, page 327 for explanation)
    phozero=fo*1.*scrnSize;

    ##we compute the phase power spectrum
    f = (ff**2.+phozero**2.)**(-11./12.)
    f*=na.sqrt(fact*((scrnSize/ro)**(5./3.))*(nfft*1.)**2.)
    ## Von Karman phase power spectrum : 0.023/ro*(f^2+fo^2)^(-11/6)

    ##FFT of white noise
    ff = na.fft.fft2(rn);

    ## coloration of white noise by Von Karman spectrum
    ff.real*=f;
    ff.imag*=f

    ##inverse 2D Fourier Transform
    phi = (na.fft.ifft2(ff)).real
    #phi=phi.astype(natype)

    #if colAdd is less than zero, we are adding now columns on the right of
    #   the array.
    #If rowAdd is less than zero, we are adding new rows at the bottom of
    #   the array
    #(note, that if plotting in Gist, this is the top of the array).

    phaseArray=na.zeros((scrnYPxls+1,scrnXPxls+1),natype)#was dpix agbc swapped
    if colAdd<0:
        if rowAdd<0:
            phaseArray[1:,1:]=phi[:scrnYPxls,:scrnXPxls]#was dpix
        else:
            phaseArray[:-1,1:]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    else:
        if rowAdd<0:
            phaseArray[1:,:-1]=phi[:scrnYPxls,:scrnXPxls]#was dpix
        else:
            phaseArray[:-1,:-1]=phi[:scrnYPxls,:scrnXPxls]#was dpix
    return phaseArray






### Some classes ###
#
# 1. mirror             - class to hold a dm set of actuators and faceplate, 
#                       and methods of making the actuators bend the faceplate
#                       to fit a given atmosphere.
# 2. linearmirror       - inherits fitting methods, and gives adds the function
#                       'interpolate' which models faceplate using linear
#                        interpolation.
# 3. splinemirror       - inherits fitting methods, and gives adds the function
#                       'interpolate' which models faceplate using cubic
#                        spline interpolation.
# 4. quicksplinemirror  - inherits fitting methods, and gives adds the function
#                       'interpolate' which models faceplate using wrapped C
#                        code for cubic spline interpolation.
# 5. gaussianmirror     - inherits fitting methods, and gives adds the function
#                       'interpolate' which models faceplate using (gaussian)
#                        influence function addition.


class mirror:
    def __init__(self,nact=8,spacing=8,shres=None,target=None):
#       nact = actuators across faceplate (inclusive of those at edge)
#       spacing = # of interpolation points between actuators (inclusive of
#       points at actuators)
#       shres = shank-hartman resolution... pixels required across phasescreen
#       target = phasescreen or other shape mirror is to be fitted to
        self.actuators=None
        self.nact       = nact
        self.spacing    = spacing
        if shres==None:
            shres = 1+spacing*(nact-1)
        self.shres=shres
        if target==None:
            self.target = makewavefront(int(
                2*(2**math.ceil(math.log(shres)/math.log(2)))),
                                        0.05)[0:shres,0:shres]
        else:
            self.target = target
        self.calls = 0 # Parameter to find how many calls to the interpolation
#                                                           routine are made
                  
    def getactuators(self,largegrid,n_actuators,int_act):
#       Function to produce an initial guess of optimal actuator heights by 
#       taking the phasescreen values at actuator positions. Messy... 
#       copies data for relevant rows and then relevant columns from that.
        reducex = numpy.zeros((n_actuators, numpy.shape(largegrid)[0]))
        reduced = numpy.zeros((n_actuators, n_actuators))
        for n in range(0,numpy.shape(largegrid)[0],int_act):
            reducex[n/int_act,:]=largegrid[n,:]
        for m in range(0,numpy.shape(largegrid)[0],int_act):
            reduced[:,m/int_act]=reducex[:,m]
        self.actuators = reduced
        return reduced

    def fit_metropolis(self,schedule):
#       Schedule is a tuple of: iterations, zerotemp iterations,
#       iterations per temp, start temp, start elasticity & base elasticity.

#       Minimisation function based on simulated annealing.
        arr_1   = self.target
        arr_2   = self.actuators
        sh_res  = self.shres
        spacing = self.spacing
        #rmslog = []
        anneallog = []
        interpolation = self.interpolate(arr_2,sh_res,spacing)
        old_rms = findrms(arr_1,interpolation)
        start_temp = (schedule[0]-schedule[1])/schedule[2]
        for t in range(start_temp,0,-1):
            for k in range(schedule[2]):
                # Pick a random actuator
                x = numpy.random.random_integers(0,numpy.shape(arr_2)[0]-1)
                y = numpy.random.random_integers(0,numpy.shape(arr_2)[1]-1)
                act_trial = numpy.copy(arr_2)
                # Poke it randomly
                act_trial[x,y]=act_trial[x,y]+((schedule[5]+(t*schedule[4])/
                                            start_temp))*numpy.random.randn()
                # Does this make a better fit?
                interpolation_trial= self.interpolate(act_trial,sh_res,spacing)
                new_rms = findrms(arr_1,interpolation_trial)
                d_rms = new_rms - old_rms
                # Do it if it does, if not then do it with probability
                #                           depending on how bad it is...
                if math.exp((-d_rms*start_temp)/(t*schedule[3])) \
                                            > numpy.random.uniform():
                    arr_2 = numpy.copy(act_trial)
                    interpolation = interpolation_trial
                    old_rms = new_rms
                rmslog.append(old_rms)
            anneallog.append([t, old_rms])

        # Zero temp iterations...    
        for k in range(schedule[2]):
            # Pick a random actuator
            x = numpy.random.random_integers(0,numpy.shape(arr_2)[0]-1)
            y = numpy.random.random_integers(0,numpy.shape(arr_2)[1]-1)
            act_trial = numpy.copy(arr_2)
            # Poke it randomly
            act_trial[x,y]=act_trial[x,y]+schedule[5]*numpy.random.randn()
            # Does this make a better fit?
            interpolation_trial= self.interpolate(act_trial,sh_res,spacing)
            new_rms = findrms(arr_1,interpolation_trial)
            d_rms = new_rms - old_rms
            # Do it ONLY if it does.
            if d_rms < 0.:
                arr_2 = numpy.copy(act_trial)
                interpolation = interpolation_trial
                old_rms = new_rms
            rmslog.append(old_rms)
            if k%10 == 0:
                anneallog.append([t, old_rms])
        self.actuators = arr_2
        return arr_2, rmslog

    def fit_steepest(self,maxitns):
#       Another working, but poor, minimisation function.
#       Works by steepest descent...
        self.calls = 0
        arr_1   = self.target
        arr_2   = self.actuators
        sh_res  = self.shres
        spacing = self.spacing
        #rmslog = []
        interpolation = self.interpolate(arr_2,sh_res,spacing)
        start_rms = findrms(arr_1,interpolation)
        rmslog.append(start_rms)
        grad = numpy.zeros(numpy.shape(arr_2))
        # An array to hold the n_act*n_act gradient of the function findrms...

#       Firstly find a gradient...
	for m in range(maxitns):
            for i in range(numpy.shape(arr_2)[0]):
                    for j in range(numpy.shape(arr_2)[0]):
                        arr_poke = numpy.copy(arr_2)
#                       Poke an actuator a little bit...
                        arr_poke[i,j] = arr_poke[i,j] * 1.0002
#                       And whatever difference it makes... record it.                        
                        grad[i,j] = (start_rms - findrms(
                            arr_1,self.interpolate(arr_poke,sh_res,spacing)))/ \
                            (arr_2[i,j] / 1000.)
            stepsize = 8000.
            for n in range(1000):
#               Then go in that direction
                new_rms = findrms(self.interpolate(arr_2+grad*stepsize,
                                                   sh_res,spacing),arr_1)
                if new_rms < start_rms:
#                   And keep going in that direction until it doesn't go
#                                                           down any more
                    arr_2 = arr_2 + grad*stepsize
                    rmslog.append(new_rms)
                    stepsize = 4000./float(m+1) # Alter stepsize to be
#                                               smaller closer to a minimum...
                    start_rms = new_rms
                else:
                    print ('Exited at grad addition number ' + str(n))
                    break
            if ((rmslog[len(rmslog)-n]) - new_rms) < 0.005:
                # Termination Criterion
                break
        self.actuators = arr_2
        return arr_2, rmslog

    def fit_polak(self, maxitns, lmintol=0.001, tol=0.01, poke=1.0002):
#       A rather more refined gradient method. This uses conjugate
#       gradients and golden section line minimisation.
        self.calls = 0
        arr_1   = self.target
        arr_2   = self.actuators
        sh_res  = self.shres
        spacing = self.spacing
        #rmslog = []
        interpolation = self.interpolate(arr_2,sh_res,spacing)
        start_rms = findrms(arr_1,interpolation)
        grad = numpy.zeros(numpy.shape(arr_2))
# Find steepest descent
        for i in range(numpy.shape(arr_2)[0]):
                for j in range(numpy.shape(arr_2)[0]):
#                   Poke an actuator a little bit...
                    arr_poke = numpy.copy(arr_2)
                    arr_poke[i,j] = arr_poke[i,j] * poke
#                   And whatever difference it makes... record it.   
                    grad[i,j] = (start_rms - findrms(
                        arr_1,self.interpolate(arr_poke,sh_res,spacing))) / \
                        (arr_2[i,j] / 1000.)
        print 'gradient found'

#       Declare vectors to hold gradent, G
        G = numpy.copy(grad)
#       and to hold conjugate gradient, H
        H = numpy.copy(G)

#       Polak-Ribiere Algorithm:
	for m in range(maxitns):        # Only do it so many times...
#       Golden Section Line Minimisation
# First off bracket the minimum...
            R = 0.61803399 # phi-1
            C = 1.-R
            x0 = 0.     #First  point on line
            x1 = 10000. #Second point on line
            # Find f(x0) and f(x1)
            f0 = findrms(self.interpolate(arr_2,sh_res,spacing),arr_1)
            f1 = findrms(self.interpolate(arr_2+grad*x1,sh_res,spacing),arr_1)
            # Ensure f is greater at x0, swap x0 and x1 if necessary.
            if f0 < f1:
                ft = f1
                f1 = f0
                f0 = ft
                xt = x1
                x1 = x0
                x0 = xt
            gold_power = (1.+R)
            x3 = x0+(gold_power*(x1-x0))
            f3 = findrms(self.interpolate(arr_2+grad*x3,sh_res,spacing),arr_1)
            # Keep steping away from x0 until a bracketing triplet is found
            while f3 < f1:
                gold_power = gold_power * (1.+R)
                x3 = x0+(gold_power*(x1-x0))
                f3 = findrms(self.interpolate
                             (arr_2+grad*x3,sh_res,spacing),arr_1)
            print 'minimum bracketed'

# Then do a golden section search for it...
            # Pick which interval to start looking in
            if abs(x3-x1) > abs(x1-x0):
                x2 = x1+C*(x3-x1)
                f2 = findrms(self.interpolate
                             (arr_2+grad*x2,sh_res,spacing),arr_1)
            else:
                x2 = x1
                x1 = x1-C*(x3-x1)
                f2 = f1
                f1 = findrms(self.interpolate
                             (arr_2+grad*x1,sh_res,spacing),arr_1)
            # While bracketing triplet is larger than specified accuracy...
            while abs(x3-x0) > lmintol*(abs(x1)+abs(x2)):
            # Depending on which side of x1 the higher point is...
                if f2 < f1:
                    # Shuffle points one way
                    x0 = x1
                    x1 = x2
                    x2 = R*x1 + C*x3
                    f0 = f1
                    f1 = f2
                    f2 = findrms(self.interpolate
                                 (arr_2+grad*x2,sh_res,spacing),arr_1)
                else:
                    # or shuffle them the other...
                    x3 = x2
                    x2 = x1
                    x1 = R*x2 + C*x0
                    f3 = f2
                    f2 = f1
                    f1 = findrms(self.interpolate
                                 (arr_2+grad*x1,sh_res,spacing),arr_1)
            # Return whichever of the resulting points is lower
            if f1 < f2:
                x = x1
            else:
                x = x2
            # Take actuators that far in the direction of gradient
            arr_2 = arr_2 + grad*x
            new_rms = findrms(self.interpolate(arr_2,sh_res,spacing),arr_1)
            print 'iteration ', m, ' : rms = ', new_rms
            if start_rms - new_rms < tol:#numpy.average(grad**2) < 1.e-9:
                break
            start_rms = new_rms

            
#new direction... the Polak-Ribiere magic
            # find gradient again (should code this as a function...)
            for i in range(numpy.shape(arr_2)[0]):
                for j in range(numpy.shape(arr_2)[0]):
                    arr_poke = numpy.copy(arr_2)
                    arr_poke[i,j] = arr_poke[i,j] * 1.0002
                    grad[i,j] = (start_rms -
                                 findrms(arr_1,self.interpolate
                                         (arr_poke,sh_res,spacing))) / \
                                         (arr_2[i,j] / 1000.)
            GG = 0.
            DGG = 0.
            # then do this...
            for i in range(numpy.shape(arr_2)[0]):
                for j in range(numpy.shape(arr_2)[0]):
                    GG = GG + G[i,j]**2
                    DGG = DGG + (grad[i,j]+G[i,j])*grad[i,j]
            GAM = DGG/GG
            G = grad
            H = G+GAM*H
            # and the new conjugate gradient to follow is:
            grad = H
            
        self.actuators = arr_2
        return arr_2, rmslog



class linearmirror(mirror):
    def interpolate(self,actuators,width,int_act):
        self.calls = self.calls + 1
        n_act = numpy.shape(actuators)[0]
        
        xinterp = numpy.zeros((width,n_act))
        yinterp = numpy.zeros((width,width))

        x = numpy.arange(n_act)
        xint = numpy.arange(0,width)/float(int_act)
        for row in range(0,n_act):
            y = actuators[:,row]
            xinterp[:,row]=linint(x,y,xint)

        x = numpy.arange(n_act)
        xint = numpy.arange(0,width)/float(int_act)
        for col in range(0,width):
            y = xinterp[col,:]
            yinterp[col,:]=linint(x,y,xint)

        return yinterp
    

class splinemirror(mirror):
    def interpolate(self,actuators,width,int_act):
        self.calls = self.calls + 1
        n_act = numpy.shape(actuators)[0]
        
        xinterp = numpy.zeros((width,n_act))
        yinterp = numpy.zeros((width,width))

        x = numpy.arange(n_act)
        xint = numpy.arange(0,width)/float(int_act)
        for row in range(0,n_act):
            y = actuators[:,row]
            xinterp[:,row]=spline.splint(x,y,xint) 
# splint is Python coding of a Numerical Recipies function
        x = numpy.arange(n_act)
        xint = numpy.arange(0,width)/float(int_act)
        for col in range(0,width):
            y = xinterp[col,:]
            yinterp[col,:]=spline.splint(x,y,xint)

        return yinterp



class quicksplinemirror(mirror):
    def interpolate(self,actuators,width,int_act):
        self.calls = self.calls + 1
	n_act = numpy.shape(actuators)[0]
        x = y = numpy.arange(width).astype("d")/(width-1)
	x2 = numpy.arange(n_act).astype("d")/(n_act-1)
	phsOut = numpy.zeros((width,width),"d")
	gslCubSplineInterp_UB(actuators,x2,x2,y,x,phsOut,4)
	return phsOut



class gaussmirror(mirror):
    # need to redefine __ini__ to have gcoeff (Gaussian width)
    def __init__(self,nact=8,spacing=8,shres=None,target=None,gcoeff=24.):
        self.actuators=None
        self.nact       = nact
        self.spacing    = spacing
        self.gcoeff     = gcoeff
        if shres==None:
            shres = 1+spacing*(nact-1)
        self.shres=shres
        if target==None:
            self.target = makewavefront(int(2*(2**math.ceil(
                math.log(shres)/math.log(2)))),0.05)[0:shres,0:shres]
        else:
            self.target = target
        self.calls = 0

    def interpolate(self,actuators,width,int_act):
        self.calls = self.calls + 1
        faceplate = numpy.zeros((width,width))
        # Need to use actuator deviations from average... so find it
        av = numpy.average(actuators)
        zactuators= actuators - av
	table = numpy.zeros((width*2-1,width*2-1))
	timevar1= time.clock()
	# This works quicker if you create a lookup table for the influence
	# function...
	if self.calls < 2:
            for m in range(width*2-1):
                for n in range(width*2-1):
                    table[m,n]=math.exp(
                        -((n-width+1)**2 + (m-width+1)**2)/self.gcoeff)
            print 'Creating the table took: ', time.clock() - timevar1, ' s'
	timevar1= time.clock()
	# Then add the influence function of each actuator (2 for loops)
	# to the faceplate at each pixel (the other 2 for loops)
        for i in range(self.shres):
            for j in range(self.shres):
                for k in range(self.nact):
                    for l in range(self.nact):
                        faceplate[i,j] = faceplate[i,j] + \
                            zactuators[k,l]*table[(i-8*k+width-1),
                                                  (j-8*l+width-1)]
        print 'Adding up took: ', time.clock() - timevar1, ' seconds.'
        faceplate = faceplate + av
        return faceplate



    
### The actual program ###

# Useful annealing Scedules
#Schedule is a tuple of: iterations, zerotemp iterations, iterations per temp,
#start temp, start elasticity & base elasticity
schedule1 = (8000,400,10,0.016,3.16,1)
schedule2 = (8000,400,20,0.016,3.16,1)
scheduleC = (4000,200,20,0.004,6,0.2) # Only takes a reasonable length of time
                                      # with C-coded spline interpolation 
scheduleq = (4000,400,10,0.016,4,1)

# Pick one...
schedule = scheduleC

rmslog=[]

phasescreen = makeInitialScreen(dpix=57,Dtel=4.2,L0=30.,
                                seed=int(time.time()),globR0=0.085)[0:57,1:58]
dm = quicksplinemirror(target=phasescreen)
#dm = gaussmirror(target=phasescreen)
#dm = linearmirror(target=phasescreen)
print 'Mirror Initialised'
actuators = dm.getactuators(dm.target,dm.nact,dm.spacing)
timevar = time.clock()
reconstruction = dm.interpolate(dm.actuators,dm.shres,dm.spacing)
timevar = time.clock() - timevar
print 'Interpolation took ', timevar, 'seconds.'
initial_rms = findrms(dm.target,reconstruction)
print 'Initial rms =   ', initial_rms
print 'Fitting Mirror'
timevar = time.clock()
result = dm.fit_polak(20)
timevar = time.clock() - timevar
fittedactuators = result[0]
fittedreconstruction = dm.interpolate(fittedactuators,dm.shres,dm.spacing)
print 'Polak  fit #1', findrms(dm.target,fittedreconstruction)
print '...in ', timevar, ' seconds, and '
print dm.calls, ' calls to the mirror interpolation.'
print 'Mirror Fitted,'
final_rms = findrms(dm.target,fittedreconstruction)
print 'Final rms =   ', final_rms


### Standard Outputs ###
run_number = 'SPLINE01' # File referencing string
FITS.Write(dm.target,
           'phasescreen' + run_number +'.FITS')
FITS.Write(actuators,
           'actuators' + run_number +'.FITS')
FITS.Write(reconstruction,
           'reconstruction' + run_number +'.FITS')
FITS.Write(reconstruction-dm.target,
           'difference' + run_number +'.FITS')
FITS.Write(fittedactuators,
           'fittedactuators' + run_number +'.FITS')
FITS.Write(fittedreconstruction,
           'fittedreconstruction' + run_number +'.FITS')
FITS.Write(fittedreconstruction-dm.target,
           'fitteddifference' + run_number +'.FITS')
pylab.plot(result[1])
pylab.savefig('rmslog' + run_number + '.ps')
os.system('ps2pdf rmslog' + run_number + '.ps')
f = open('/home/dmw/AO/log' + run_number, 'w')
f.write('FIT PARAMETERS\n')
f.write('Starting rms:\n')
f.write(str(initial_rms)+'\n')
f.write('Finishing rms:\n')
f.write(str(final_rms)+'\n')
f.close()
print 'Done!'

#pylab.imshow(phasescreen)
#pylab.imshow(reconstruction)
#pylab.imshow(phasescreen-reconstruction)
