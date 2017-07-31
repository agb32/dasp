#Import important packages
import numpy
import scipy
import subprocess
import os
import glob
import re
import random
import shutil
from scipy import interpolate
from astropy.io import fits
import re
import aplpy
from images2gif import writeGif
from PIL import Image
import matplotlib.pyplot as plt
from skimage.transform import resize
from decimal import getcontext, Decimal
import signal
from numpy.matlib import repmat
import operator

#Photosensor Module
class Photosensor:
    def __init__(self,uin,telescopeprimirror,telescopefocalength,telescopereflectcoef,exposuretime,wavelength,pixel_size,sensortype='CCD'):
        #All the photosensor with initialization of:
        #cccd=Photosensor(psf,telescopeprimirror,telescopefocalength,telescopereflectcoef,\
        #exposuretime,wavelength,ccdpixelsize,sensortype='CCD')
        #And run with:
        #ccdimg=cccd.outputimg(psf)
        #Reference to paper: High leverl numerical simulations of noise in CCD and 
        #CMOS photosensors: review and tutorial By Konnik and Welsh
        #Basic information get from input
        ccdpixelnumw,ccdpixelnuml=numpy.shape(uin)
        #Basic Physical parameter definition
        #<-----define metric
        self.m=1
        self.cm=1e-2*self.m
        self.mm=1e-3*self.m
        self.mum=1e-6*self.m
        self.nm=1e-9*self.m
        self.rad=1
        self.mrad=1e-3*self.rad
        #<-----define physical constant that would be useful
        self.h=6.62606896*10**(-34) #Plank's constant, in [Joule*s]
        self.c=2.99792458*10**8     #speed of light, in [m/s]
        self.Boltzman_Constant=8.617343*10**(-5) #Boltzman constant, [eV/K].
        self.Boltzman_Constant_JK=1.3806504*10**(-23) #Boltzman constant, [J/K].
        self.q=1.602176487*10**(-19)   #a charge of an electron [C], Cylon
        self.k1=10.909*10**(-15)
        self.Signal_CCD_photons= numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        #Used to control whether we need particular data
        darkframe=0
        Venonlinearity=0
        VVnonlinearity=0
        ADCnonlinearity=0
        irradiance=1
        electrons=1
        volts=1 
        fDN=1
        plots=[irradiance,electrons,volts,fDN]
        #Basic definition part, change if you really know these parameters
        #Noise that considered 1 for calculation and 0 for no
        photonshotnoise=1
        fPRNU=1
        darkcurrent=1
        darkcurrent_Dshot=1
        darkcurrent_DarkFPN_pixel=1
        darkcurrent_offsetFPN=1
        sensenoderesetnoise=1
        #Signal parameters need to take into consideration
        self.telescopeaper=telescopeprimirror*self.m
        self.telescopefocal=telescopefocalength*self.m
        self.wavelength=wavelength*self.mum
        self.telescopereflectcoef=telescopereflectcoef
        #Noise model and parameter for nPRNU
        modle='Janesick-Gaussian'
        nparameters=[]
        factor=0.01
        n_noisematrix=numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        self.nPRNU=[modle,nparameters,factor,n_noisematrix]
        #Noise model and parameter for nDN
        nDN=0.3
        model='Janesick-Gaussian'
        dparameters=[]
        noisematrix=0
        self.darkFPN=[nDN,model,dparameters,noisematrix]
        #Another possible model is logNormal model and you can change that with the commonents below
        """
        model='LogNormal'
        parameters=[0,0.4] 
        model='Wald'
        parameters=2
        """
        #Dark FPN_offset
        model='Janesick-Gaussian'
        parameters=[]
        DNcolumn=0.0005
        dfnoisematrix=0
        self.darkFPN_offset=[model,parameters,DNcolumn,dfnoisematrix]
        Factor=0.8
        sigma_ktc=0
        sr_noisematrix=0
        self.sn_reset=[Factor,sigma_ktc,sr_noisematrix]
        #Nonlinearity
        self.A_SNratio=0.05
        self.A_SFratio=1.05
        self.ADCratio=1.1
        #Flag for all the noise we considered 
        self.flag=[photonshotnoise,fPRNU,darkcurrent,darkcurrent_Dshot,\
                   darkcurrent_DarkFPN_pixel,darkcurrent_offsetFPN,sensenoderesetnoise,\
                   plots,electrons,volts,fDN,darkframe,Venonlinearity,VVnonlinearity,ADCnonlinearity]
        #Used for noise model set
        self.noise=[self.nPRNU,self.darkFPN,self.darkFPN_offset,self.sn_reset]
        #Photo detector parameters
        self.SensorType=sensortype
        #Do it later
        pixel_size=pixel_size*self.telescopefocal*numpy.pi/180/3600*self.mum
        self.pixel_size=[pixel_size, pixel_size] #pixels size
        self.t_I=exposuretime #Exposure/Integration time
        self.QE_I=0.8  #Quantum Efficiency of the photo sensor
        self.FillFactor=0.5 #Pixel Fill Factor for CMOS photo sensors.
        self.QuantumYield=1  #quantum yield (number of electrons per one photon interaction)
        self.FW_e=2*10**4    #full well of the pixel
        self.V_REF=3.1        #Reference voltage to reset the sense node(3-10V)
        #Sense Nose 
        self.A_SN=5*10**(-6) #Sense node gain, A_SN [V/e]
        self.A_SF=1         #Source follower gain, [V/V]
        self.A_CDS=1         #Correlated Double Sampling gain, [V/V]
        self.N_bits=12        #noise is more apparent on high Bits
        self.S_ADC_OFFSET=0   #Offset of the ADC, in DN
        #Temperature and dark current noise parameters
        self.T=300  #operating temperature, [K]
        self.DFM=1  #dark current figure of merit
        self.Eg_0=1.1557 #Sensor material constants 
        self.alpha= 7.021*10**(-4) # material parameter, [eV/K]
        self.beta= 1108            # material parameter, [K]
        #Define some maps/images that may be useful
        self.Signal_CCD_photons=numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        self.Signal_CCD_electrons=numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        self.dark_signal=numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        self.nonlinearity=[self.A_SNratio,self.A_SFratio,self.ADCratio]
        sensor_size=[ccdpixelnumw,ccdpixelnuml]
        self.light_signal=numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        self.Signal_CCD_voltage=numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        self.Signal_CCD_DN=numpy.zeros([ccdpixelnumw,ccdpixelnuml])
        
    #Define a circ tool for calculation
    def tool_circ(self,x,y,D):
        r=numpy.sqrt(x**2+y**2)
        z=(r<D/2).astype(numpy.float)
        z[r==D/2]=0.5
        return z
        
    #PART I Photon to electrons process
    #Define ccd illumination energy here
    def ccd_illumination_prepare(self,N,M):
        #What calculate here is the photon output percentage of the telescope
        delta_x=self.primirror*1.0/N
        delta_y=self.primirror*1.0/M
        dx=numpy.arange(-N/2*delta_x,N/2*delta_x,delta_x)
        dy=numpy.arange(-N/2*delta_y,N/2*delta_y,delta_y)
        [x1,y1]=numpy.meshgrid(dx,dy)
        uin = self.telescope.reflectcoef*self.tool_circ(x1,y1,self.primirror)
        if Cccd.flag[8][0]==1:
            uin_irradiance=numpy.abs(uin)**2
        return uin
                
        #Define different random number generator here
    def tool_rand_distributions_generator(self,distribName, distribParams, sampleSize):
        funcName ='tool_rand_distributions_generator'
        distribNameInner = distribName.lower()
        out=[]
        if numpy.prod(sampleSize)>0:
            if distribNameInner in {'exp','exponential'}:
                self.checkParamsNum(funcName,'Exponential','exp',distribParams,[1]) 
                lambda2 =distribParams[0]
                validateParam(funcName,'Exponential','exp','lambda','lambda',lambda2,{'> 0'})
                out=-numpy.log(numpy.random.rand(sampleSize))/lambda2
            elif distribNameInner in {'lognorm','lognormal','cobbdouglas','antilognormal'}:
                self.checkParamsNum(funcName, 'Lognormal', 'lognorm', distribParams, [0, 2]);
                if len(distribParams)==2:
                    mu=distribParams[0]
                    sigma=distribParams[1]
                    self.validateParam(funcName,'Lognormal','lognorm','[mu,sigma]', 'sigma',sigma,{'> 0'})
                else:
                    mu = 0
                    sigma = 1
                out=numpy.exp(mu+sigma*numpy.ramdom.randn( sampleSize))
            elif distribNameInner in {'ig', 'inversegauss', 'invgauss'}:
                self.checkParamsNum(funcName,'Inverse Gaussian','ig',distribParams,[2])
                theta=distribParams[0]
                chi=distribParams[1]
                self.validateParam(funcName,'Inverse Gaussian','ig','[theta,chi]','theta', theta, {'> 0'})
                self.validateParam(funcName, 'Inverse Gaussian', 'ig', '[theta, chi]', 'chi', chi, {'> 0'})
                chisq1=numpy.random.randn(sampleSize)**2;
                out=theta+0.5*theta/chi*( theta*chisq1-numpy.sqrt(4*theta*chi*chisq1 + theta**2*chisq1**2) )
                l=numpy.random.rand(sampleSize)>=theta/(theta+out)
                out[l]=theta**2/out[l]
            elif distribNameInner in {'logistic'}:
                self.checkParamsNum(funcName, 'Logistic', 'logistic', distribParams, [0, 2]);
                if len(distribParams)==2:  #numel
                    a=distribParams[0]
                    k=distribParams[1]
                    self.validateParam(funcName, 'Laplace', 'laplace', '[a, k]', 'k', k, {'> 0'})
                else:
                    a=0
                    k=1
                u1=numpy.random.rand(sampleSize)
                out=a-k*numpy.log( 1/u1 -1)
            elif distribNameInner in {'wald'}:
                self.checkParamsNum(funcName, 'Wald', 'wald', distribParams, [1])
                chi = distribParams[0]
                self.validateParam(funcName, 'Wald', 'wald', 'chi', 'chi', chi, {'> 0'})
                out=tool_rand_distributions_generator( 'ig', [1,chi], sampleSize)
            else:
                print('\nRANDRAW: Unknown distribution name: ', distribName)
        return out
        
    #Define parameter checker for random parameter check
    def  checkParamsNum(self,funcName, distribName, runDistribName, distribParams, correctNum):
        if ~any(len(distribParams) == correctNum):
            print(distribName,'Variates Generation:\n','Wrong numebr of \
            parameters (run',funcName,'(',runDistribName,') for help)')
            return 0
            
    #Define parameter validate method here
    def validateParam(self,funcName,distribName,runDistribName,distribParamsName,paramName,param,conditionStr):
        condLogical = 1
        eqCondStr = []
        for nn in range(0,len(conditionStr)):
            if nn==1:
                eqCondStr = [eqCondStr or conditionStr]
            else:
                eqCondStr = [eqCondStr and conditionStr]         
            eqCond = conditionStr[0]
            if eqCond=={'<'}:
                condLogical = condLogical & (param< float(conditionStr[2:]))
            elif eqCond=={'<='}:
                condLogical = condLogical & (param<=float(conditionStr[2:]))              
            elif eqCond=={'>'}:
                condLogical = condLogical & (param> float(conditionStr[2:])) 
            elif eqCond=={'>='}:
                condLogical = condLogical & (param>=float(conditionStr[2:]))
            elif eqCond=={'~='}:
                condLogical = condLogical & (param!=float(conditionStr[2:]))
            elif eqCond=={'=='}:
                if cmp(conditionStr[2:],'integer')==0:
                    condLogical = condLogical & (param==numpy.floor(param))          
                else:
                    condLogical = condLogical & (param==float(conditionStr[2:]))
        if condLogical==0:
            print(distribName,'Variates Generation:tool_rand_distributions_generator,(',runDistribName,distribParamsName,\
                    'SampleSize)\n Parameter paramName should be eqCondStr\n (run',funcName,'(',runDistribName,')for help)')
        return 0
        
    #Set constants here for easier calculation
    def ccd_set_photosensor_constants(self,uin):
        #Maybe duplicate
        fN=numpy.shape(uin)[0]
        fM=numpy.shape(uin)[1]
        self.Signal_CCD_photons=numpy.zeros([fN,fM])
        self.Signal_CCD_electrons=numpy.zeros([fN,fM])
        self.dark_signal=numpy.zeros([fN,fM])
        self.nonlinearity=[self.A_SNratio,self.A_SFratio,self.ADCratio]
        self.sensor_size=[fN,fM]
        self.light_signal=numpy.zeros([fN,fM])
        self.Signal_CCD_voltage=numpy.zeros([fN,fM])
        self.Signal_CCD_DN=numpy.zeros([fN,fM])
        return 0
                
    #Define pattern of the FPN (Fixed Pattern Noise)
    def ccd_FPN_models(self, sensor_signal_rows, sensor_signal_columns, \
                        noisetype, noisedistribution, noise_params):
        noiseout=[]
        if noisedistribution=='AR-ElGamal':
            if operator.eq(noisetype,'pixel')==0:
                x2=numpy.random.randn(sensor_signal_rows,sensor_signal_columns)
                noiseout=signal.lfilter([1,0],noise_params,x2)
            elif operator.eq(noisetype, 'column')==0:
                x=signal.lfilter([1,0],noise_params,numpy.random.randn(1,sensor_signal_columns))
                noiseout=repmat(x,sensor_signal_rows,1)
        elif noisedistribution=='Janesick-Gaussian':
            if operator.eq(noisetype, 'pixel')==0:
                noiseout =numpy.random.randn(sensor_signal_rows,sensor_signal_columns)
            elif operator.eq(noisetype, 'column')==0:
                x=numpy.random.randn(1,sensor_signal_columns)
                noiseout = repmat(x,sensor_signal_rows,1)
            elif operator.eq(noisetype, 'row')==0:
                x=numpy.random.randn(sensor_signal_rows,1)
                noiseout = repmat(x,1,sensor_signal_columns)
        elif noisedistribution=='Wald':
            if operator.eq(noisetype,'pixel')==0:
                noiseout = tool_rand_distributions_generator\
                ('wald',noise_params[0],[sensor_signal_rows, sensor_signal_columns])\
                +numpy.random.rand(sensor_signal_rows, sensor_signal_columns)
            elif operator.eq(noisetype, 'column')==0:
                x = tool_rand_distributions_generator\
                ('lognorm',[noise_params[0],noise_params[1]],[1, sensor_signal_columns])
                noiseout = repmat(x,sensor_signal_rows,1)
        return noiseout
        
        #Define light FPN 
    def ccd_photosensor_lightFPN(self):
        self.noise[0][3]=self.ccd_FPN_models(self.sensor_size[0],self.sensor_size[1],\
                                        'pixel',self.noise[0][0],self.noise[0][1])
        self.light_signal=self.light_signal*(1+self.noise[0][3]*self.noise[0][2])
        return  self.light_signal
        
    #Define light noise for ccd 
    def ccd_photosensor_lightnoises(self,uin):
        if operator.eq('CMOS',self.SensorType)==0:
            PA=self.FillFactor*self.pixel_size[0]*self.pixel_size[1]
        else:
            PA=self.pixel_size[0]*self.pixel_size[1]
        uin_irradiance= PA*abs(uin)**2
        P_photon=(self.h * self.c)/self.wavelength
        self.Signal_CCD_photons=numpy.round(uin_irradiance*self.t_I/P_photon)
        if self.flag[0]==1:
            self.Signal_CCD_photons=self.Signal_CCD_photons-numpy.min(self.Signal_CCD_photons)
            self.Signal_CCD_photons=numpy.random.poisson(self.Signal_CCD_photons)
        QE=self.QE_I*self.QuantumYield
        self.light_signal=self.Signal_CCD_photons*QE
        if self.flag[1]==1:
            self.light_signal=self.ccd_photosensor_lightFPN()
        return self.light_signal
        
    #Define the dark FPN noise for ccd
    def ccd_photosensor_darkFPN(self):
        self.noise[1][3] = self.ccd_FPN_models(self.sensor_size[0],self.sensor_size[1], \
                                                'pixel', self.noise[1][1], self.noise[1][3])
        self.dark_signal = self.dark_signal*(1 + (self.noise[1][0])*(self.noise[1][3]))
        return self.dark_signal
        
    #Define the photosnsor dar noise here
    def ccd_photosensor_darknoises(self):
        PA=self.pixel_size[0]*self.pixel_size[1]*10**4
        self.Eg=self.Eg_0-(self.alpha*(self.T**2))/(self.beta+self.T)
        self.DARK_e = (self.t_I)*2.55*10**15*PA*self.DFM*(self.T**(1.5))\
        *numpy.exp(-self.Eg/(2*self.Boltzman_Constant*self.T))
        self.dark_signal = (self.DARK_e)*numpy.ones([len(self.Signal_CCD_electrons),\
                                                  len(self.Signal_CCD_electrons[0])])
        #<----- ### Start:: adding Dark Shot noise
        if self.flag[3]==1:   #Cccd.glag[3]=ccd.flag.darkcurrent_Dshot
            self.dark_signal=numpy.random.poisson(self.dark_signal) 
        #<----- ### END:: adding Dark Shot noise
        #<----- ### Start:: adding Dark FPN  %%% being added to dark current, it is too small.
        if self.flag[4]==1:   #Cccd.flag[4]=Cccd.flag.darkcurrent
            self.dark_signal=self.ccd_photosensor_darkFPN()
        #<----- ### END:: adding Dark FPN  %%% being added to dark current, it is too small.
        return self.dark_signal
        
    #Part II Sensor noise part charge to voltage
    #node sense part
    def ccd_sense_node_chargetovoltage(self):
        self.C_SN=self.q/self.A_SN
        self.V_FW=self.FW_e*self.q/self.C_SN
        self.V_min=self.q*self.A_SN/self.C_SN
        if operator.eq('CMOS',self.SensorType)==0:
            if self.flag[6]==1:  #Cccd.flag[6]=ccd.flag.sensenoderesetnoise
                if self.noise[3][0]>1: #Cccd.noise[3][0]=ccd.noise.sn_reset.Factor
                    self.noise[3][0]=1
                    print('Sensor Simulator::: Warning! The compensation factor you entered',\
                            self.noise[3][0],' for \n the Sense Node Reset Noise cannot be more \
                            than 1! The factor is set to 1.\n')
                elif self.noise[3][0]<0:
                    self.noise[3][0]=0
                    print('Sensor Simulator::: Warning! The compensation factor you entered ',\
                            self.noise[3][0],'for the Sense Node Reset Noise cannot be negative! \
                            The factor is set to 0, SNReset noise is not simulated.')
                self.noise[3][1]=numpy.sqrt((self.Boltzman_Constant_JK)*(self.T)/(self.C_SN))
                self.noise[3][2]=numpy.exp(self.noise[3][1]*numpy.random.randn(self.sensor_size[0],\
                                                                        self.sensor_size[1] ))-1
                if self.flag[12]==1:             
                    self.Signal_CCD_voltage =(self.V_REF+self.noise[3][0]*self.noise[3][2])*\
                    (numpy.exp(-self.nonlinearity[0]*self.q*self.Signal_CCD_electrons/self.k1))
                else:
                    self.Signal_CCD_voltage=(self.V_REF+self.noise[3][0]*self.noise[3][2])\
                    -(self.Signal_CCD_electrons*self.A_SN)
            else:
                if self.flag[12]==1:
                    self.Signal_CCD_voltage=self.V_REF*(numpy.exp(-self.nonlinearity[0]*\
                                                               self.q*self.Signal_CCD_electrons/self.k1))
                else:
                    self.Signal_CCD_voltage=self.V_REF-(self.Signal_CCD_electrons*self.A_SN)
        else: #<---The sensor is CCD
            if self.flag[12]==1:
                self.Signal_CCD_voltage=self.V_REF*(numpy.exp(-self.nonlinearity[0]*self.q*\
                                                           self.Signal_CCD_electrons/self.k1))
            else: 
                self.Signal_CCD_voltage=self.V_REF-self.Signal_CCD_electrons*self.A_SN
        return self.Signal_CCD_voltage
    
    #Follow noise simulation
    def ccd_source_follower(self):
        if self.flag[13]==1: #Cccd.flag[10]=ccd.flag.VVnonlinearity
            nonlinearity_alpha = (self.A_SF*(self.nonlinearity[1]-1))/(self.V_FW)
            self.A_SF_new= nonlinearity_alpha*((self.V_REF-self.Signal_CCD_voltage)\
                                                /(self.V_REF))+(self.A_SF)*numpy.ones\
            ([self.sensor_size[0],self.sensor_size[1]])
            self.Signal_CCD_voltage=(self.Signal_CCD_voltage)*(self.A_SF_new)
        else:
            self.Signal_CCD_voltage=self.Signal_CCD_voltage*self.A_SF
        return self.Signal_CCD_voltage
    
    #ccd cds part
    def ccd_cds(self):
        if operator.eq('CMOS',self.SensorType)==0:
            if self.flag[5] == 1:  #ccd.flag[5]=ccd.flag.darkcurrent_offsetFPN
                self.noise[2][3]=self.ccd_FPN_models(self.sensor_size[0], self.sensor_size[1],\
                                                'column',self.noise[2][0], self.noise[2][1])
                self.Signal_CCD_voltage = self.Signal_CCD_voltage*(1+ self.noise[2][3]*\
                                                                    (self.V_FW*self.noise[2][2]))
        self.Signal_CCD_voltage=self.Signal_CCD_voltage*self.A_CDS
        return self.Signal_CCD_voltage
    
    #Part III Adc Transfrom part voltage to digital signal
    def ccd_adc(self):
        N_max=2**self.N_bits
        self.A_ADC=N_max/(self.V_FW-self.V_min)
        if self.flag[14]==1:    #Cccd.flag[14]=ccd.flag.ADCnonlinearity
            A_ADC_NL=self.nonlinearity[2]*self.A_ADC
            nonlinearity_alpha=(numpy.log(A_ADC_NL)/numpy.log(self.A_ADC)-1)/self.V_FW
            signal=self.V_REF-self.Signal_CCD_voltage
            A_ADC_new = (self.A_ADC)*numpy.ones([len(signal),len(signal[0])])
            A_ADC_new = A_ADC_new**(1-nonlinearity_alpha*signal)
            S_DN =numpy.round(self.S_ADC_OFFSET+A_ADC_new*signal)
        else:
            S_DN =numpy.round(self.S_ADC_OFFSET+self.A_ADC*(self.V_REF-self.Signal_CCD_voltage))
            S_DN[S_DN<=0]=0
            S_DN[S_DN>=N_max]=N_max
            self.Signal_CCD_DN=S_DN
        return self.Signal_CCD_DN
        
    #Overall function that is used to generate a figure from the un-sampled optical image
    def outputimg(self,uin):
        self.ccd_set_photosensor_constants(uin)
        if self.flag[11]!=1:   # Cccd.flag[11]=ccd.flag.darkframe
            #<-This routine for adding light noise (photon shot noise and photo response non-uniformity
            self.light_signal=self.ccd_photosensor_lightnoises(uin)
        if self.flag[2] ==1:  # Cccd.flag[2]=ccd.flag.darkcurrent
            #<-adding dark current noises that consist of Dark FPN and Dark shot noise
            self.dark_signal=self.ccd_photosensor_darknoises()           
        self.Signal_CCD_electrons = self.light_signal + self.dark_signal
        idx = (self.Signal_CCD_electrons>=self.FW_e)
        self.Signal_CCD_electrons[idx] = self.FW_e
        self.Signal_CCD_electrons=numpy.floor(self.Signal_CCD_electrons)
        self.Signal_CCD_voltage = self.ccd_sense_node_chargetovoltage()
        #<-- Signal's Voltage amplification by Source Follower
        self.Signal_CCD_voltage = self.ccd_source_follower()
        #<-- Signal's amplification and de-noising by Correlated Double Sampling
        self.Signal_CCD_voltage = self.ccd_cds()
        #<-- Analogue-To-Digital Converter
        self.Signal_CCD_DN = self.ccd_adc()
        return self.Signal_CCD_DN
