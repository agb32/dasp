"""power spectral density models, and generation of random time series from these"""
import numpy

def psd1():
    psd=(numpy.arange(1000)**2+2**2)**(-2.)
    return psd

def timeSeries(psd=None):
    if psd is None:
        psd=psd1()
    a=numpy.sqrt(2*psd)
    phi=numpy.random.random(a.size)*2*numpy.pi
    z=a*numpy.exp(1j*phi)
    seq=numpy.fft.ifft(z)
    return seq
