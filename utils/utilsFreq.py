# utility functions for frequency related stuff
import numpy as np
import numpy.fft as fft
import math

def getFrequencyArray(fs, samples):
	# frequencies go from to nyquist		
	nyquist = fs/2
	return np.linspace(0, nyquist, samples)

# use this function for all FFT calculations
# then if change FFT later (i.e. FFTW), just replace one function
def forwardFFT(data, **kwargs):
	if "norm" in kwargs and not kwargs["norm"]:
		return fft.rfft(data, axis=0)
	return fft.rfft(data, norm='ortho', axis=0)

def inverseFFT(data, length, **kwargs):
	if "norm" in kwargs and not kwargs["norm"]:
		return fft.irfft(data, n=length)
	return fft.irfft(data, n=length, norm='ortho')

def padNextPower2(size):
	next2Power = math.ceil(math.log(size,2))
	next2Size = math.pow(2, int(next2Power))
	return int(next2Size) - size


