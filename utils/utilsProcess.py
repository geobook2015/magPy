# utility functions for frequency related stuff
import numpy as np
import math
import scipy.signal as signal
import scipy.interpolate as interp
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
# import utils
from utilsIO import *

###########################
### SMOOTHING FUNCTIONS
###########################
def smooth1d(x, winLen=11, window='hanning'):
	"""
	numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
	scipy.signal.lfilter
	
	TODO: the window parameter could be the window itself if an array instead of a string
	NOTE: length(output) != length(input), to correct this: return y[(winLen/2-1):-(winLen/2)] instead of just y.
	""" 
	
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	
	if x.size < winLen:
		raise ValueError, "Input vector needs to be bigger than window size."
		
	if winLen<3:
		return x
	
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman', 'cosine', 'parzen']:
		raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	
	s = np.pad(x, (winLen, winLen), mode="edge")
	if window == 'flat': #moving average
		w=np.ones(winLen, 'd')
	else:
		w=eval('signal.'+window+'(winLen)')
	
	off = winLen + (winLen - 1)/2
	if winLen%2 == 0:
		off = winLen + (winLen/2)
	
	y = np.convolve(s, w/w.sum(), mode='full')
	return y[off:off+x.size] 
	
def smooth2d(x, winLen=11, window='hanning'):
	# winLen[0] is smoothing along windows
	# winLen[1] is smoothing in a single window
	kernel = np.outer(signal.hanning(winLen[0], 8), signal.gaussian(winLen[1], 8))
	# pad to help the boundaries
	padded = np.pad(x, ((winLen[0], winLen[0]),(winLen[1], winLen[1])), mode="edge")
	# 2d smoothing
	blurred = signal.fftconvolve(padded, kernel, mode='same')
	return blurred[winLen[0]:winLen[0]+x.shape[0], winLen[1]:winLen[1]+x.shape[1]]	


###########################
### STANDARD FILTERING
### AND DECIMATION
###########################
def downsampleTime(data, downsampleFactor):
	# if downsampleFactor is 1, nothing to do
	if downsampleFactor == 1:
		return data

	# a downsample factor should not be greater than 13
	# hence factorise downsampleFactors that are greater than this
	if downsampleFactor > 13:
		downsamples = factorise(downsampleFactor)
		generalPrint("Decimation", "Downsample factor {} greater than 13. Downsampling will be performed in multiple steps of {}".format(
			downsampleFactor, arrayToStringInt(downsamples)))
	else:
		downsamples = [downsampleFactor]

	# downsample for each factor in downsamples
	for factor in downsamples:
		for c in data:
			data[c] = signal.decimate(data[c], factor, zero_phase=True)
	return data

def factorise(number):
	import primefac
	factors = list(primefac.primefac(number))
	downsamples = []
	# there's a few pathological cases here that are being ignored
	# what if the downsample factor is the product of two primes greater than 13
	# let's ignore this situation for the time being
	val = 1
	for f in factors:
		test = val*f
		if test > 13:
			downsamples.append(val)
			val = 1
		val = val*f
	# logic: on the last value of f, val*f is tested
	# if this is greater than 13, the previous val is added, which leaves one factor leftover
	# if not greater than 13, then this is not added either
	# so append the last value. the only situation in which this fails is the last factor itself is over 13.
	downsamples.append(val)
	return downsamples


def resample(data, fs, fsNew):
	# resample the data using the polyphase method which does not assume periodicity
	# need to calculate the upsample and then the downsample
	# using polyphase filtering, the final sample rate = up / down * original sample rate
	# need to calculate up and down
	from fractions import Fraction
	frac = Fraction(1.0/fs).limit_denominator() # because this is most probably a float
	frac = Fraction(frac*int(fsNew))
	frac.limit_denominator()
	# now do the resampling
	# if frac.numerator == 1:
	# 	# then decimate instead of resample
	# 	return downsampleTime(data, frac.denominator)

	# otherwise, normal polyphase filtering	
	resampleData = {}
	for c in data:
		resampleData[c] = signal.resample_poly(data[c], frac.numerator, frac.denominator)
	return resampleData

# lowpass butterworth filter
def lpFilter(data, fs, cutoff, order=5):
	# create the filter
	normalisedCutoff = 2.0*cutoff/fs
	b, a = signal.butter(order, normalisedCutoff, btype="lowpass", analog=False)
	# filter each channel
	return filterData(data, b, a)

# highpass butterworth filter
def hpFilter(data, fs, cutoff, order=5):
	# create the filter
	normalisedCutoff = 2.0*cutoff/fs
	b, a = signal.butter(order, normalisedCutoff, btype="highpass", analog=False)
	return filterData(data, b, a)

def bpFilter(data, fs, cutoffLow, cutoffHigh, order=5):
	# create the filter
	normalisedCutoffLow = 2.0*cutoffLow/fs
	normalisedCutoffHigh = 2.0*cutoffHigh/fs
	b, a = signal.butter(order, [normalisedCutoffLow, normalisedCutoffHigh], btype="bandpass", analog=False)
	return filterData(data, b, a)

def filterData(data, b, a, padLen=10000):
	# filter each channel
	filteredData = {}
	for c in data:
		# filteredData[c] = signal.filtfilt(b, a, data[c], method="pad", padtype="odd", padlen=padLen)
		filteredData[c] = signal.filtfilt(b, a, data[c], method="gust", irlen=500)
	return filteredData

# Required input defintions are as follows;
# time:   Time between samples
# band:   The bandwidth around the centerline freqency that you wish to filter
# freq:   The centerline frequency to be filtered
# ripple: The maximum passband ripple that is allowed in db
# order:  The filter order.  For FIR notch filters this is best set to 2 or 3,
#         IIR filters are best suited for high values of order.  This algorithm
#         is hard coded to FIR filters
# filter_type: 'butter', 'bessel', 'cheby1', 'cheby2', 'ellip'
# data:         the data to be filtered
def notchFilter(data, fs, freq, band):
	nyq  = fs/2.0
	low  = freq - band/2.0
	high = freq + band/2.0
	low  = low/nyq
	high = high/nyq

	# some options
	#ripple = 
	order = 2
	filter_type = "bessel"

	#b, a = signal.iirfilter(order, [low, high], rp=ripple, btype='bandstop', analog=False, ftype=filter_type)
	b, a = signal.iirfilter(order, [low, high], btype='bandstop', analog=False, ftype=filter_type)
	filteredData = signal.lfilter(b, a, data)
	return filteredData

###########################
### SHIFTING
###########################
# this actually shifts the time by some amount of time
def timeShift(data, fs, shift, **kwargs):
	shiftMode = "samples"
	if "mode" in kwargs:
		shiftMode = kwargs["mode"]

	# calculate out shift

	# apply shift
	return data

# need a function to interpolate the sampling so that it coincides with full seconds
# the function also shifts the start point to the next full second
# TODO: this function needs to be more robust for low (< 1Hz) sample frequencies as the use of microseconds and seconds makes no sense for this
# THIS FUNCTION WILL TRUNCATE THE DATA TO THE NEXT SECOND
def interpolateToSecond(fs, startTime, data):
	# data properties
	chans = data.keys()
	samplePeriod = 1.0/fs
	# set initial vals
	numSamples = data[chans[0]].size

	# now caluclate the interpolation
	microseconds = startTime.time().microsecond
	# check if the dataset already begins on a second
	if microseconds == 0:
		return startTime, numSamples, data # do nothing, already on the second

	# now turn microseconds into a decimal
	microseconds = microseconds/1000000.0
	# now calculate the number of complete samples till the next second
	eps = 0.000000001
	test = microseconds
	samplesToDrop = 0	
	# this loop will always either calculate till the full second or the next sample passed the full second		
	while test < 1.0 - eps:
		test += samplePeriod
		samplesToDrop += 1

	# if this is exact, i.e. integer number of samples to next second, just need to drop samples
	multiple = (1.0 - microseconds)/samplePeriod
	if np.absolute(multiple - samplesToDrop) < eps: # floating point arithmetic
		dataInterp = {} # create a new dictionary for data
		for chan in chans:
			dataInterp[chan] = data[chan][samplesToDrop:]		
		# update the other data
		numSamplesInterp = numSamples - samplesToDrop
		startTimeInterp = startTime + timedelta(seconds=1.0*samplesToDrop/fs)	
		return startTimeInterp, numSamplesInterp, dataInterp	
	
	# if here, then we have calculated one extra for samplesToDrop
	samplesToDrop -= 1 

	# now the number of samples to the next full second is not an integer
	# interpolation will have to be performed
	shift = (multiple - samplesToDrop)*samplePeriod
	sampleShift = shift/samplePeriod
	x = np.arange(0, numSamples)
	xInterp = np.arange(samplesToDrop, numSamples - 1) + sampleShift
	# calculate return vars
	numSamplesInterp = xInterp.size	
	startTimeInterp = startTime + timedelta(seconds=1.0*samplesToDrop/fs) + timedelta(seconds=shift)

	# do the interpolation
	dataInterp = {}
	for chan in chans:
		interpFunc = interp.InterpolatedUnivariateSpline(x, data[chan])
		dataInterp[chan] = interpFunc(xInterp)

	# need to calculate how much the 
	return startTimeInterp, numSamplesInterp, dataInterp

###########################
### REMOVE BAD DATA
###########################
def removeZeros(data):
	# this function finds a stretch of zeros and tries to fill them in with better data
	# i.e. interpolated data or some such
	# first, identify zeros
	for chan in data:
		data[chan] = removeZerosSingle(data[chan])
	return data

def removeZerosSingle(data):
	eps = 0.000000001 # use this because of floating point precision	
	# set an x array
	x = np.arange(data.size)
	# find zero locations
	zeroLocs = np.where(np.absolute(data) < eps)[0] # this returns a tuple, take the first index
	if len(zeroLocs) == 0:
		return data # no zeros to remove

	# now want to find consecutive zeros
	grouped = groupConsecutive(zeroLocs)
	indicesToFix = []
	# now find groups of 3+
	for g in grouped:
		if g.size >= 20:
			indicesToFix = indicesToFix + list(g)
	# now have the indices we want to fix
	# can go about interpolating values there
	indicesToFix = np.array(sorted(indicesToFix))
	mask = np.ones(data.size, np.bool)
	mask[indicesToFix] = 0
	data[indicesToFix] = np.interp(indicesToFix, x[mask], data[mask])
	return data

def removeNans(data):
	# find nan in the dataset and removes the values
	for chan in data:
		data[chan] = removeNansSingle(data[chan])
	return data

def removeNansSingle(data):
	# set an x array
	x = np.arange(data.size)
	# find locations of nans - this is a bool array with True in locations with nan values
	nanLocs = np.isnan(data)
	# if no nans, do nothing
	if not np.any(nanLocs):
		return data # no nans to remove

	# create mask
	mask = np.ones(data.size, np.bool)
	mask[nanLocs] = 0 # using numpy indexing with bool arrays	
	# no need to group, want to remove every nan
	data[nanLocs] = np.interp(x[nanLocs], x[mask], data[mask])
	return data

def groupConsecutive(vals, stepsize=1):
    """Return list of consecutive lists of numbers from vals (number list)."""
    return np.split(vals, np.where(np.diff(vals) != stepsize)[0]+1)








