# get evaluation frequencies
import numpy as np
import math

def getEvaluationFreq(fs, minFreq):
	# highest frequency is Nyquist/4
	# lowest frequency is user definable
	# this uses
	# f_i = f_max / pow(2,(i-1)/2)
	fmax = fs/4
	freq = []
	i = 1
	f = fmax
	while f > minFreq:
		f = fmax/math.pow(2,(i-1.0)/2.0)
		freq.append(f)
		i = i + 1
	return np.array(freq)

def getEvaluationFreqSize(fs, numFreq):
	# highest frequency is Nyquist/4
	# lowest frequency is user definable
	# this uses
	# f_i = f_max / pow(2,(i-1)/2)
	fmax = fs/4
	freq = []
	i = 1
	f = fmax
	while i <= numFreq:
		f = fmax/math.pow(2,(i-1.0)/2.0)
		freq.append(f)
		i = i + 1
	return np.array(freq)	

def getCustomEvaluationFreq():
	# this will come from elsewhere
	# for now, simply return above
	return getEvaluationFreq()