# utility functions for calculating statistics
import numpy as np
import scipy.stats as stats
import scipy.interpolate as interp
# import utils
import sys
import os
sys.path.append(os.path.join('..', 'utils'))
from utilsProcess import *
from utilsRobust import *


#################
### Basic stat info
### used in projectStatCalc
### used in projectStatPlot
#################
def statNameMap(stat, **kwargs):
	statNames = {"absvalEqn": "absvalEqn", "coherence": "cohStat", "psd": "psdStat", "poldir": "polStat", "transFunc": "tfStat", "resPhase": "resPhaseStat", "partialcoh": "pcohStat",
		"coherenceRR": "cohStat_RR", "coherenceRREqn": "cohStatEqn_RR", "absvalRREqn": "absvalEqn_RR", "transFuncRR": "tfStat_RR", "resPhaseRR": "resPhaseStat_RR"
	}
	statElements = {
		"absvalEqn": ["absExEx", "absHyEx", "absExEy", "absHyEy", "absExHx", "absHyHx", "absExHy", "absHyHy", "absEyEx", "absHxEx", "absEyEy", "absHxEy", "absEyHx", "absHxHx", "absEyHy", "absHxHy"],
		"coherence": ["cohExHx", "cohExHy", "cohEyHx", "cohEyHy"],
		"psd": ["psdEx", "psdEy", "psdHx", "psdHy"],
		"poldir": ["polExEy", "polHxHy"],
		"transFunc": ["ExHxReal", "ExHxImag", "ExHyReal", "ExHyImag", "EyHxReal", "EyHxImag", "EyHyReal", "EyHyImag"],
		"resPhase": ["ExHxRes", "ExHxPhase", "ExHyRes", "ExHyPhase", "EyHxRes", "EyHxPhase", "EyHyRes", "EyHyPhase"],
		"partialcoh": ["bivarEx", "bivarEy", "parExHx", "parExHy", "parEyHx", "parEyHy"],
		# "coherenceRR": ["ExExRR", "ExEyRR", "ExHxRR", "ExHyRR", "EyExRR", "EyEyRR", "EyHxRR", "EyHyRR", "HxExRR", "HxEyRR", "HxHxRR", "HxHyRR", "HyExRR", "HyEyRR", "HyHxRR", "HyHyRR"],
		"coherenceRR": ["ExHxRR", "ExHyRR", "EyHxRR", "EyHyRR", "HxHxRR", "HxHyRR", "HyHxRR", "HyHyRR"],
		"coherenceRREqn": ["ExHxR-HyHxR", "ExHyR-HyHyR", "EyHxR-HxHxR", "EyHyR-HxHyR"],
		"absvalRREqn": ["absHyHxR", "absExHxR", "absHyHyR", "absExHyR", "absHxHxR", "absEyHxR", "absHxHyR", "absEyHyR"], 
		"transFuncRR": ["ExHxRealRR", "ExHxImagRR", "ExHyRealRR", "ExHyImagRR", "EyHxRealRR", "EyHxImagRR", "EyHyRealRR", "EyHyImagRR"],
		"resPhaseRR": ["ExHxResRR", "ExHxPhaseRR", "ExHyResRR", "ExHyPhaseRR", "EyHxResRR", "EyHxPhaseRR", "EyHyResRR", "EyHyPhaseRR"]		
	}
	return statNames[stat], statElements[stat]

#################
### Window statistics
### These are not currently in use
#################
# signal to noise ratio is calculated as mean over std.
# this might be more useful in spectral domain
# can do this in the frequecy domain, on the amplitude
def calcSNR(specData):
	output = {}
	for c in data:
		tmp = np.absolute(data[c])
		output[c] = np.average(tmp)/np.std(tmp)
	return output

# The Pearson correlation coefficient measures the linear relationship
# between two datasets. Strictly speaking, Pearson's correlation requires
# that each dataset be normally distributed.
def pearsonCoefficient(data):
	# construct input matrices
	# this needs to be columns for observations, rows for variables
	# and this needs to be output as a dictionary
	chans = sorted(list(data.keys()))
	numChans = len(chans) # the channels 
	output = {}
	for i in xrange(0, numChans):
		for j in xrange(i, numChans):
			key = "{}{}".format(chans[i], chans[j])
			# now calculate the pearson correlation coefficient
			pcc, pval = stats.pearsonr(data[chans[i]], data[chans[j]])			
			output[key] = pcc
	return output

def pearsonCoefficientSpec(data):
	# calcaulates the PCC for magnitude and phase separately
	mag = {}
	phase = {}
	for c in data:
		mag[c] = np.absolute(data[c])
		phase[c] = np.unwrap(np.angle(data[c]))
	magPCC = pearsonCoefficient(mag)
	phasePCC = pearsonCoefficient(phase)
	return magPCC, phasePCC		
				



		
