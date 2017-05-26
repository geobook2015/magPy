"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
Decimation parameter calculator
Given an evaluation frequency, the following is calculated:

"""
import numpy as np
# utils
from utilsEvalFreq import *
from utilsIO import *

class DecimationParams(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, fs):
		# data reader
		self.fs = fs
		# set a variable used for calculating from frequecy
		self.divFactor = 2
		# calculate some initial values based on default vals
		# 7 levels with 7 frequencies per level
		self.calcFrequencyParams(7, 7)	

	###################
	### GET GENERAL INFO
	##################
	def getSampleFreq(self):
		return self.fs
		
	def getSampleFreqLevel(self, declevel):
		return self.fs/self.getDecFactor(declevel)
		
	def getNumLevels(self):
		return self.numLevels

	def getDecFactors(self):
		return self.decFactors

	def getDecFactor(self, decLevel):
		return self.decFactors[decLevel]

	def getIncrementalFactor(self, decLevel):
		if decLevel == 0:
			return int(self.decFactors[decLevel])
		else:
			return int(self.decFactors[decLevel]/self.decFactors[decLevel-1])

	def getDecFrequencies(self):
		return self.decFrequencies

	def getEvalFrequencies(self):
		return self.evalFreq

	def getNumFreqPerLevel(self):
		return self.freqPerLevel

	def getEvalFrequenciesForLevel(self, level):
		return self.evalFreqPerLevel[level,:]

	def getEvalFrequenciesAllLevels(self):
		return self.evalFreqPerLevel

	###################
	### GET DECIMATED DATA
	### FOR THIS LEVEL
	##################
	def getData(self):
		return self.data

	###################
	### SET FREQUENCIES 
	### PER LEVEL
	##################
	def setFrequencyParams(self, evalFreq, freqPerLevel, maxLevel):
		self.sortFreq(evalFreq)
		self.calcDecimationParams(evalFreq, freqPerLevel, maxLevel)

	def setDecimationParams(self, numLevels, freqPerLevel):
		self.calcFrequencyParams(numLevels, freqPerLevel)

	###################
	### CALCULATE DECIMATION LEVELS
	### BASED ON AMOUNT OF DATA
	##################
	def calcDecimationParams(self, evalFreq, maxLevel, freqPerLevel):
		# in case list
		evalFreq = np.array(evalFreq)
		# calculating decimation parameters from evaluation frequencies
		maxf = self.fs/4
		# find the first evaluation frequency less than or equal to maxf
		fHigh = evalFreq[0]
		for ifreq in xrange(0, evalFreq.size):
			if evalFreq[ifreq] <= maxf:
				fHigh = evalFreq[ifreq]
				break
		iHigh = evalFreq.tolist().index(fHigh)
		evalFreqSub = evalFreq[iHigh:]
		# calculate number of levels
		numLevels = maxLevel
		# check if enough evaluation frequencies
		if len(evalFreqSub) < freqPerLevel*maxLevel:
			#numLevels = int(math.floor(len(evalFreqSub)/freqPerLevel))
			numLevels = int(math.ceil(1.0*len(evalFreqSub)/freqPerLevel))

		# do another subslice
		evalFreqSub = evalFreqSub[:numLevels*freqPerLevel]

		# now create an array of evalation frequencies per decimation level
		# evalFreqPerLevel = np.ones(shape=(numLevels, freqPerLevel))
		evalFreqPerLevel = np.ones(shape=(numLevels, freqPerLevel)) * -1
		for ilevel in xrange(0, numLevels):
			for ifreq in xrange(0, freqPerLevel):
				if ilevel*freqPerLevel + ifreq >= len(evalFreqSub):
					break
				evalFreqPerLevel[ilevel, ifreq] = evalFreqSub[ilevel*freqPerLevel + ifreq]

		# now calculate decimation factors
		decFactors = np.ones(shape=(numLevels))		
		decFrequencies = np.ones(shape=(numLevels))
		for ilevel in xrange(0, numLevels):
			decFactors[ilevel], decFrequencies[ilevel] = self.calcNearestFactor(evalFreqPerLevel[ilevel][0])

		# finally, set all parameters
		self.evalFreq = evalFreqSub
		self.freqPerLevel = freqPerLevel
		self.numLevels = numLevels
		self.evalFreqPerLevel = evalFreqPerLevel
		self.decFactors = decFactors
		self.decFrequencies = decFrequencies

	def calcFrequencyParams(self, numLevels, freqPerLevel):
		# calculating evaluation frequency parameters from decimation factors
		# takes number of decimation levels
		# takes number of frequencies per level 
		# and assigns evaluation frequencies accordingly
		# along with decimation factors 
		numFreq = numLevels*freqPerLevel
		evalFreq = getEvaluationFreqSize(self.fs, numFreq)
		self.calcDecimationParams(evalFreq, numLevels, freqPerLevel)

	def sortFreq(self, freq):
		# sorted in descending order
		# sort in place
		freq[::-1].sort()

	def calcNearestFactor(self, freq):
		# want sampling frequency to be 4 times greater than highest freq
		fsMin = freq*4
		# set to initial sampling frequency
		f = float(self.fs)
		fac = 1
		while f > fsMin*self.divFactor:
			f = f/self.divFactor
			fac = fac*self.divFactor
		return fac, f

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("DECIMATION PARAMETER INFO BEGIN")		
		self.printText("####################")
		self.printText("Sampling frequency = {:f}".format(self.fs))		
		self.printText("Number of decimation levels = {:d}".format(self.getNumLevels()))
		for il in xrange(0, self.getNumLevels()):
			self.printText("Level = {:d}\tsample freq. [Hz] = {:.6f}\tsample rate [s] = {:.6f}\tdec. factor = {:07d}\tinc. factor = {:d}".format(
				il, self.decFrequencies[il], 1.0/self.decFrequencies[il], int(self.decFactors[il]), self.getIncrementalFactor(il)))
		self.printText("####################")
		self.printText("DECIMATION PARAMETER INFO END")		
		self.printText("####################")
		self.printEvalFreq()

	def printEvalFreq(self):
		self.printText("####################")
		self.printText("DECIMATION PARAMETER FREQUENCIES BEGIN")		
		self.printText("####################")
		self.printText("Evaluation frequencies [Hz]")
		for il in xrange(0, self.getNumLevels()):
			freqForLevel = self.getEvalFrequenciesForLevel(il)
			eFreqStr = arrayToString(freqForLevel)
			self.printText("Level = {:d}: {}".format(il, eFreqStr))
		self.printText("####################")
		self.printText("DECIMATION PARAMETER FREQUENCIES END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Decimation Parameters Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Decimation Parameters Warning", warnStr)
