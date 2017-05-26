"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
The decimator takes the atsReader and returns the decimation levels
"""
import numpy as np
import scipy.signal as signal
# utils
from utilsIO import *
from utilsProcess import *

class Decimator(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, data, fs, decParams):
		# data reader
		self.data = data
		self.fs = fs*1.0
		self.chans = list(self.data.keys())		
		self.numSamples = self.data[self.chans[0]].size 
		self.decParams = decParams
		self.minSamples = 100
		self.padWidth = 51
		self.level = -1
		self.maxDownsampleFactor = 8

	###################
	### GET INFORMATION ABOUT 
	### CURRENT STATUS
	##################
	def getCurrentLevel(self):
		return self.level

	def getSampleFreq(self):
		return self.fs

	def getNumSamples(self):
		return self.numSamples

	def getMinSamples(self):
		return self.minSamples

	def getMaxDownsampleFactor(self):
		return self.maxDownsampleFactor

	###################
	### GET DECIMATED DATA
	### FOR THIS LEVEL
	##################
	def getData(self):
		return self.data

	###################
	### SET MINIMUM NUMBER OF
	### SAMPLES FOR DECIMATED DATA
	##################			
	def setMinSamples(self, minSamples):
		self.minSamples = minSamples

	def setMaxDownsampleFactor(self, maxFactor):
		self.maxDownsampleFactor = maxFactor

	###################
	### DOWNSAMPLE TO NEXT LEVEL
	##################
	def incrementLevel(self):
		# increment level
		# 0 is the first level
		self.level = self.level+1
		downsampleFactor = self.decParams.getIncrementalFactor(self.level)

		# downsample factor should not be too high in one go
		# if greater than 8, downsample if more goes
		numDownsamples = 1
		downsampleList = [downsampleFactor]
		if downsampleFactor > self.getMaxDownsampleFactor():
			# this should give an integer
			numDownsamples = downsampleFactor/self.getMaxDownsampleFactor()
			downsampleList = [self.getMaxDownsampleFactor(), numDownsamples]
			# print info
			self.printText("Downsample factor of {:d} greater than max decimation factor {:d}.".format(downsampleFactor, self.getMaxDownsampleFactor()))
			self.printText("Downsampling in multiple decimations given by factors: {}".format(arrayToStringInt(downsampleList)))
		for iDS in xrange(0, numDownsamples):
			check = self.downsample(downsampleList[iDS])		
			if not check: # check outcome of decimation
				return False	
			# set sampling frequency and number of samples
			self.fs = self.fs/downsampleList[iDS]					
			self.numSamples = self.data[self.chans[0]].size

		return True

	###################
	### RESAMPLE FUNCTION
	##################	
	def downsample(self, downsampleFactor):
		# check to see not at max level
		if self.level >= self.decParams.getNumLevels():
			self.printWarning("Error, number of decimation levels exceeded, returning no data")
			return
		# if downsample factor is 1, nothing to do
		if downsampleFactor == 1:
			# nothing to do
			return True
		
		# need to check if it makes sense for this decimation level
		# downsampling reduce the number of samples by downsample factor
		# if new number of samples is too small, return False
		if self.getNumSamples()/downsampleFactor < self.minSamples:
			self.printWarning("Next decimation level has less than {} samples. Decimation is exiting.\nSet minimum of samples required using decimator.setMinSamples().".format(self.minSamples))
			return False

		## finally, do the downsampling
		#for c in self.chans:
		#	# decimate
		#	self.data[c] = signal.decimate(self.data[c], downsampleFactor, zero_phase=True)

		downsampleTime(self.data, downsampleFactor)
		return True


	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("DECIMATOR INFO BEGIN")		
		self.printText("####################")		
		self.printText("Current level = {:d}".format(self.getCurrentLevel()))
		if self.getCurrentLevel() == -1:
			self.printText("This is the initial level - no decimation has occured")
		self.printText("Current sample freq. [Hz] = {:.6f}".format(self.getSampleFreq()))
		self.printText("Current sample rate [s] = {:.6f}".format(1.0/self.getSampleFreq()))
		self.printText("Current number of samples = {:d}".format(self.getNumSamples()))
		self.printText("####################")
		self.printText("DECIMATOR INFO END")	
		self.printText("####################")	

	def printText(self, infoStr):
		generalPrint("Decimator Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Decimator Warning", warnStr)

