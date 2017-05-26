"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
Given a window size and a reference time
Window parameters calculates out the 
"""
import numpy as np
# utils
from utilsIO import *
from utilsWindow import *

class WindowParams(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, decParams):
		self.decParams = decParams
		self.minSize = getDefaultMinWindowSize()
		self.minOlap = getDefaultMinOverlapSize()
		self.calcParameters(2, 4)

	###################
	### GET GENERAL INFO
	##################
	def getWindowSizeAll(self):
		return self.windows

	def getWindowSize(self, decLevel):
		return self.windows[decLevel]

	def getOverlapsAll(self):
		return self.overlaps

	def getOverlap(self, decLevel):
		return self.overlaps[decLevel]

	def getMinWindowSize(self):
		return self.minSize

	def getMinOverlap(self):
		return self.minOlap

	def getRefTime(self):
		return self.getRefTime

	###################
	### SET WINDOW PARAMETERS
	##################
	# directly set window size and parameters
	def setWindowParameters(self, windowSizes, windowOverlaps):
		if len(windowSizes) != self.decParams.getNumLevels() or len(windowOverlaps) != self.decParams.getNumLevels():
			print "Error: not enough window sizes given. Must be equal to number of decimation levels"
			return
		self.windows = windowSizes
		self.overlaps = windowOverlaps

	def setMinWindowSize(self, minSize):
		self.minSize = minSize

	def setMinOverlap(self, minOlap):
		self.minOlap = minOlap
		
	def setMinParams(self, minSize, minOlap):
		self.setMinWindowSize(minSize)
		self.setMinOverlap(minOlap)
		self.calcParameters(2, 4)

	###################
	### CALCULATE WINDOW SIZES
	##################
	def calcParameters(self, windowFactor, overlapFactor):
		# this calculates window sizes and some indexing paramters
		# with regards to reference time
		self.windows = np.ones(shape=(self.decParams.getNumLevels()), dtype=int)
		self.overlaps = np.ones(shape=(self.decParams.getNumLevels()), dtype=int)		
		decFreq = self.decParams.getDecFrequencies()
		for il in xrange(0, self.decParams.getNumLevels()):
			self.windows[il] = int(decFreq[il]/windowFactor)
			if self.windows[il] < self.minSize:
				self.windows[il] = self.minSize
			self.overlaps[il] = int(self.windows[il]/overlapFactor)
			if self.overlaps[il] < self.minOlap:
				self.overlaps[il] = self.minOlap

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("WINDOW PARAMETER INFO BEGIN")		
		self.printText("####################")		
		self.printText("Number of decimation levels = {:d}".format(self.decParams.getNumLevels()))
		decFrequencies = self.decParams.getDecFrequencies()
		for il in xrange(0, self.decParams.getNumLevels()):
			self.printText("Level = {:d}, sample freq. [Hz] = {:.6f}, sample rate [s] = {:6f}".format(
				il, decFrequencies[il], 1.0/decFrequencies[il]))		
			self.printText("Window size = {:d}, window duration [s] = {:f}".format(
				self.windows[il], (self.windows[il]-1)/decFrequencies[il]))
			self.printText("Window overlap = {:d}, overlap duration [s] = {:f}".format(
				self.overlaps[il], (self.overlaps[il]-1)/decFrequencies[il]))			
		self.printText("####################")
		self.printText("WINDOW PARAMETER INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Window Parameters Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Window Parameters Warning", warnStr)
