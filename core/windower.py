"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
The windower takes the data for a decimation level 
And a reference time to reference the window to
And then calculates number of windows
And gives access to windows
The important thing is that windows can be compared to each other
Can also choose to save a window
"""
import numpy as np
import math
from datetime import datetime, timedelta
# utils
from utilsIO import *


class Windower(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, refTime, dataTime, data, fs, winSize, winOlap):
		# data reader
		self.data = data
		self.fs = fs
		self.winSize = winSize
		self.winDuration = (winSize-1)/fs
		self.winOlap = winOlap
		self.chans = list(self.data.keys())
		self.numSamples = len(self.data[self.chans[0]])
		# refTime and dataTime are already datetime objects
		self.refTime = refTime
		self.dataTime = dataTime
		# min window warning setting
		self.minWindows = 3
		# initialise
		self.initialiseWindows()

	###################
	### GET GENERAL INFO
	##################
	def getWindowSize(self):
		return self.winSize

	def getWindowDuration(self):
		return self.winDuration

	def getOverlap(self):
		return self.winOlap

	def getNumWindows(self):
		return self.numWindows		

	def getGlobalWindowOffset(self):
		return self.winOffset

	def getFirstWindowTime(self):
		return self.firstWindowTime

	def getSampleFreq(self):
		return self.fs

	def getNumSamples(self):
		return self.numSamples

	def getRefTime(self):
		return self.refTime

	def getDataTime(self):
		return self.dataTime

	def getWindowSamples(self):
		return self.winSamples

	def getWindowTimes(self):	
		winTimes = []
		iW = 0
		for samples in self.getWindowSamples():
			start = samples[0]
			stop = samples[1]
			win = []
			# global index
			win.append(self.getGlobalWindowOffset() + iW)
			# start time
			deltaStart = timedelta(seconds=start/self.getSampleFreq())
			timeStart = self.getDataTime() + deltaStart	
			deltaEnd = timedelta(seconds=stop/self.getSampleFreq())
			timeEnd = self.getDataTime() + deltaEnd					
			# samples2end = self.getWindowSize() - 1 # need to remove initial sample
			# timeEnd = timeStart + timedelta(seconds=samples2end/self.getSampleFreq())	
			win.append(timeStart)
			win.append(timeEnd)
			winTimes.append(win)
			iW = iW + 1
		return winTimes

	def getActiveWindows(self):
		return self.winActive

	def getWindowActive(self, iW):
		return self.winActive[iW]

	def getChans(self):
		return self.chans

	def getGlobalIndex(self, iW):
		return  iWindow + self.getGlobalWindowOffset()

	def getMinWindows(self):
		return self.minWindows

	###################
	### GET DECIMATED DATA
	### FOR THIS LEVEL
	##################
	def setActiveWindows(self, winActive):
		self.winActive = winActive

	def setMinWindows(self, minWindows):
		self.minWindows = minWindows

	###################
	### GET DECIMATED DATA
	### FOR THIS LEVEL
	##################
	def getData(self, iWindow):
		winSamples = self.getWindowSamples()[iWindow]
		winData = {}
		for c in self.getChans():
			winData[c] = self.data[c][winSamples[0]:winSamples[1] + 1] # add 1 because numpy indexing like this is not inclusive
		return winData

	def getDataGlobal(self, iGlobal):
		iWindow = iGlobal - self.getGlobalWindowOffset()
		return self.getData(iWindow)

	###################
	### CALCULATE DECIMATION LEVELS
	### BASED ON AMOUNT OF DATA
	##################
	def initialiseWindows(self):
		# have a reference time
		# the first window starts there
		deltaRefStart = self.getDataTime() - self.getRefTime()
		if deltaRefStart.total_seconds() < 0:
			self.printInfo("Error: reference time is after start of recording. Stuff may go wrong!")
		# increment of window start times
		# -1 because inclusive of sample at start
		winStartIncrement = 1.0*(self.getWindowSize() - self.getOverlap())/self.getSampleFreq()
		# calculate number of windows started before reference time
		# and then by taking the ceiling, find the global index of the first window in the data
		self.winOffset = int(math.ceil(deltaRefStart.total_seconds()/winStartIncrement))
		# calculate start time of first global window
		offsetSeconds = self.winOffset*winStartIncrement
		# calculate the first window time
		self.firstWindowTime = self.getRefTime() + timedelta(seconds=offsetSeconds)
		# calculate first sample
		deltaStart = self.firstWindowTime - self.getDataTime()
		sampleStart = deltaStart.total_seconds()*self.getSampleFreq()	
		# next calculate number of windows
		# sample start is the first sample
		# window size is window size inclusive of first sample
		winStart = sampleStart
		winEnd = sampleStart + self.getWindowSize() - 1
		winStartOff = self.getWindowSize() - self.getOverlap()
		winSamples = []
		while winEnd < self.getNumSamples():
			winSamples.append([winStart, winEnd])
			winStart = winStart + winStartOff
			winEnd = winStart + self.getWindowSize() - 1
		self.numWindows = len(winSamples)
		# warn if number of windows is small
		if self.getNumWindows() < self.getMinWindows():
			self.printWarning("Number of windows in data is small - consider stopping decimation")		
		# save winSamples as numpy list in class			
		self.winSamples = np.array(winSamples, dtype=int)
		# set all windows initially to active
		self.winActive = np.ones(shape=(self.numWindows), dtype=bool)

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("WINDOWER INFO BEGIN")		
		self.printText("####################")
		self.printText("Sample freq. [Hz] = {:f}".format(self.getSampleFreq()))		
		self.printText("Window size = {:d}".format(self.getWindowSize()))
		self.printText("Overlap size = {:d}".format(self.getOverlap()))
		self.printText("Window duration [s] = {:.3f}".format(self.getWindowDuration()))
		self.printText("Reference time {}".format(self.getRefTime()))
		self.printText("Data start time {}".format(self.getDataTime()))
		self.printText("Number of complete windows in data = {:d}".format(self.getNumWindows()))
		if self.getNumWindows() < self.getMinWindows():
			self.printText("Number of windows in data is small - consider stopping decimation")
		if self.getNumWindows() > 0:		
			self.printText("Global index of first window from reference time = {}".format(self.getGlobalWindowOffset()))
			self.printText("First window starts at time {}, sample {:d}".format(self.getFirstWindowTime(), self.getWindowSamples()[0,0]))		
		self.printText("####################")
		self.printText("WINDOWER INFO END")		
		self.printText("####################")

	def printWindowTimes(self):
		winTimes = self.getWindowTimes()
		winSamples = self.getWindowSamples()
		self.printText("####################")
		self.printText("WINDOWER TIMES BEGIN")		
		self.printText("####################")
		self.printText("NOTE: Sample ranges are inclusive, to get number of samples, use: sample end - sample start + 1")
		for win, winS in zip(winTimes, winSamples):
			self.printText("Global index = {:d}, start time = {}, end time = {}, start sample = {:d}, end sample = {:d}".format(win[0], win[1], win[2], winS[0], winS[1]))
		self.printText("####################")
		self.printText("WINDOWER TIMES END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Windower Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Windower Warning", warnStr)
