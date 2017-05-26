"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
Univariate Statistics
Controls saving and access to various statistics
The number of statistics per window has to be provided by the user
This is better for long-term usage
"""
import os
import numpy as np
import math
# utils
from utilsIO import *

# frequency statistics
# have a single statistics for each evaluation frequency
# and for each window
# meaning, in total, there are nwindow*nfreq statistics
class StatisticFrequency(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, statName):
		self.statName = statName	
		self.numWindows = 0
		# details about the statistics
		self.winStats = []	
		self.numStatsPerWindow = 0
		# details about the evaluation frequencies
		self.evalFreq = []
		self.freq2index = {}	
		# data type
		self.dtype = "float"	
		
	###################
	### GET STATISTIC
	##################
	def getStatName(self):
		return self.statName

	def getStatLocal(self, localIndex):
		return self.stats[localIndex, :]

	def getStatGlobal(self, globalIndex):
		# find the local index of the global index
		localIndices = np.where(self.stats == globalIndex)
		if localIndices.size > 1:
			self.printWarning("There are two statistics with the same global window index")
			self.printWarning("This should never happen")
			self.printWarning("Exiting")
			exit()
		return self.stats[localIndices[0]]
		
	def getNumWindows(self):
		return self.numWindows

	def getGlobalIndex(self, iW):
		return self.globalIndices[iW]
		
	def getWindowStats(self):
		return self.winStats

	def getNumStatsPerWindow(self):
		return self.numStatsPerWindow

	def getEvalFreq(self):
		return self.evalFreq

	def getNumEvalFreq(self):
		return self.evalFreq.size
	
	def getGlobalIndices(self):
		return self.globalIndices
	
	def getStats(self, **kwargs):
		return self.stats

	def getStatsForSet(self, winSet):
			statWin = set(np.arange(self.getGlobalOffset(), self.getGlobalOffset() + self.getNumWindows()))
			winInData = winSet.intersection(statWin)
			winIndices = np.array(sorted(list(winInData))) - self.getGlobalOffset()
			return winIndices, self.stats[winIndices]		
		
	###################
	### SET NUMBER OF WINDOWS
	##################	
	def setStatParams(self, numWindows, winStats, evalFreq, **kwargs):
		self.numWindows = numWindows
		# details about the statistics
		self.winStats = winStats	
		self.numStatsPerWindow = len(winStats)
		# details about the evaluation frequencies
		self.evalFreq = evalFreq
		for idx, eFreq in enumerate(evalFreq):
			self.freq2index[eFreq] = idx		
		# data type
		if "dtype" in kwargs:
			dtype = kwargs["dtype"]
		self.stats = np.empty(shape=(self.numWindows, self.evalFreq.size, self.numStatsPerWindow), dtype=self.dtype)	
		# an array to hold global indices
		self.globalIndices = np.empty(shape=(self.numWindows), dtype=int)		
		
	###################
	### ADD STATISTIC
	##################	
	def addStat(self, localIndex, globalIndex, stat):
		# add the data
		# the data is in the format [eFreq][key][val]
		for idx, eFreq in enumerate(self.getEvalFreq()):
			# do it in this order because sorted
			for ist, st in enumerate(self.getWindowStats()):
				self.stats[localIndex, idx, ist] = stat[eFreq][st]
		# then add the global index
		self.globalIndices[localIndex] = globalIndex

	###################
	### READ STAT/WINDOW FILE
	###################			
	def readStatFile(self, datapath, inc):
		statFile, infoFile = self.statFileName(datapath, self.getStatName(), inc)
		statFile = statFile + ".npy"		
		# want a channel ordering
		if not os.path.exists(infoFile) or not os.path.exists(statFile):
			return False
		f = open(infoFile, 'r')
		lines = f.readlines()
		f.close()
		self.numWindows = int(lines[0].strip())
		self.winStats = lines[1].strip().split()
		self.numStatsPerWindow = len(self.winStats)
		self.evalFreq = np.fromstring(lines[2].strip(), dtype=float, sep=',')
		for idx, eFreq in enumerate(self.evalFreq):
			self.freq2index[eFreq] = idx	
		self.dtype = lines[3].strip()
		# now deal with the global indices
		indexInformation = lines[5:] # this is the information about local to global map
		self.globalIndices = np.empty(shape=(self.numWindows), dtype=int)
		localIndices = []
		globalIndices = []
		for inInfo in indexInformation:
			inInfo = inInfo.strip()
			if inInfo == "":
				continue
			# now get the local and global indices
			split = inInfo.split("-")
			localIndices.append(int(split[0].strip()))
			globalIndices.append(int(split[1].strip()))
		# now fill the global indices array
		for idx, localI in enumerate(localIndices):
			if idx == len(localIndices)-1:
				numWindowsFromLocal = self.numWindows - localI
				self.globalIndices[localI:] = np.arange(globalIndices[idx], globalIndices[idx] + numWindowsFromLocal)
				break
			numWindowsFromLocal = localIndices[idx+1] - localI
			self.globalIndices[localI:localI + numWindowsFromLocal] = np.arange(globalIndices[idx], globalIndices[idx] + numWindowsFromLocal)
		# load the statistics
		self.stats = np.load(statFile)
		return True

	###################
	### WRITE STAT FILE
	###################
	def writeStatFile(self, datapath, inc):
		# write the info file - numWindows and channel ordering
		checkAndMakeDir(os.path.join(datapath, "{}".format(self.getStatName())))
		statFile, infoFile = self.statFileName(datapath, self.getStatName(), inc)
		# want a channel ordering
		f = open(infoFile, 'w')
		f.write("{:d}\n".format(self.getNumWindows()))
		f.write("{}\n".format(" ".join(self.getWindowStats())))	
		# now write out the evaluation frequencies
		f.write("{}\n".format(arrayToString(self.evalFreq)))
		# finally, write out the datatype
		f.write("{}\n".format(self.dtype))
		# write out the window information
		# only write out when not consecutive
		f.write("Window map: localIndex - globalIndex\n")
		prevI = -10
		for idx, gIndex in enumerate(self.globalIndices):
			if gIndex != prevI + 1:
				f.write("{:d} - {:d}\n".format(idx, gIndex))
			prevI = gIndex
		f.close()
		# and save binary stat file
		np.save(statFile, self.getStats())	

	###################
	### GET FILE NAME
	##################
	def statFileName(self, datapath, filename, inc):
		# inc is meant to be the decimation parameter	
		# exclude extension on statfile because .npy added by numpy
		statFile = ""
		infoFile = ""
		winFile = ""
		if inc < 0:
			statFile = os.path.join(datapath, "{}".format(filename), "{}".format(filename))
			infoFile = os.path.join(datapath, "{}".format(filename), "{}.info".format(filename))
		else:
			statFile = os.path.join(datapath, "{}".format(filename), "{}{:02d}".format(filename, inc))
			infoFile = os.path.join(datapath, "{}".format(filename), "{}{:02d}.info".format(filename, inc))				
		return statFile, infoFile	
		
	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("STATISTIC INFO BEGIN")		
		self.printText("####################")	
		self.printText("Statistic Name = {}".format(self.getStatName()))	
		self.printText("Number of windows = {:d}".format(self.getNumWindows()))	
		self.printText("Statistics per Window = {}".format(", ".join(self.getWindowStats())))	
		self.printText("Evalutation frequencies = {}".format(arrayToString(self.getEvalFreq())))	
		self.printText("####################")
		self.printText("STATISTIC INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Statistic Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Statistic Warning", warnStr)		
		