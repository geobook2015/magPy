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

# this is a statistic relating to only a window
class StatisticWindow(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self):	
		self.global2local = {}
		self.numWindows = 0
		# details about the values
		self.statsPerWindow = []
		self.numStatsPerWindow = 0
		# constraints
		self.constrainedWindows = set()
		self.constraint = []
		# win bool
		self.winBool = np.array([0])
		
	###################
	### GET STATISTIC
	##################	
	def getStatLocal(self, localIndex):
		return self.stats[localIndex, :]

	def getStatGlobal(self, globalIndex):
		return self.stats[self.global2local(globalIndex), :]
		
	def getNumWindows(self):
		return self.numWindows
		
	def getWindowStats(self):
		return self.statsPerWindow

	def getNumStatsPerWindow(self):
		return self.numStatsPerWindow
		
	def getGlobalIndices(self):
		return self.gIndices
	
	def getStats(self):
		return self.stats

	def getWindowConstraint(self, gIndex):
		return (gIndex in self.constrainedWindows)

	def getConstraint(self):
		return self.constraint
		
	# return set of windows which match constraint	
	def getConstrainedWindows(self):
		if len(self.getConstraint()) == 0:
			#return all windows
			return set(self.global2local.keys())
		return self.constrainedWindows

	# get win bool
	def getWinBool(self):
		return self.winBool
		
	###################
	### SET NUMBER OF WINDOWS
	##################	
	def setStatParams(self, numWindows, winStats):
		self.numWindows = numWindows
		# stats per window - this is a list which describes each stat
		# an example would be [Zxx, Zxy, Zyx, Zyy]
		self.winStats = winStats	
		self.numStatsPerWindow = len(winStats)	
		# set the win bool array - all true for now
		# need to set constraints to improve this
		self.winBool = np.ones(shape=(numWindows), dtype=bool)
	
	###################
	### SET CONSTRAINTS
	##################	
	def setConstraint(self, constraint, **kwargs):
		# constraint is a dictionary - only constrain on what we are interested in
		self.constraint = constraint
		# now get a set of constrained windows
		if len(constraint) != self.getNumStatsPerWindow():
			self.printWarning("Number of constraints is not sufficient - should equal the number of stats per window ({})".format(self.getNumStatsPerWindow()))
			self.printWarning("Setting all windows match the constraint")
			self.constrainedWindows = set(self.global2local.keys())
			return

		# apply any options - i.e. to use the absolute of the statistic for constraint
		statsC = np.array(self.stats)
		if "abs" in kwargs and kwargs["abs"] and kwargs["abs"]:
			statsC = np.absolute(statsC)

		# now create a set of windows that meet the constraints	
		self.winBool = np.ones(shape=(self.getNumWindows()), dtype="bool")
		for i in xrange(0, self.getNumWindows()):
			row = statsC[i,:]
			chk = True
			for j in xrange(0, self.getNumStatsPerWindow()):
				# if min and max are equal, do not use this stat for constraining
				if self.constraint[j][0] == self.constraint[j][1]:
					continue
				chk = chk and (row[j] > self.constraint[j][0] and row[j] < self.constraint[j][1]) 
			# then set true or false for the window
			self.winBool[i] = chk
		# find indices where this is true
		trueIndices = np.where(self.winBool)
		# now get the global indices
		globalIndices = self.gIndices[trueIndices]
		# and make the set
		self.constrainedWindows = set(globalIndices)
		
	###################
	### ADD STATISTIC
	##################	
	def addStat(self, localIndex, globalIndex, stat):
		data = np.zeros(shape=(self.getNumStatsPerWindow()), dtype='float')
		# do it in this order because sorted
		for idx, c in enumerate(self.getStatsPerWindow()):
			data[idx] = stat[c]
		self.global2local[globalIndex] = localIndex
		self.gIndices[localIndex] = globalIndex
		self.stats[localIndex] = data

	###################
	### READ STAT/WINDOW FILE
	###################			
	def readStatFile(self, datapath, filename, inc):
		statFile, infoFile = self.statFileName(datapath, filename, inc)
		statFile = statFile + ".npy"		
		# want a channel ordering
		if not os.path.exists(infoFile) or not os.path.exists(statFile):
			return False
		f = open(infoFile, 'r')
		lines = f.readlines()
		f.close()
		self.numWindows = int(lines[0].strip())
		self.numChannels = int(lines[1].strip())
		self.channels = lines[2].strip().split()
		self.statsPerWindow = lines[3].strip().split()
		self.numStatsPerWindow = len(self.statsPerWindow)		
		# load the statistics
		data = np.load(statFile)
		self.gIndices = np.array(data[:,0], dtype='int')
		self.stats = data[:,1:]
		for i in xrange(0, self.numWindows):
			self.global2local[self.gIndices[i]] = i	
		
		return True

	###################
	### WRITE STAT FILE
	###################
	def writeStatFile(self, datapath, filename, inc):
		# write the info file - numWindows and channel ordering
		checkAndMakeDir(datapath)
		statFile, infoFile, winFile = self.statFileName(datapath, filename, inc)
		# want a channel ordering
		f = open(infoFile, 'w')
		f.write("{:d}\n".format(self.getNumWindows()))
		f.write("{:d}\n".format(self.getNumChannels()))
		f.write("{}\n".format(" ".join(self.getChannels())))
		f.write("{}\n".format(" ".join(self.getStatsPerWindow())))			
		f.close()
		# and save binary stat file
		# concatenate gindices and stats
		data = np.hstack((self.getIndices(), self.getStats()))
		np.save(statFile, data)

	###################
	### WRITE WINDOW FILE
	###################	
	# writes a boolean file if you want to use a statistic 
	def writeWindowFile(self, datapath, filename, inc):
		# write the info file - numWindows and channel ordering
		checkAndMakeDir(datapath)
		statFile, infoFile, winFile = self.statFileName(datapath, filename, inc)
		# and save binary window file
		# concatenate gindices and stats
		data = np.hstack((self.getIndices(), self.getWinBool()))
		np.save(winFile, data)		

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
			statFile = os.path.join(datapath, "{}".format(filename))
			infoFile = os.path.join(datapath, "{}.info".format(filename))
			winFile = os.path.join(datapath, "{}.win".format(filename))
		else:
			statFile = os.path.join(datapath, "{}{:02d}".format(filename, inc))
			infoFile = os.path.join(datapath, "{}{:02d}.info".format(filename, inc))
			winFile = os.path.join(datapath, "{}{:02d}.win".format(filename, inc))					
		return statFile, infoFile, winFile	
		
	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("STATISTIC INFO BEGIN")		
		self.printText("####################")		
		self.printText("Number of windows = {:d}".format(self.getNumWindows()))	
		self.printText("Number of channels = {:d}".format(self.getNumChannels()))
		self.printText("Channels = {}".format(", ".join(self.getChannels())))
		self.printText("Number of statistics per window = {:d}".format(self.getNumStatsPerWindow()))	
		self.printText("Statistics per Window = {}".format(", ".join(self.getStatsPerWindow())))		
		constraints = self.getConstraint()
		stats = self.getStatsPerWindow()
		for i in xrange(0, len(constraints)):
			self.printText("\tConstraint: {} < {} < {}".format(constraints[i][0], stats[i], constraints[i][1]))
		self.printText("####################")
		self.printText("STATISTIC INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Statistic Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Statistic Warning", warnStr)		
		