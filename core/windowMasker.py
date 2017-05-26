"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
The window masker reads in statistics and writes out a file
with windows which match the give constraints for each statistic
The masked file will then go into the window WindowSelector
Which passes on shared windows to the processor
"""
import os
from datetime import date, time, datetime, timedelta
# my classes
from spectrumReader import SpectrumReader
from statisticFrequency import StatisticFrequency 
# utils
from utilsIO import *
from utilsWindow import *

class WindowMasker(object):

	def __init__(self, proj, site, fs, decParams, winParams):
		self.proj = proj
		self.site = site
		self.fs = fs
		# stats is the stat to use
		self.stats = []
		# get the evaluation frequencies from decParams
		self.evalFreq = decParams.getEvalFrequenciesAllLevels()
		self.numLevels = decParams.getNumLevels()
		# constraints are the constraints
		# indexed by eFreq, stat and component
		# inside out is the same, but signals whether value should fall in or outside
		self.constraints = {}
		self.insideOut = {}
		# finally, a dictionary for storing the windows for each evaluation frequency
		# these are
		self.maskWindows = {}	
		# save the name of the mask
		self.maskName = ""	

	###################
	### GET FUNCTIONS
	##################
	def getMaskName(self):
		return self.maskName

	def getSite(self):
		return self.site

	def getSampleFreq(self):
		return self.fs

	def getStats(self):
		return self.stats

	def getEvalFreq(self):
		return self.evalFreq

	def getNumLevels(self):
		return self.numLevels

	def getConstraints(self):
		return self.constraints

	def getConstraintFreq(self, level, eIdx):
		eFreq = self.getEvalFreq()[level][eIdx]		
		return self.constraints[eFreq]

	def getConstraintFreqStat(self, level, eIdx, stat):
		eFreq = self.getEvalFreq()[level][eIdx]		
		return self.constraints[eFreq][stat]
		
	# return set of windows which should be removed
	# windows that do not macth constraint
	def getMaskWindows(self):
		return self.maskWindows

	# this will be used by the window selector
	# get windows for certain frequency
	def getMaskWindowsFreq(self, iDec, eIdx):
		eFreq = self.getEvalFreq()[iDec][eIdx]		
		return self.maskWindows[eFreq]

	###################
	### ADD CONSTRAINTS
	##################	
	def addStats(self, stats):
		self.stats = stats
		self.prepareDicts()

	def addConstraint(self, stat, constraint, **kwargs):
		for ilevel in xrange(0, self.getNumLevels()):
			self.addConstraintLevel(stat, constraint, ilevel, **kwargs)

	def addConstraintLevel(self, stat, constraint, level, **kwargs):
		for eIdx, eFreq in enumerate(self.getEvalFreq()[level]):
			self.addConstraintFreq(stat, constraint, level, eIdx, **kwargs)

	# don't take evaluation frequency 
	# because this sort of float matching is probably not a good ideas
	def addConstraintFreq(self, stat, constraint, level, eIdx, **kwargs):
		eFreq = self.getEvalFreq()[level][eIdx]
		insideOut = []
		if "insideOut" in kwargs:
			insideOut = kwargs["insideOut"]
		for key in constraint:
			self.constraints[eFreq][stat][key] = constraint[key]
			if key in insideOut:
				self.insideOut[eFreq][stat][key] = True
			else:
				self.insideOut[eFreq][stat][key] = False

	###################
	### PREPARE DICTIONARIES TO STORE INFO
	##################	
	def prepareDicts(self):
		for ilevel in xrange(0, self.getNumLevels()):
			for eFreq in self.evalFreq[ilevel]:
				# empty set for maskWindows to start off with
				# i.e. remove none
				self.maskWindows[eFreq] = set()
				self.constraints[eFreq] = {}
				self.insideOut[eFreq] = {}
				for stat in self.stats:
					self.constraints[eFreq][stat] = {}
					self.insideOut[eFreq][stat] = {}		

	###################
	### APPLY CONSTRAINTS
	##################	
	# the whole thing here is to find windows that fail the constraint check
	# these will be removed by the windowSelector
	def applyConstraints(self):
		# so now, run through all the statfiles listed
		rawMeas = self.proj.getSiteTimeFilesFs(self.getSite(), self.getSampleFreq())
		evalFreqs = self.getEvalFreq()					
		# loop over stats
		for stat in self.getStats():
			# loop over decimation levels
			for iDec in xrange(0, self.numLevels):	
				# create the stat handler
				# one instance can real files over and over again
				statHandler = StatisticFrequency(stat)							
				
				# loop over measurement folders
				for meas in rawMeas:

					# try and open the file
					check = statHandler.readStatFile(self.proj.getStatDataPathMeas(self.getSite(), meas), iDec)		
					if not check:
						continue # to next measurement

					numWindows = statHandler.getNumWindows()
					windowStats = statHandler.getWindowStats()
					for iW in xrange(0, numWindows):
						# now check each window
						windowVals = statHandler.getStatLocal(iW)
						# now loop over evaluation frequencies
						for eIdx, eFreq in enumerate(evalFreqs[iDec]):
							constraintCheck = True
							freqVal = windowVals[eIdx]						
							for component in self.constraints[eFreq][stat]:
								index = windowStats.index(component)
								componentVal = freqVal[index]
								test = componentVal > self.constraints[eFreq][stat][component][0] and componentVal < self.constraints[eFreq][stat][component][1]
								if self.insideOut[eFreq][stat][component]:
									test = not test
								constraintCheck = constraintCheck and test
							if not constraintCheck: # if the stat fails, add it - this is a list of windows to exclude
								self.maskWindows[eFreq].add(statHandler.getGlobalIndex(iW))

	
	#################
	### WINDOW FILES
	#################
	def writeWindowFile(self, filename):
		# write the window file
		self.maskName = filename
		infoname, winname = self.getFileNames(filename)
		infofile = open(infoname, "w")
		# first write out constraints
		infofile.write("{}\n".format(self.getSite()))
		infofile.write("{:.6f}\n".format(self.getSampleFreq()))		
		infofile.write("{}\n".format(", ".join(self.getStats())))
		evalFreq = sorted(list(self.getConstraints().keys()))
		for eFreq in evalFreq:
			infofile.write("Frequency = {:.6f}\n".format(eFreq))
			for stat in self.stats:
				infofile.write("Statistic = {}\n".format(stat))
				for component in self.constraints[eFreq][stat]:
					minVal = self.constraints[eFreq][stat][component][0]
					maxVal = self.constraints[eFreq][stat][component][1]
					infofile.write("{}\t{}\t{}\t{}\n".format(component, minVal, maxVal, self.insideOut[eFreq][stat][component]))
		# then loop through each 
		infofile.close()
		# how can I write out
		maskSize = 0
		for eFreq in evalFreq:
			if len(self.maskWindows[eFreq]) > maskSize:
				maskSize = len(self.maskWindows[eFreq])
		# create window mask array and initalise to -1
		winMaskArray = np.ones(shape=(len(evalFreq), maskSize), dtype=int)*-1
		# now fill the array
		for eIdx, eFreq in enumerate(evalFreq):
			lst = list(self.maskWindows[eFreq])
			winMaskArray[eIdx,0:len(lst)] = lst 
		np.save(winname, winMaskArray)			

	def readWindowFile(self, filename):
		self.maskName = filename
		# read the window file
		infoname, winname = self.getFileNames(filename)
		infofile = open(infoname, "r")
		lines = infofile.readlines()
		infofile.close()
		# this is all passed into the constructor
		# self.site = lines[0].strip()
		# self.fs = float(lines[1].strip())
		self.stats = lines[2].strip().split(",")
		for idx, s in enumerate(self.stats):
			self.stats[idx] = self.stats[idx].strip()
		self.prepareDicts() # prepare dictionaries to hold statistics
		lines = lines[3:]
		# now read in the statistics
		evalFreq = sorted(list(self.constraints.keys()))	
		eIdx = -1
		for l in lines:
			if "Frequency" in l:
				eIdx = eIdx + 1
			elif "Statistic" in l:
				stat = l.strip().split("=")[1]
				stat = stat.strip()
			else:
				split = l.strip().split() # component is 0, min is 1, max is 2, inout is 3
				self.constraints[evalFreq[eIdx]][stat][split[0]] = [float(split[1]), float(split[2])]
				inOut = False
				if split[3] == "True":
					inOut = True
				self.insideOut[evalFreq[eIdx]][stat][split[0]] = inOut
		# now read the window mask file
		winMaskArray = np.load(winname + ".npy")
		for eIdx, eFreq in enumerate(evalFreq):
			self.maskWindows[eFreq] = set(winMaskArray[eIdx])
			# remove -1 from the set
			self.maskWindows[eFreq] = self.maskWindows[eFreq] - set([-1])

	def getFileNames(self, filename):
		maskPath = os.path.join(self.proj.getStatDataPathSite(self.site), "masks")
		checkAndMakeDir(maskPath)
		name = filename + "_{:d}".format(int(self.fs))
		infofile = os.path.join(maskPath, name + ".info")
		winfile = os.path.join(maskPath, name) # no need for extension here, numpy adds one
		return infofile, winfile

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("WINDOW MASKER INFO BEGIN")		
		self.printText("####################")	
		self.printText("Site = {}".format(self.getSite()))	
		self.printText("Sampling frequency = {:.6f}".format(self.getSampleFreq()))	
		self.printText("Statistics to use for constraints = {}".format(", ".join(self.getStats())))	
		self.printEvalFreq()	
		self.printText("####################")
		self.printText("WINDOW MASKER INFO END")		
		self.printText("####################")

	def printEvalFreq(self):
		self.printText("Evaluation frequencies [Hz]")
		for il in xrange(0, self.getNumLevels()):
			freqForLevel = self.getEvalFreq()[il]
			eFreqStr = arrayToString(freqForLevel)
			self.printText("Level = {:d}: {}".format(il, eFreqStr))

	def printConstraints(self):
		self.printText("####################")
		self.printText("WINDOW MASKER CONSTRAINTS INFO BEGIN")		
		self.printText("####################")		
		evalFreq = sorted(list(self.getConstraints().keys()))
		for eFreq in evalFreq:
			self.printText("Frequency = {:.6f} [Hz]".format(eFreq))
			for stat in self.getStats():
				self.printText("\tStatistic = {}".format(stat))
				# check to see if any constraints
				if len(self.constraints[eFreq][stat].keys()) == 0:
					self.printText("\tNo constraints for this statistic")
					continue
				# if there are, print out
				for component in self.constraints[eFreq][stat]:
					minVal = self.constraints[eFreq][stat][component][0]
					maxVal = self.constraints[eFreq][stat][component][1]
					self.printText("\t{}\t{}\t{}\t{}".format(component, minVal, maxVal, self.insideOut[eFreq][stat][component]))
		self.printText("####################")
		self.printText("WINDOW MASKER CONSTRAINTS INFO END")		
		self.printText("####################")	

	def printMaskedWindows(self):
		self.printText("####################")
		self.printText("WINDOW MASKER WINDOWS INFO BEGIN")		
		self.printText("####################")		
		evalFreq = sorted(list(self.getConstraints().keys()))
		for eFreq in evalFreq:
			self.printText("Frequency = {:.6f} [Hz], number of masked windows = {:d}".format(eFreq, len(self.maskWindows[eFreq])))
			if len(self.maskWindows[eFreq]) == 0:
				self.printText("No masked windows")
				continue
			self.printText("{}".format(list2ranges(self.maskWindows[eFreq])))
		self.printText("####################")
		self.printText("WINDOW MASKER WINDOWS INFO END")		
		self.printText("####################")				

	def printText(self, infoStr):
		generalPrint("Window Masker", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Window Masker", warnStr)
			
