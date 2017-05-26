"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
The window selector calculates which global windows to use
Calculates overlapping windows between given sites
This removes the burden on the upcoming processor
Upcoming processor can then simply get the data for the windows
And process them
"""
import os
from datetime import date, time, datetime, timedelta
# my classes
from spectrumReader import SpectrumReader
from windowMasker import WindowMasker
# utils
from utilsIO import *
from utilsWindow import *

class WindowSelector(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, proj, fs, decParams, winParams):
		self.proj = proj
		self.fs = fs 
		self.decParams = decParams
		self.winParams = winParams
		self.sites = []
		# shared indices
		self.sharedIndices = {}
		# the masks to use for each site
		self.siteMasks = {}		
		# weights to use for each site
		self.siteWeights = {}
		# the spec files for eac site at fs
		self.siteSpecFolders = {}		
		self.siteSpecReaders = {}
		# global index ranges for all of the spec files
		self.siteSpecRanges = {}
		# set of all global indices for each site
		self.siteGlobalIndices = {}
		self.prepend = "spectra"
		# time constraints
		# priority is datetimes
		# then dates
		# then times
		self.datetimes = {}		
		self.dates = {}
		self.times = {}
		# final constraints saved in
		self.datetimeConstraints = {}
		# add a list for each decimation level		
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):
			self.datetimes[iDec] = []
			self.dates[iDec] = []
			self.times[iDec] = []
			self.datetimeConstraints[iDec] = []			


	###################
	### GET FUNCTIONS
	##################
	def getSites(self):
		return self.sites

	def getSharedWindows(self):
		return self.sharedIndices

	def getSharedWindowsLevel(self, iDec):
		return self.sharedIndices[iDec]

	def getNumSharedWindows(self, iDec):
		return len(self.sharedIndices[iDec])		

	def getWindowsForFreq(self, iDec, eIdx):
		sharedIndices = self.getSharedWindowsLevel(iDec)
		# now mask for the particular frequency - mask for each given site
		for s in self.getSites():
			for mask in self.getMasks()[s]:
				# remove the masked windows from shared indices
				sharedIndices = sharedIndices - mask.getMaskWindowsFreq(iDec, eIdx)	
		return sharedIndices

	# do other helper, which calculates the number of non masked windows for the whole level
	# this should significantly speed up calculation when constraints are applied
	def getUnmaskedWindowsLevel(self, iDec):
		indices = set()
		evalFreq = self.getDecParams().getEvalFrequenciesForLevel(iDec)
		for eIdx, eFreq in enumerate(evalFreq):
			indices.update(self.getWindowsForFreq(iDec, eIdx))
		return indices

	def getSpecReaders(self):
		return self.siteSpecReaders				

	def getSpecRanges(self):
		return self.siteSpecRanges

	def getGlobalIndices(self):
		return self.siteGlobalIndices		

	def getSampleFreq(self):
		return self.fs

	def getPrepend(self):
		return self.prepend

	def getDecParams(self):
		return self.decParams

	def getWinParams(self):
		return self.winParams		

	def getDatetimeConstraints(self):
		self.calcDatetimeConstraints()
		return self.datetimeConstraints	

	def getLevelDatetimeConstraints(self, iDec):
		self.calcDatetimeConstraints()
		return self.datetimeConstraints[iDec]

	def getMasks(self):
		return self.siteMasks

	def getSpecReaderForWindow(self, site, iDec, iWin):
		specRanges = self.getSpecRanges()[site][iDec]
		specReaders = self.getSpecReaders()[site][iDec]
		for sF in specRanges:
			if iWin >= specRanges[sF][0] and iWin <= specRanges[sF][1]:
				return sF, specReaders[sF]

		# if here, no window found
		self.printWarning("Shared window {}, decimation level {} does not appear in any files given the constraints applied".format(iWin, iDec))
		return False, False

	def getDataSize(self, iDec):
		# return data size of first file
		dataSize = -1
		site = self.getSites()[0]
		specReaders = self.getSpecReaders()[site][iDec]		
		for sF in specReaders:
			return specReaders[sF].getDataSize()

	###################
	### SET FUNCTIONS
	##################
	def setSites(self, sites):
		# first remove repeated sites
		sitesSet = set(sites)
		sites = list(sitesSet)
		# now continue
		self.sites = sites	
		for s in self.sites:
			self.siteMasks[s] = []			
			self.siteSpecFolders[s] = []
			self.siteSpecReaders[s] = {}			
			self.siteSpecRanges[s] = {}
			# use sets to hold gIndices
			# optimised to find intersections
			self.siteGlobalIndices[s] = {}
		# at the end, calculate global indices
		self.calcGlobalIndices()

	# this is the prepend for the spectra files
	def setPrepend(prepend):
		self.prepend = prepend	

	###################
	### ADD CONSTRAINTS
	##################					
	# for datetime constrains, dates take priority
	def addDatetimeConstraint(self, start, stop):
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):	
			self.addLevelDatetimeConstraint(start, stop, iDec)

	def addLevelDatetimeConstraint(self, start, stop, iDec):
		datetimeStart = datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
		datetimeStop = datetime.strptime(stop, '%Y-%m-%d %H:%M:%S')		
		self.datetimes[iDec].append([datetimeStart, datetimeStop])		

	def addDateConstraint(self, dateC):
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):	
			self.addLevelDateConstraint(dateC, iDec)

	def addLevelDateConstraint(self, dateC, iDec):
		datetimeC = datetime.strptime(dateC, '%Y-%m-%d').date()		
		self.dates[iDec].append(datetimeC)		

	def addTimeConstraint(self, start, stop):
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):	
			self.addLevelTimeConstraint(start, stop, iDec)

	def addLevelTimeConstraint(self, start, stop, iDec):
		timeStart = datetime.strptime(start, '%H:%M:%S').time()
		timeStop = datetime.strptime(stop, '%H:%M:%S').time()			
		self.times[iDec].append([timeStart, timeStop])		

	# this is a mask for with values for each evaluation frequency
	def addWindowMask(self, site, maskName, **kwargs): 
		winMasker = WindowMasker(self.proj, site, self.getSampleFreq(), self.getDecParams(), self.getWinParams())
		winMasker.readWindowFile(maskName)
		self.siteMasks[site].append(winMasker)

	###################
	### GET SHARED GLOBAL WINDOWS
	### THIS DOES NOT INCLUDE ANY MASKS WHICH MIGHT BE APPLIED
	##################		
	def calcSharedWindows(self):
		if len(self.getSites()) == 0:
			self.printWarning("No sites given to Window Selector. At least one site needs to be given.")
			return False

		# calculate datetime constraints
		self.calcDatetimeConstraints()
		
		# initialise the sharedIndices with a set from one site
		sites = self.getSites()
		siteInit = sites[0]
		numLevels = self.getDecParams().getNumLevels()
		for iDec in xrange(0, numLevels):
			self.sharedIndices[iDec] = self.getGlobalIndices()[siteInit][iDec]
			
		# now for each decimation level
		# calculate the shared ones
		for iDec in xrange(0, numLevels):			
			for s in self.getSites():
				self.sharedIndices[iDec] = self.sharedIndices[iDec].intersection(self.getGlobalIndices()[s][iDec])

		# apply time constraints
		# time constraints should be formulate as a set
		# and then, find the intersection again
		for iDec in xrange(0, numLevels):
			constraints = self.getLevelDatetimeConstraints(iDec)
			if len(constraints) != 0:				
				datetimeIndices = set()
				for dC in constraints:
					gIndexStart, firstWindowStart = datetime2gIndex(self.proj.getRefTime(), dC[0], self.decParams.getSampleFreqLevel(iDec), self.winParams.getWindowSize(iDec), self.winParams.getOverlap(iDec))
					gIndexEnd, firstWindowEnd = datetime2gIndex(self.proj.getRefTime(), dC[1], self.decParams.getSampleFreqLevel(iDec), self.winParams.getWindowSize(iDec), self.winParams.getOverlap(iDec))
					gIndexEnd = gIndexEnd - 1 # as the function returns the next window starting after time
					if gIndexEnd < gIndexStart:
						gIndexEnd = gIndexStart
					datetimeIndices.update(range(gIndexStart, gIndexEnd))
					self.printText("Decimation level = {}. Applying date constraint {} - {}, global index constraint {} - {}".format(iDec, dC[0], dC[1], gIndexStart, gIndexEnd))
				self.sharedIndices[iDec] = self.sharedIndices[iDec].intersection(datetimeIndices)	
		

	###################
	### GET WINDOW RANGES
	##################	
	def calcGlobalIndices(self):
		# get all the spectra files with the correct sampling frequency
		for s in self.getSites():
			timeFilesFs = self.proj.getSiteTimeFilesFs(s, self.getSampleFreq())
			specFiles = self.proj.getSiteSpectraFiles(s)
			specFilesFs = []
			for sF in specFiles:
				if sF in timeFilesFs:
					specFilesFs.append(sF)

			self.siteSpecFolders[s] = specFilesFs
			# for each decimation level			
			# loop through each of the spectra folders
			# and find the global indices ranges for each decimation level
			numLevels = self.decParams.getNumLevels()
			for iDec in xrange(0, numLevels):
				# get the dictionaries ready
				self.siteSpecReaders[s][iDec] = {}
				self.siteSpecRanges[s][iDec] = {}
				self.siteGlobalIndices[s][iDec] = set()	
				# loop through spectra folders and figure out global indices							
				for sF in self.siteSpecFolders[s]:
					specReader = SpectrumReader(os.path.join(self.proj.getSpecDataPathSite(s), sF))
					check = specReader.openBinaryForReading(self.getPrepend(), iDec)
					# see if file exists
					# if not, continue
					if not check:
						continue	
					self.siteSpecReaders[s][iDec][sF] = specReader
					globalRange = specReader.getGlobalRange()
					self.siteSpecRanges[s][iDec][sF] = globalRange
					# and save set of global indices
					self.siteGlobalIndices[s][iDec].update(range(globalRange[0], globalRange[1]+1))

	# Datetime constraints: priority is datetime, then dates, then times
	def calcDatetimeConstraints(self):
		# calculate site dates if required
		siteDates = self.calcSiteDates()

		# datetime constraints are for each decimation level
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):
			# calculate date and time constraints for each level
			# begin with the datetime constraints - these have highest priority			
			self.datetimeConstraints[iDec] = self.datetimes[iDec]	
			
			# check to see whether any date and time constraints
			if len(self.dates[iDec]) == 0 and len(self.times[iDec]) == 0:
				continue

			dateConstraints = []
			if len(self.dates[iDec]) != 0:
				# apply time constraints only on specified days
				dateConstraints = self.dates[iDec]
			else:
				dateConstraints = siteDates

			# finally, add the time constraints to the dates
			# otherwise add the whole day
			dateAndTimeConstraints = []
			if len(self.times[iDec]) == 0:
				# add whole days
				for dC in dateConstraints:
					start = datetime.combine(dC, time(0,0,0))
					stop = datetime.combine(dC, time(23,59,59))
					dateAndTimeConstraints.append([start, stop])
			else:
				# add each time for each day
				for tC in self.times[iDec]:
					for dC in dateConstraints:
						start = datetime.combine(dC, tC[0])
						stop = datetime.combine(dC, tC[1])
						# check if this goes over a day
						if tC[1] < tC[0]:
							# then need to increment the day
							dCNext = dC + timedelta(days=1)
							stop = datetime.combine(dCNext, tC[1])
						# append to constraints
						dateAndTimeConstraints.append([start, stop])	

			# finally, combine datetimes and dateAndTimeConstraints
			self.datetimeConstraints[iDec] = self.datetimeConstraints[iDec] + dateAndTimeConstraints	
			self.datetimeConstraints[iDec] = sorted(self.datetimeConstraints[iDec])	

	def calcSiteDates(self):
		starts = []
		stops = []
		for s in self.getSites():
			starts.append(self.proj.getSiteStart(s))
			stops.append(self.proj.getSiteStop(s))	
		# need all the dates between
		d1 = max(starts).date()
		d2 = min(stops).date()
		if d1 > d2:
			self.printWarning("A site passed to the window selector does not overlap with any other sites. There will be no shared windows")
			return
		# now with d2 > d1
		siteDates = []
		delta = d2 - d1
		# + 1 because inclusive of stop and start days
		for i in range(delta.days + 1):
			siteDates.append(d1 + timedelta(days=i))

		return siteDates

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("WINDOW SELECTOR INFO BEGIN")		
		self.printText("####################")		
		self.printText("Sampling frequency [Hz] = {:.6f}".format(self.getSampleFreq()))
		self.printText("Sites = {}".format(", ".join(self.getSites())))
		self.printText("####################")
		self.printText("WINDOW SELECTOR INFO END")		
		self.printText("####################")

	def printAllSiteInfo(self):
		for s in self.getSites():
			self.printSiteInfo(s)

	def printSiteInfo(self, site):
		self.printText("####################")
		self.printText("WINDOW SELECTOR SITE INFO BEGIN")		
		self.printText("####################")		
		self.printText("Sampling frequency [Hz] = {:.6f}".format(self.getSampleFreq()))
		self.printText("Site = {}".format(site))
		self.printText("Site global index information")
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):	
			self.printText("Decimation Level = {:d}".format(iDec))
			ranges = self.getSpecRanges()
			for sF in sorted(list(ranges[site][iDec].keys())):
				startTime1, endTime1 = gIndex2datetime(ranges[site][iDec][sF][0], self.proj.getRefTime(), self.getSampleFreq()/self.decParams.getDecFactor(iDec), self.winParams.getWindowSize(iDec), self.winParams.getOverlap(iDec))
				startTime2, endTime2 = gIndex2datetime(ranges[site][iDec][sF][1], self.proj.getRefTime(), self.getSampleFreq()/self.decParams.getDecFactor(iDec), self.winParams.getWindowSize(iDec), self.winParams.getOverlap(iDec))				
				self.printText(
					"Measurement file = {}\ttime range = {} - {}\tGlobal Indices Range = {:d} - {:d}".format(
						sF, startTime1, endTime2, ranges[site][iDec][sF][0], ranges[site][iDec][sF][1]
					)
				)
		self.printText("####################")
		self.printText("WINDOW SELECTOR SITE INFO END")		
		self.printText("####################")		

	def printSharedIndices(self):
		self.printText("####################")
		self.printText("WINDOW SELECTOR SHARED INDICES INFO BEGIN")					
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):
			self.printText("####################")				
			self.printText("Decimation Level = {:d}".format(iDec))	
			self.printText("Number of shared windows = {:d}".format(self.getNumSharedWindows(iDec)))				
			self.printText("Shared Window Indices: {}".format(list2ranges(self.getSharedWindows()[iDec])))
			self.printText("NOTE: These are the shared windows at each decimation level. Windows for each evaluation frequency might vary depending on masks")
		self.printText("####################")
		self.printText("WINDOW SELECTOR SHARED INDICES INFO END")		
		self.printText("####################")		

	def printDatetimeConstraints(self):
		# calculate datetime constraints
		self.calcDatetimeConstraints()
		# print out the info
		self.printText("####################")
		self.printText("WINDOW SELECTOR CONSTRAINT INFO BEGIN")		
		self.printText("####################")	
		self.printText("Datetime constraints")
		numLevels = self.decParams.getNumLevels()
		for iDec in xrange(0, numLevels):				
			self.printText("Decimation Level = {:d}".format(iDec))	
			for d in self.getLevelDatetimeConstraints(iDec):
				self.printText("Constraint {} - {}".format(d[0], d[1]))	
		self.printText("####################")
		self.printText("WINDOW SELECTOR CONSTRAINT INFO END")		
		self.printText("####################")	

	def printWindowMasks(self):
		self.printText("####################")
		self.printText("WINDOW SELECTOR MASK INFO BEGIN")		
		self.printText("####################")		
		for s in self.getSites():
			self.printText("Site = {}".format(s))
			if len(self.getMasks()[s]) == 0:
				self.printText("\tNo masks for this site")
			else:
				for mask in self.getMasks()[s]:
					self.printText("\tMask = {}".format(mask.getMaskName()))
		self.printText("####################")
		self.printText("WINDOW SELECTOR MASK INFO END")		
		self.printText("####################")	

	def printWindowsForFrequency(self):
		self.printText("####################")
		self.printText("WINDOW SELECTOR FREQUENCY WINDOWS INFO BEGIN")		
		self.printText("####################")		
		for iDec in xrange(0, self.getDecParams().getNumLevels()):
			evalFreq = self.getDecParams().getEvalFrequenciesForLevel(iDec)
			unmaskedWindows = self.getNumSharedWindows(iDec)
			for eIdx, eFreq in enumerate(evalFreq):
				maskedWindows = self.getWindowsForFreq(iDec, eIdx)
				self.printText("Evaluation frequency = {:.6f}, shared windows = {:d}, windows after masking = {:d}".format(eFreq, unmaskedWindows, len(maskedWindows)))
				self.printText("{}".format(list2ranges(maskedWindows)))
		self.printText("####################")
		self.printText("WINDOW SELECTOR FREQUENCY WINDOWS INFO END")		
		self.printText("####################")	

	def printText(self, infoStr):
		generalPrint("Window Selector Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Window Selector Warning", warnStr)


