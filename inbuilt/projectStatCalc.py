#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('..', 'core'))
sys.path.append(os.path.join('..', 'stats'))
sys.path.append(os.path.join('..', 'utils'))
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
# import my classes
from project import Project
from decimationParameters import DecimationParams
from windowParameters import WindowParams
from spectrumReader import SpectrumReader
from windowSelector import WindowSelector
from statisticCalculator import StatisticCalculator
from statisticFrequency import StatisticFrequency
# utils
from utilsWindow import *
from utilsIO import *
from utilsStats import *

# calculate various statistics
# and save for either plotting or thresholding
def projectStatCalc(proj, **kwargs):
	generalPrint("ProjectStatCalc", "Calculating project statistics with options: {}".format(kwargs))
	# default options
	options = parseKeywords(getDefaultOptions(proj), kwargs)			
	# get the reference time
	datetimeRef = proj.getRefTime()	
	# create the statistic calculator
	statCalculator = StatisticCalculator()

	# loop over sites
	for s in options["sites"]:
		# print site info
		proj.printSiteInfo(s)

		# get measurement directories for site s
		timeMeas = proj.getSiteTimeFiles(s)

		# loop over measurement folders and calculate spectra for each one
		for meas in timeMeas:
			# get measurement sample frequency
			fs = proj.getMeasSampleFreq(s, meas)			
			# check to see if in given frequency list
			if fs not in options["freqs"]:
				continue

			# print info
			generalPrint("Statistic Calculation", "Calculating stats for site {}, measurement {}".format(s, meas))	

			# define decimation parameters
			decParams = DecimationParams(fs)
			if len(options["evalfreq"]) == 0:
				decParams.setDecimationParams(options["declevels"], options["freqlevel"])
			else:
				decParams.setFrequencyParams(options["evalfreq"], options["declevels"], options["freqlevel"])
			#decParams.printInfo()
			numLevels = decParams.getNumLevels()

			# create the spectrum reader
			specReader = SpectrumReader(proj.getSpecDataPathMeas(s, meas))
			specReader.printInfo()

			# loop through decimation levels
			for iDec in xrange(0, numLevels):	
				# open the spectra file for the current decimation level
				check = specReader.openBinaryForReading('spectra', iDec)
				if not check:
					continue # probably because this decimation level not calculated
				specReader.printInfo()

				# get the channels
				dataChans = specReader.getChannels()
				if len(options["chans"]) > 0:
					dataChans = options["chans"]

				# get windows
				numWindows = specReader.getNumWindows()
				evalFreq = decParams.getEvalFrequenciesForLevel(iDec)
				fsDec = specReader.getSampleFreq()
				globalOffset = specReader.getGlobalOffset()
				fArray = specReader.getFrequencyArray()

				statHandlers = {}
				# create the statistic handlers
				for stat in options["stats"]:
					statName, statElements = statNameMap(stat)
					statHandlers[stat] = StatisticFrequency(statName)
					statHandlers[stat].setStatParams(numWindows, statElements, evalFreq)
					# statHandlers[stat].printInfo()										

				# loop over windows and calculate the relevant statistics
				for iW in xrange(0, numWindows):
					# get data
					winSpec = specReader.readBinaryWindowLocal(iW)
					globalIndex = iW + specReader.getGlobalOffset()
					# give the statistic calculator the spectra
					statCalculator.setSpectra(fArray, winSpec, evalFreq)
					# get the desired statistics
					for sH in statHandlers:
						data = statCalculator.getDataForStatName(sH)
						statHandlers[sH].addStat(iW, globalIndex, data)

				# save statistic
				for sH in statHandlers:
					statHandlers[sH].writeStatFile(proj.getStatDataPathMeas(s, meas), iDec)


# calculate statistics which use a remote reference
# for example, remote reference transfer functions
def projectRemoteStatCalc(proj, remoteSite, **kwargs):
	# default options
	options = parseKeywords(getDefaultOptions(proj), kwargs)			
	# get the reference time
	datetimeRef = proj.getRefTime()	
	# create the statistic calculator
	statCalculator = StatisticCalculator()	

	# loop over sites
	for s in options["sites"]:
		# print site info
		proj.printSiteInfo(s)

		# get measurement directories for site s
		rawMeas = proj.getSiteTimeFiles(s)

		# loop over measurement folders and calculate spectra for each one
		for meas in rawMeas:
			# get measurement sample frequency
			fs = proj.getMeasSampleFreq(s, meas)			
			# check to see if in given frequency list
			if fs not in options["freqs"]:
				continue	

			# print info
			generalPrint("Statistic Remote Calculation", "Calculating stats for site {}, measurement {} with reference {}".format(s, meas, remoteSite))

			# define decimation parameters
			decParams = DecimationParams(fs)
			if len(options["evalfreq"]) == 0:
				decParams.setDecimationParams(options["declevels"], options["freqlevel"])
			else:
				decParams.setFrequencyParams(options["evalfreq"], options["declevels"], options["freqlevel"])
			#decParams.printInfo()
			numLevels = decParams.getNumLevels()

			# create the window parameters for the win selector
			winParams = WindowParams(decParams)					
			# create the window selector and find the shared windows
			winSelector = WindowSelector(proj, fs, decParams, winParams)
			winSelector.setSites([s, remoteSite])
			winSelector.calcSharedWindows()	# calc shared windows between site and remote	

			# create the spectrum reader for the local site
			specReader = SpectrumReader(proj.getSpecDataPathMeas(s, meas))
			specReader.printInfo()	

			# loop through decimation levels
			for iDec in xrange(0, numLevels):	
				# open the spectra file for the current decimation level
				check = specReader.openBinaryForReading('spectra', iDec)
				# check = check and 
				if not check:
					continue # probably because this decimation level not calculated
				specReader.printInfo()

				# get the channels
				dataChans = specReader.getChannels()
				if len(options["chans"]) > 0:
					dataChans = options["chans"]

				# get a set of the shared windows at this decimation level
				# these are the global indices
				sharedWindows = winSelector.getSharedWindowsLevel(iDec)
				# get other information regarding only this spectra file
				numWindows = specReader.getNumWindows()				
				evalFreq = decParams.getEvalFrequenciesForLevel(iDec)
				fsDec = specReader.getSampleFreq()
				globalOffset = specReader.getGlobalOffset()
				fArray = specReader.getFrequencyArray()
				# now want to find the size of the intersection between the windows in this file and the shared windows
				sharedWindowsMeas = sharedWindows.intersection(set(np.arange(globalOffset, globalOffset+numWindows)))
				sharedWindowsMeas = sorted(list(sharedWindowsMeas))
				numSharedWindows = len(sharedWindowsMeas)

				statHandlers = {}
				# create the statistic handlers
				for stat in options["remotestats"]:
					statName, statElements = statNameMap(stat)
					statHandlers[stat] = StatisticFrequency(statName)
					# remember, this is with the remote reference, so the number of windows is number of shared windows
					statHandlers[stat].setStatParams(numSharedWindows, statElements, evalFreq)										

				# loop over windows and calculate the relevant statistics
				# loop over the shared windows between the remote station and local station
				for iW, globalWindow in enumerate(sharedWindowsMeas):
					# get data and set in the statCalculator
					winSpec = specReader.readBinaryWindowGlobal(globalWindow)
					statCalculator.setSpectra(fArray, winSpec, evalFreq)					
					# for the remote site, use the reader in win selector
					remoteSF, remoteReader = winSelector.getSpecReaderForWindow(remoteSite, iDec, globalWindow)
					winSpecRR = remoteReader.readBinaryWindowGlobal(globalWindow)
					statCalculator.addRemoteSpec(winSpecRR)	
					
					for sH in statHandlers:
						data = statCalculator.getDataForStatName(sH)
						statHandlers[sH].addStat(iW, globalWindow, data)					

				# save statistic
				for sH in statHandlers:
					statHandlers[sH].writeStatFile(proj.getStatDataPathMeas(s, meas), iDec)	


def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()	
	default["chans"] = []
	default["evalfreq"] = []	
	default["declevels"] = 7
	default["freqlevel"] = 6
	default["stats"] = ["absvalEqn", "coherence", "psd", "poldir", "transFunc", "resPhase", "partialcoh"]	
	default["remotestats"] = ["coherenceRR", "coherenceRREqn", "absvalRREqn", "transFuncRR", "resPhaseRR"]
	return default

def parseKeywords(default, keywords):
	# check keyword arguments
	for w in default:
		if w in keywords:
			default[w] = keywords[w]		
	return default		