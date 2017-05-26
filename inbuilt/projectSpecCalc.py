#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('..', 'core'))
sys.path.append(os.path.join('..', 'utils'))
import numpy as np
from datetime import datetime
# import my classes
from project import Project
#from dataReaderATS import DataReaderATS
#from dataReaderSpam import DataReaderSPAM
#from dataReaderInternal import DataReaderInternal
from calibrator import Calibrator
from decimationParameters import DecimationParams
from windowParameters import WindowParams
from decimator import Decimator
from windower import Windower
from spectrumCalculator import SpectrumCalculator
from spectrumWriter import SpectrumWriter
# utilities
from utilsPlotter import *
from utilsProcess import *


# calculate spectra for the project
# the philosophy here is that spectra is calculated out for all data
# and then later limited using statistics and time constraints
def projectSpecCalc(proj, **kwargs):
	generalPrint("ProjectSpecCalc", "Calculating spectra for project with options: {}".format(kwargs))
	# default options
	options = parseKeywords(getDefaultOptions(proj), kwargs)

	# calculate the spectra
	# get the reference time
	datetimeRef = proj.getRefTime()	
	# prepare the calibrator
	cal = Calibrator(proj.getCalDataPath())
	if options["calibrate"]:
		cal.printInfo()			

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
			if int(fs) not in options["freqs"]:
				continue
			
			# print measurement info
			proj.printMeasInfo(s, meas)

			# get measurement start and end times
			datetimeStart = proj.getMeasStartTime(s, meas)
			datetimeEnd = proj.getMeasEndTime(s, meas)			

			# get data, sensor info, serial info, chopper info for calibration
			reader = proj.getMeasDataReader(s, meas)
			# get data start and end times - these may not be equal to startDate and endDate
			# there is the issue of ats data recording end time as one sample late
			# hence get the actual end time
			dataStartTime, dataEndTime = reader.getDataTimes(datetimeStart, datetimeEnd)				
			dataChans = reader.getChannels()
			if len(options["chans"]) > 0:
				dataChans = options["chans"]
			# alternatively, could simply do getPhysicalSamples() and get all data that way
			data = reader.getPhysicalData(dataStartTime, dataEndTime, chans=dataChans)
			
			if options["calibrate"]:
				# do the calibration here				
				sensors = reader.getSensors(dataChans)
				serials = reader.getSerials(dataChans)
				choppers = reader.getChoppers(dataChans)	
				data = cal.calibrate(data, fs, sensors, serials, choppers)

			# notch filter if required
			for n in options["notch"]:
				for c in data:
					data[c] = notchFilter(data[c], fs, n, n/5.0)				

			# define decimation parameters
			decParams = DecimationParams(fs)
			if len(options["evalfreq"]) == 0:
				decParams.setDecimationParams(options["declevels"], options["freqlevel"])
			else:
				decParams.setFrequencyParams(options["evalfreq"], options["declevels"], options["freqlevel"])
			decParams.printInfo()
			numLevels = decParams.getNumLevels()

			# now do window parameters
			winParams = WindowParams(decParams)
			# winParams.printInfo()

			# create the decimator
			dec = Decimator(data, fs, decParams)
			# dec.printInfo()

			# loop through decimation levels
			for iDec in xrange(0, numLevels):	
				# get the data for the current level
				check = dec.incrementLevel()
				if not check:
					break # not enough data
				#dec.printInfo()
				data = dec.getData()

				# create the windower and give it window parameters for current level
				fsDec = dec.getSampleFreq()
				win = Windower(datetimeRef, dataStartTime, data, fsDec, winParams.getWindowSize(iDec), winParams.getOverlap(iDec))
				numWindows = win.getNumWindows()	
				# win.printInfo()
				if numWindows < 2:
					break # do no more decimation	

				# create the spectrum calculator and statistics calculators
				specCalc = SpectrumCalculator(fsDec, winParams.getWindowSize(iDec))	
				
				# get ready a file to save the spectra
				specWrite = SpectrumWriter(proj.getSpecDataPathMeas(s, meas), datetimeRef)
				specWrite.openBinaryForWriting("spectra", iDec, fsDec, winParams.getWindowSize(iDec), 
					winParams.getOverlap(iDec), win.getGlobalWindowOffset(), numWindows, dataChans)

				# loop though windows, calculate spectra and save
				for iW in xrange(0, numWindows):
					# get the window data
					winData = win.getData(iW)			

					# calculate spectra
					f, specData = specCalc.calcFourierCoeff(winData)

					# write out spectra
					specWrite.writeBinary(specData, iW)
				
				# close spectra and stat files
				specWrite.closeFile()

def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()	
	default["chans"] = []
	default["evalfreq"] = []
	default["declevels"] = 7
	default["freqlevel"] = 6
	default["calibrate"] = True
	default["notch"] = []
	return default

def parseKeywords(default, keywords):
	# check user options
	for w in default:
		if w in keywords:
			default[w] = keywords[w]	
	return default

