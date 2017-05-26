#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('..', 'core'))
sys.path.append(os.path.join('..', 'utils'))
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
# import my classes
from project import Project
from calibrator import Calibrator
from dataWriterInternal import DataWriterInternal
# utils
from utilsProcess import *
from utilsIO import *


# resample data in the project
def projectInterpResamp(proj, **kwargs):
	generalPrint("Project Interp Resamp", "Resampling / Interpolating project time files with options: {}".format(kwargs))
	# resample info is a dictionary
	# {sampleRateToResample: sampleRateToResampleTo}
	# this function then does this for all measurement directories of that sample rate in the project
	options = parseKeywords(getDefaultOptions(proj), kwargs)	

	# create a data calibrator and writer instance
	cal = Calibrator(proj.getCalDataPath())
	if options["calibrate"]:
		cal.printInfo()		
	writer = DataWriterInternal()	

	# loop over sites
	for s in options["sites"]:
		# print site info
		proj.printSiteInfo(s)		

		for fs in options["freqs"]:
			timeFiles = proj.getSiteTimeFilesFs(s, fs)
			
			if len(timeFiles) == 0:
				continue # nothing to process

			# otherwise, resample
			for tF in timeFiles:
				# get the reader
				reader = proj.getMeasDataReader(s, tF)
				reader.printInfo()
				# get the dates for the data
				timeStart = reader.getStartDatetime()
				timeStop = reader.getStopDatetime()
				
				# now check the user provided dates
				# don't change timeStart, timeEnd yet because that breaks the checking
				if options["start"]:
					startUser = datetime.strptime(options["start"], "%Y-%m-%d %H:%M:%S")
					if startUser > timeStop: # this data has nothing to contribute in the optional date range
						continue
				if options["stop"]:
					stopUser = datetime.strptime(options["stop"], "%Y-%m-%d %H:%M:%S")
					if stopUser < timeStart: # this data has nothing to contribute in the optional date range
						continue

				# if the data contributes, copy in the data if relevant
				if options["start"]:
					timeStart = datetime.strptime(options["start"], "%Y-%m-%d %H:%M:%S")
				if options["stop"]:
					timeStop = datetime.strptime(options["stop"], "%Y-%m-%d %H:%M:%S")				
				
				# calculate the samples
				sampleStart, sampleEnd = reader.time2sample(timeStart, timeStop)
				dataStartTime, dataStopTime = reader.sample2time(sampleStart, sampleEnd) # need to do it this way round to protect against fractional sampling
				# now get the data
				data = reader.getPhysicalSamples(startSample=sampleStart, endSample=sampleEnd)		
				numSamples = reader.getNumSamples()	
				chans = reader.getChannels()				
				headers = reader.getHeaders()
				chanHeaders, chanMap = reader.getChanHeaders()	
				dataFs = fs						

				# if calibration is on, calibrate the data
				if options["calibrate"]:
					generalPrint("Project Interp Resamp", "Calibrating time data: site {}, time data {}".format(s, tF))
					sensors = reader.getSensors(reader.getChannels())
					serials = reader.getSerials(reader.getChannels())
					choppers = reader.getChoppers(reader.getChannels())	
					data = cal.calibrate(data, dataFs, sensors, serials, choppers)				

				# Interpolation to the second - do this first
				# make sure all the data starts on a full second
				if options["interp"]:
					if dataStartTime.microsecond != 0:
						generalPrint("Project Interp Resamp", "Interpolating to second: site {}, time data {} at {} Hz".format(s, tF, dataFs))
						# the recording is not on the second - NOTE, this will fail with longer sample periods (i.e. greater than a second)
						startTimeInterp, numSamplesInterp, dataInterp = interpolateToSecond(dataFs, dataStartTime, data)
						numSamples = numSamplesInterp
						dataStartTime = startTimeInterp
						data = dataInterp

				# Check if fs belongs to options["resamp"] keys and resample if does
				# resampling does not change the start time - hence dataStartTime is unchanged
				if dataFs in options["resamp"]: # then need to resample this data
					generalPrint("Project Interp Resamp", "Resampling site = {}, time data {} at {} Hz to {} Hz".format(s, tF, fs, options["resamp"][fs]))
					data = resample(data, dataFs, options["resamp"][dataFs])
					# update info for saving file
					dataFs = options["resamp"][fs]
					numSamples = data[chans[0]].size

				# recall, start times stay the same with resampling - only the sample rate changes and the number of samples
				# write out data - the data writer automatically deals with the end date
				outPath = os.path.join(proj.getTimeDataPathSite(s), options["prepend"] + tF + options["postpend"])
				writer.setOutPath(outPath)
				writer.writeData(headers, chanHeaders, data, 
					start_time=dataStartTime.strftime("%H:%M:%S.%f"), 
					start_date=dataStartTime.strftime("%Y-%m-%d"), 
					numSamples=numSamples,
					sample_freq=dataFs, 
					lsb_applied=True
				)
				writer.printInfo()				

def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()
	default["start"] = False
	default["stop"] = False
	default["calibrate"] = False	
	default["resamp"] = {}
	default["interp"] = False	
	default["prepend"] = "resamp_"
	default["postpend"] = ""
	return default

def parseKeywords(default, keywords):
	# check user options
	for w in default:
		if w in keywords:
			default[w] = keywords[w]	
	return default		