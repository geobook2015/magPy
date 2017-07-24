# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
"""
import os
import glob
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import numpy as np
# import dataReader
from dataReader import DataReader
# utils
from utilsIO import *
from utilsChecks import *
from utilsProcess import removeZeros, removeZerosSingle

class DataReaderSPAM(DataReader):
	"""
	One major difference with SPAM data is that you can have multiple data files in the same measurement directory
	Therefore, have to read all the XTR/XTRX files
	Read in the start and end times
	NOTE: For spam data, the start and end times in the XTR files
    are the time of the first sample and the time of the last sample
	"""

	###################
	### PREPARE
	##################
	def prepare(self):
		# get a list of the header and data files in the folder
		self.headerF = glob.glob(os.path.join(self.dataPath,"*.XTR"))
		if len(self.headerF) == 0:
			self.headerF = glob.glob(os.path.join(self.dataPath,"*.XTRX"))
		self.dataF = glob.glob(os.path.join(self.dataPath,"*.RAW"))
		# data byte information
		self.dataByteOffset = {} # because this potentially might be different for each file
		self.recChannels = {}
		self.dataByteSize = 4
		# data type
		self.dtype = np.float32
		# get the number of data files and header files
		# this should be equal
		self.numHeaderFiles = len(self.headerF)
		self.numDataFiles = len(self.dataF)

	###################
	### COUNT DATA
	##################
	# SPAM data is already scaled by LSB and in single precision float in volts
	# instead of going back to counts, which seems a bit pointless (as usually, we are only interested in the physical data)
	def getUnscaledSamples(self, **kwargs):
		# initialise chans, startSample and endSample with the whole dataset
		options = self.parseGetDataKeywords(kwargs)

		# get the files to read and the samples to take from them, in the correct order
		dataFilesToRead, samplesToRead, scalings = self.getDataFilesForSamples(options["startSample"], options["endSample"])
		# set up the dictionary to hold the data
		data = {}
		for chan in options["chans"]:
			data[chan] = np.zeros(shape=(options["endSample"] - options["startSample"] + 1), dtype=self.dtype)

		# loop through chans and get data
		sampleCounter = 0
		for dFile, sToRead, scalar in zip(dataFilesToRead, samplesToRead, scalings):
			# get samples - this is inclusive
			dSamples = sToRead[1] - sToRead[0] + 1
			dSamplesRead = dSamples*self.recChannels[dFile] # because spam files always record 5 channels
			# read the data
			byteOff = self.dataByteOffset[dFile] + sToRead[0]*self.recChannels[dFile]*self.dataByteSize
			dFilePath = os.path.join(self.getDataPath(), dFile)
			dataRead = np.memmap(dFilePath, dtype=self.dtype, mode="r", offset=byteOff, shape=(dSamplesRead))
			# now need to unpack this
			for chan in options["chans"]:
				# check to make sure channel exists
				self.checkChan(chan)
				# get the channel index - the chanIndex should give the right order in the data file
				# as it is the same order as in the header file
				chanIndex = self.chanMap[chan]
				# use the range sampleCounter -> sampleCounter +  dSamples, because this actually means sampleCounter + dSamples - 1
				# scale by the lsb scalar here - note that these can be different for each file in the run
				data[chan][sampleCounter : sampleCounter + dSamples] = dataRead[chanIndex:dSamplesRead:self.recChannels[dFile]]*scalar[chan]
				# previous implementation had no scaling
				# data[chan][sampleCounter : sampleCounter + dSamples] = dataRead[chanIndex:dSamplesRead:self.recChannels[dFile]]
			# increment sample counter
			sampleCounter = sampleCounter + dSamples # get ready for the next data read

		# return the data
		return data

	def getDataFilesForSamples(self, startSample, endSample):
		# have the datafiles saved in sample order beginning with the earliest first
		# go through each datafile and find the range to be read
		dataFilesToRead = []
		samplesToRead = []
		scalings = []
		for idx, dFile in enumerate(self.dataFileList):
			fileStartSamp = self.dataRanges[idx][0]
			fileEndSamp = self.dataRanges[idx][1]
			if fileStartSamp > endSample or fileEndSamp < startSample:
				continue # nothing to read from this file
			# in this case, there is some overlap with the samples to read
			dataFilesToRead.append(dFile)
			readFrom = 0 # i.e. the first sample in the datafile
			readTo = fileEndSamp - fileStartSamp # this the last sample in the file
			if fileStartSamp < startSample:
				readFrom = startSample - fileStartSamp
			if fileEndSamp > endSample:
				readTo = endSample - fileStartSamp
			# this is an inclusive number readFrom to readTo including readTo
			samplesToRead.append([readFrom, readTo])
			scalings.append(self.scalings[idx])
		return dataFilesToRead, samplesToRead, scalings

	###################
	### PHYSICAL DATA
	##################
	# get data in physical units
	# getUnscaled does actually scale to remove the gain - as this might be different for each individual data file
	# get physicalSamples puts everything in field units
	# using field units
	# electric field in mV/km
	# magnetic fields is mV - need to calibrate to get magnetic field in nT
	def getPhysicalSamples(self, **kwargs):
		# initialise chans, startSample and endSample with the whole dataset
		options = self.parseGetDataKeywords(kwargs)
		# get data
		data = self.getUnscaledSamples(chans=options["chans"], startSample=options["startSample"], endSample=options["endSample"])
		# the LSB is applied in getUnscaledSamples - this is for ease of calculation and because each data file in the run might have a separate lsb
		# so all that is left is to divide by the dipole length in km and remove the average
		for chan in options["chans"]:
			if chan == 'Ex':
				# multiply by 1000/self.getChanDx same as dividing by dist in km
				data[chan] = 1000*data[chan]/self.getChanDx(chan)
			if chan == 'Ey':
				# multiply by 1000/self.getChanDy same as dividing by dist in km
				data[chan] = 1000*data[chan]/self.getChanDy(chan)

			# if remove zeros - False by default
			if options["remzeros"]:
				data[chan] = removeZerosSingle(data[chan])
			# if remove nans - False by default
			if options["remnans"]:
				data[chan] = removeNansSingle(data[chan])
			# remove the average from the data - True by default
			if options["remaverage"]:
				data[chan] = data[chan] - np.average(data[chan])

		# if magnetic channel, just return
		return data

	###################
	### READ AND MERGE HEADERS
	###################
	def spamHeaders(self):
		sections = ['STATUS', 'TITLE', 'PROJECT', 'FILE', 'SITE', 'CHANNAME', 'DATA']
		sectionHeaders = {}
		sectionHeaders['STATUS'] = ['STATUS']
		sectionHeaders['TITLE'] = ['AUTHOR', 'VERSION', 'DATE', 'COMMENT']
		sectionHeaders['FILE'] = ['NAME', 'FREQBAND', 'DATE']
		sectionHeaders['CHANNAME'] = ['ITEMS', 'NAME']
		sectionHeaders['DATA'] = ['ITEMS', 'CHAN']
		return sections, sectionHeaders

	# some default chan headers
	def chanDefaults(self):
		chanH = {}
		chanH["gain_stage1"] = 1
		chanH["gain_stage2"] = 1
		chanH["hchopper"] = 0 # this depends on sample frequency
		chanH["echopper"] = 0
		# channel output information (sensor_type, channel_type, ts_lsb, pos_x1, pos_x2, pos_y1, pos_y2, pos_z1, pos_z2, sensor_sernum)
		chanH["ats_data_file"] = ""
		chanH["num_samples"] = 0
		chanH["sensor_type"] = ""
		chanH["channel_type"] = ""
		chanH["ts_lsb"] = 1
		# the lsb/scaling is not applied. data is raw voltage which needs to be scaled
		# an lsb is constructed from the scaling in the XTR/XTRX file to take the data to mV
		chanH["lsb_applied"] = False # check this
		chanH["pos_x1"] = 0
		chanH["pos_x2"] = 0
		chanH["pos_y1"] = 0
		chanH["pos_y2"] = 0
		chanH["pos_z1"] = 0
		chanH["pos_z2"] = 0
		chanH["sensor_sernum"] = 0
		return chanH

	# read in header file
	# there may be more than a single file in the folder
	# in that case - want to merge header files
	# they should all be the same sampling frequency
	def readHeader(self):
		# read header files
		self.headersList = []
		self.chanHeadersList = []
		for headerFile in self.headerF:
			if "xtrx" in headerFile.lower():
				headers, chanHeaders = self.readHeaderXTRX(headerFile)
			else:
				headers, chanHeaders = self.readHeaderXTR(headerFile)
			self.headersList.append(headers)
			self.chanHeadersList.append(chanHeaders)

		# check to make sure no gaps
		# calculate out the sample ranges
		# and list the data files for each sample
		self.mergeHeaders(self.headersList, self.chanHeadersList)

	# read xtr format
	def readHeaderXTR(self, headerFile):
		with open(headerFile, "r") as f:
			lines = f.readlines()
		sections, sectionHeaders = self.spamHeaders()
		sectionLines = {}
		# let's get data
		for line in lines:
			line = line.strip()
			line = line.replace("'", " ")
			# continue if line is empty
			if line == "":
				continue
			if "[" in line:
				sec = line[1:-1]
				sectionLines[sec] = []
			else:
				sectionLines[sec].append(line)
		# the base class is built around a set of headers based on ATS headers
		# though this is a bit more work here, it saves lots of code repetition
		headers = {}
		# recording information (start_time, start_date, stop_time, stop_date, ats_data_file)
		fileLine = sectionLines["FILE"][0]
		fileSplit = fileLine.split()
		headers["sample_freq"]	= np.absolute(float(fileSplit[-1]))
		timeLine = sectionLines["FILE"][2]
		timeSplit = timeLine.split()
		# these are the unix time stamps
		startDate = float(timeSplit[1] + "." + timeSplit[2])
		datetimeStart = datetime.utcfromtimestamp(startDate)
		stopDate = float(timeSplit[3] + "." + timeSplit[4])
		datetimeStop = datetime.utcfromtimestamp(stopDate)
		headers["start_date"] = datetimeStart.strftime("%Y-%m-%d")
		headers["start_time"] = datetimeStart.strftime("%H:%M:%S.%f")
		headers["stop_date"] = datetimeStop.strftime("%Y-%m-%d")
		headers["stop_time"] = datetimeStop.strftime("%H:%M:%S.%f")
		# here calculate number of samples
		deltaSeconds = (datetimeStop - datetimeStart).total_seconds()
		# calculate number of samples - have to add one because the time given in SPAM recording is the actual time of the last sample
		numSamples = int(deltaSeconds*headers["sample_freq"]) + 1
		# put these in headers for ease of future calculations in merge headers
		headers["num_samples"] = numSamples
		headers["ats_data_file"] = fileSplit[1] # spam datasets only have the one data file for all channels
		# data information (meas_channels, sample_freq)
		chanLine = sectionLines["CHANNAME"][0]
		headers["meas_channels"] = chanLine.split()[1] # this gets reformatted to an int later
		numChansInt = int(headers["meas_channels"])
		# deal with the channel headers
		chanHeaders = []
		for iChan in xrange(0, numChansInt):
			chanH = self.chanDefaults()
			# set the sample frequency from the main headers
			chanH["sample_freq"] = headers["sample_freq"]
			# line data - read through the data in the correct channel order
			chanLine = sectionLines["CHANNAME"][iChan + 1]
			chanSplit = chanLine.split()
			dataLine = sectionLines["DATA"][iChan + 1]
			dataSplit = dataLine.split()
			# channel input information (gain_stage1, gain_stage2, hchopper, echopper)
			chanH["gain_stage1"] = 1
			chanH["gain_stage2"] = 1
			# channel output information (sensor_type, channel_type, ts_lsb, pos_x1, pos_x2, pos_y1, pos_y2, pos_z1, pos_z2, sensor_sernum)
			chanH["ats_data_file"] = fileSplit[1]
			chanH["num_samples"] = numSamples

			# channel information
			chanH["channel_type"] = consistentChans(chanSplit[2]) # spams often use Bx, By - use H within the software as a whole
			# the sensor number is a bit of a hack - want MFSXXe or something - add MFS in front of the sensor number
			# this is liable to break
			# at the same time, set the chopper
			calLine = sectionLines["200{}003".format(iChan + 1)][0]
			calSplit = calLine.split()
			if isMagnetic(chanH["channel_type"]):
				chanH["sensor_sernum"] = calSplit[2] # the last three digits is the serial number
				sensorType = calSplit[1].split("_")[1][-2:]
				chanH["sensor_type"] = "MFS{:02d}".format(int(sensorType))
				if "LF" in calSplit[1]:
					chanH["hchopper"] = 1
			else:
				chanH["sensor_type"] = "ELC00"
				if "LF" in calLine:
					chanH["echopper"] = 1

			# the scaling - recall, the data is raw voltage of sensors
			# gain needs to be removed - this is in the scaling = 1/1000*total_gain (gain1*gain2)
			# what we want here is 1000/total_gain = 1000000*scaling
			# the E fields also need their polarity reversed (from email with Reinhard)
			# NOTE: The scaling can be different for each dataset in the directory
			scaling = float(dataSplit[-2])
			if isElectric(chanH["channel_type"]):
				# the factor of 100000 is not entirely clear
				lsb = 1000000.0*scaling
				lsb = -1*lsb # reverse polarity
			else:
				# the spam scaling in the xtr file for magnetic fields includes the static gain correction
				# however, a static gain correction is applied in the calibration
				# in order to avoid duplication, the scaling in the xtr file ignored for the magnetic channels
				lsb = -1000.0 # volts to millivolts and a minus to switch polarity
			chanH["ts_lsb"] = lsb

			# the distances
			if chanSplit[2] == "Ex":
				chanH["pos_x1"] = float(dataSplit[4])/2
				chanH["pos_x2"] = chanH["pos_x1"]
			if chanSplit[2] == "Ey":
				chanH["pos_y1"] = float(dataSplit[4])/2
				chanH["pos_y2"] = chanH["pos_y1"]
			if chanSplit[2] == "Ez":
				chanH["pos_z1"] = float(dataSplit[4])/2
				chanH["pos_z2"] = chanH["pos_z1"]

			# append chanHeaders to the list
			chanHeaders.append(chanH)

		# check information from raw file headers
		self.headersFromRawFile(headers["ats_data_file"], headers)

		# return the headers and chanHeaders from this file
		return headers, chanHeaders

	# read the newer xtrx format
	def readHeaderXTRX(self, headerFile):
		self.printWarning("Reading of XTRX files has not been implemented yet")
		headers = {}
		chanHeaders = []
		return headers, chanHeaders

	# read the headers from the raw file
	# and figure out the data byte offset
	def headersFromRawFile(self, rawFile, headers):
		# in particular, are interested in number of samples
		# and the size of the header
		dFile = open(os.path.join(self.dataPath, rawFile), "r") # the headers are in ascii
		generalHeaderString = dFile.read(1000) # this should be long enough
		generalSplit = generalHeaderString.split()
		# read GENERAL HEADER
		generalHeader = {}
		generalHeader["recLength"] = int(generalSplit[0])
		generalHeader["fileType"] = generalSplit[1]
		generalHeader["wordLength"] = int(generalSplit[2])
		generalHeader["version"] = generalSplit[3]
		generalHeader["procId"] = generalSplit[4]
		generalHeader["numCh"] = int(generalSplit[5])
		generalHeader["totalRec"] = int(generalSplit[6])
		generalHeader["firstEvent"] = int(generalSplit[7])
		generalHeader["numEvent"] = int(generalSplit[8])
		generalHeader["extend"] = int(generalSplit[9])
		# NOTE: extended header ignored

		# read EVENT HEADER - there can be multiple of these, but normally only the one
		# events are largely deprecated. Only a single event is used
		eventHeaders = []
		fileSize = os.path.getsize(os.path.join(self.getDataPath(), rawFile))
		record = generalHeader["firstEvent"]
		for ir in xrange(0, generalHeader["numEvent"]):
			seekPt = (record-1)*generalHeader["recLength"]
			if not seekPt > fileSize:
				dFile.seek(seekPt, 0) # seek from beginning of file
    			eventString = dFile.read(1000) # read extra to make sure
    			eventSplit = eventString.split()
    			eH = {}
		        eH["start"] = int(eventSplit[0])
		        eH["startms"] = int(eventSplit[1])
		        eH["stop"] = int(eventSplit[2])
		        eH["stopms"] = int(eventSplit[3])
		        eH["cvalue1"] = float(eventSplit[4])
		        eH["cvalue2"] = float(eventSplit[5])
		        eH["cvalue3"] = float(eventSplit[6])
		        eH["EHInfile"] = int(eventSplit[7])
		        eH["nextEH"] = int(eventSplit[8])
		        eH["previousEH"] = int(eventSplit[9])
		        eH["numData"] = int(eventSplit[10])
		        eH["startData"] = int(eventSplit[11])
		        eH["extended"] = int(eventSplit[12])
		        eventHeaders.append(eH)
		        if eH["nextEH"] < generalHeader["totalRec"]:
		        	record = eH["nextEH"] # set to go to next eH
		        else:
		            break # otherwise break out of for loops
		# close the data file
		dFile.close()
		# now compare number of samples with that calculated previously
		if eventHeaders[0]["numData"] != headers["num_samples"]:
			self.printWarning("Data file: {}".format(dFile))
			self.printWarning("Number of samples in raw file header {} does not equal that calculated from data {}".format(eventHeaders[0]["numData"], headers["num_samples"]))
			self.printWarning("Number of samples calculated from data will be used")
		# set the byte offset for the file
		self.dataByteOffset[rawFile] = (eventHeaders[0]["startData"] - 1)*generalHeader["recLength"]
		self.recChannels[rawFile] = generalHeader["numCh"]

	# merge headers
	# this will check all the header files to see if there are any gaps
	# get the sample ranges for each file
	# calculate total number of samples
	# set the start and end time of the recording
	# and set class variables datetimeStart and datetimeStop
	def mergeHeaders(self, headersList, chanHeadersList):
		# take the first header as an example
		self.headers = headersList[0]
		self.chanHeaders = chanHeadersList[0]
		if len(headersList) == 1:
			# just fill in the data file list and data ranges
			self.dataFileList = [self.headers["ats_data_file"]]
			self.dataRanges = [[0, self.headers["num_samples"]-1]]
			self.scalings = []
			tmp = {}
			for cHeader in self.chanHeaders:
				tmp[cHeader["channel_type"]] = cHeader["ts_lsb"]
			self.scalings.append(tmp)
			return # then there was only one file - no need to do all the below

		# make sure that all headers have the same sample rate
		# and save the start and stop times and dates
		startTimes = []
		stopTimes = []
		numSamples = []
		for idx, header in enumerate(headersList):
			if header["sample_freq"] != self.headers["sample_freq"]:
				self.printWarning("Not all datasets in {} have the same sample frequency".format(self.dataPath))
				self.printWarning("Exiting")
				exit()
			if header["meas_channels"] != self.headers["meas_channels"]:
				self.printWarning("Not all datasets in {} have the same number of channels".format(self.dataPath))
				self.printWarning("Exiting")
				exit()
			# now store startTimes, stopTimes and numSamples
			# do this as datetimes, will be easier
			startString = "{} {}".format(header["start_date"], header["start_time"])
			stopString = "{} {}".format(header["stop_date"], header["stop_time"])
			datetimeStart = datetime.strptime(startString, "%Y-%m-%d %H:%M:%S.%f")
			datetimeStop = datetime.strptime(stopString, "%Y-%m-%d %H:%M:%S.%f")
			startTimes.append(datetimeStart)
			stopTimes.append(datetimeStop)
			numSamples.append(header["num_samples"])
		# check the start and end times
		sampleTime = timedelta(seconds=1.0/self.headers["sample_freq"])
		# sort by start times
		sortIndices = sorted(range(len(startTimes)), key=lambda k: startTimes[k])
		# now sort stop times by the same indices
		check = True
		for i in xrange(1, self.numHeaderFiles):
			# get the stop time of the previous dataset
			stopTimePrev = stopTimes[sortIndices[i-1]]
			startTimeNow = startTimes[sortIndices[i]]
			if startTimeNow != stopTimePrev + sampleTime:
				self.printWarning("There is a gap between the datafiles in {}".format(self.dataPath))
				self.printWarning("Please separate out datasets with gaps into separate folders")
				# print out where the gap was found
				self.printWarning("Gap found between datafiles:")
				self.printWarning("1. {}".format(headersList[sortIndices[i-1]]["ats_data_file"]))
				self.printWarning("2. {}".format(headersList[sortIndices[i]]["ats_data_file"]))
				# set check as false
				check = False
		# if did not pass check, then exit
		if not check:
			exit()

		# make sure there are no gaps
		totalSamples = sum(numSamples)

		# get a list of all the datafiles, scalings and the sample ranges
		self.dataFileList = []
		self.dataRanges = []
		self.scalings = []
		sample = -1
		# now need some sort of lookup table to say where the sample ranges are
		for i in xrange(0, self.numHeaderFiles):
			iSort = sortIndices[i] # get the sorted index
			self.dataFileList.append(headersList[iSort]["ats_data_file"])
			startSample = sample + 1
			endSample = startSample + numSamples[iSort] - 1 # -1 because this is inclusive of the start sample
			self.dataRanges.append([startSample, endSample])
			# increment sample
			sample = endSample
			# save the scalings for each chan
			tmp = {}
			for cHeader in self.chanHeadersList[iSort]:
				tmp[cHeader["channel_type"]] = cHeader["ts_lsb"]
			self.scalings.append(tmp)

		# now set the LSB information for the chanHeaders
		# i.e. if they change, this should reflect that
		for i in xrange(0, len(self.chanHeaders)):
			chan = self.chanHeaders[i]["channel_type"]
			lsbSet = set()
			for scalar in self.scalings:
				lsbSet.add(scalar[chan])
			if len(lsbSet) == 1:
				self.chanHeaders[i]["ts_lsb"] = list(lsbSet)[0]
			else:
				self.printWarning("Multiple different LSB values found for chan {}: {}".format(chan, list(lsbSet)))
				self.printWarning("This is handled, but the header information given will show only a single LSB value")
				self.chanHeaders[i]["ts_lsb"] = list(lsbSet)[0]

		# set start and end time for headers and chan headers
		# do the same with number of samples
		datetimeStart = min(startTimes)
		datetimeStop = max(stopTimes)
		self.headers["start_date"] = datetimeStart.strftime("%Y-%m-%d")
		self.headers["start_time"] = datetimeStart.strftime("%H:%M:%S.%f")
		self.headers["stop_date"] = datetimeStop.strftime("%Y-%m-%d")
		self.headers["stop_time"] = datetimeStop.strftime("%H:%M:%S.%f")
		self.headers["num_samples"] = totalSamples
		# set datafiles = the whole list of datafiles
		self.headers["ats_data_file"] = self.dataFileList
		for iChan in xrange(0, len(self.chanHeaders)):
			self.chanHeaders[iChan]["start_date"] = datetimeStart.strftime("%Y-%m-%d")
			self.chanHeaders[iChan]["start_time"] = datetimeStart.strftime("%H:%M:%S.%f")
			self.chanHeaders[iChan]["stop_date"] = datetimeStop.strftime("%Y-%m-%d")
			self.chanHeaders[iChan]["stop_time"] = datetimeStop.strftime("%H:%M:%S.%f")
			self.chanHeaders[iChan]["num_samples"] = totalSamples
			self.chanHeaders[iChan]["ats_data_file"] = self.dataFileList

	###################
	### DEBUG
	##################
	def printInfoBegin(self):
		self.printText("####################")
		self.printText("EMERALD READER INFO BEGIN")
		self.printText("####################")

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("EMERALD READER INFO BEGIN")
		self.printText("####################")

	def printDataFileList(self):
		self.printText("####################")
		self.printText("EMERALD READER DATA FILE LIST BEGIN")
		self.printText("####################")
		self.printText("Data File\t\tSample Ranges")
		for dFile, sRanges in zip(self.dataFileList, self.dataRanges):
			self.printText("{}\t\t{} - {}".format(dFile, sRanges[0], sRanges[1]))
		self.printText("Total samples = {}".format(self.getNumSamples()))
		self.printText("####################")
		self.printText("EMERALD READER DATA FILE LIST END")
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Emerald Reader Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Emerald Reader Warning", warnStr)
