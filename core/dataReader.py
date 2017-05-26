# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop

class attributes
self.chans = list of the channels i.e. 'Ex', 'Ey'
self.chanMap = mapping channel to index

NOTE: Different data formats give the end time as either the time of the last sample (SPAM)
Or the time of the sample after the final one (ATS)
For consistency purposes, the time of the last sample will be used internally
"""
import os
import glob
from datetime import datetime, timedelta
import numpy as np
from copy import deepcopy
# utils
from utilsIO import *
from utilsChecks import *

class DataReader(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, dataPath):
		self.dataPath = dataPath
		# get a list of the xml files in the folder
		self.prepare()
		if not self.checkFiles():
			self.printWarning("No header or data files found...exiting")
			exit()
		self.readHeader()
		self.formatHeaderData()
		self.initialise()

	# prepare sets the options for the individual readers
	def prepare(self):
		self.headerF = glob.glob(os.path.join(self.dataPath,"*.hdr"))
		self.dataF = glob.glob(os.path.join(self.dataPath,"*.npy"))
		self.dataByteOffset = 0
		self.dataByteSize = 4

	# check to make sure some files found
	def checkFiles(self):
		check = True
		if len(self.headerF) == 0:
			check = check and False
			self.printWarning("No header files found")
		if len(self.dataF) == 0:
			check = check and False
			self.printWarning("No data files found")
		return check

	###################
	### GET CLASS VARS
	##################
	def getDataPath(self):
		return self.dataPath

	###################
	### GET GLOBAL HEADERS
	### THESE SHOULD BE PUBLICALLY AVAILABLE
	##################
	# get the channels in the data
	def getChannels(self):
		return deepcopy(self.chans)

	# get all the headers
	def getHeaders(self):
		return deepcopy(self.headers)

	# get a header value
	def getHeaderVal(self, headerName):
		return self.headers[headerName]		

	# these could all be put in a parent class
	def getNumChannels(self):
		return self.getHeaderVal("meas_channels")

	def getSampleFreq(self):
		return self.getHeaderVal("sample_freq")

	def getSampleFreqInt(self):
		return int(self.getHeaderVal("sample_freq"))

	def getSampleRate(self):
		return 1.0/self.getSampleFreq()

	def getNumSamples(self):
		return self.getHeaderVal("num_samples")

	def getStartDate(self):
		return self.getHeaderVal("start_date")

	def getStartTime(self):
		return self.getHeaderVal("start_time")

	def getStartDatetime(self):
		return deepcopy(self.datetimeStart)

	def getStartDatetimeTuple(self):
		return self.datetimeStart.utctimetuple()

	def getStartDatetimeString(self):
		return self.datetimeStart.strftime("%Y-%m-%d %H:%M:%S.%f")		

	def getStopTime(self):
		return self.getHeaderVal("stop_time")

	def getStopDate(self):
		return self.getHeaderVal("stop_date")

	def getStopDatetime(self):
		return deepcopy(self.datetimeStop)

	def getStopDatetimeTuple(self):
		return self.datetimeStop.utctimetuple()		

	def getStopDatetimeString(self):
		return self.datetimeStop.strftime("%Y-%m-%d %H:%M:%S.%f")			

	def getGain1(self):
		gain1 = np.zeros(shape=(self.numChannels), dtype=bool)
		for iChan in xrange(0, self.numChannels):
			gain1[iChan] = self.getChanGain1(iChan)
		return gain1

	def getGain2(self):
		gain2 = np.zeros(shape=(self.numChannels), dtype=bool)
		for iChan in xrange(0, self.numChannels):
			gain2[iChan] = self.getChanGain2(iChan)
		return gain2

	###################
	### GET CHANNEL HEADERS
	### NO REAL NEED TO EXPOSE THESE
	##################
	def getChanHeaders(self, **kwargs):
		chans = self.getChannels()
		if "chans" in kwargs:
			chans = kwargs["chans"]
			chanHeaders = []
			chanMap = {}
			for idx, c in enumerate(chans):
				headerI = self.getChanMap[c]
				chanHeaders.append(self.getChanHeaders[headerI])
				chanMap[c] = idx
			return chanHeaders, chanMap
		else: 
			return deepcopy(self.chanHeaders), deepcopy(self.chanMap)

	def getChanMap(self):
		return self.chanMap

	def getChanHeader(self, chan, header):
		# should check to see if chan exists
		self.checkChan(chan)
		iChan = self.chanMap[chan]
		return self.chanHeaders[iChan][header]	

	def getChanType(self, chan):
		return self.getChanHeader(chan, "channel_type")						

	def getChanGain1(self, chan):		
		return self.getChanHeader(chan, "gain_stage1")

	def getChanGain2(self, chan):		
		return self.getChanHeader(chan, "gain_stage2")

	def getChanSamples(self, chan):		
		return self.getChanHeader(chan, "num_samples")
	
	def getChanSampleFreq(self, chan):
		return self.getChanHeader(chan, "sample_freq")

	def getChanLSB(self, chan):		
		return self.getChanHeader(chan, "ts_lsb")

	def getChanLSBApplied(self, chan):
		return self.getChanHeader(chan, "lsb_applied")

	def getChanDataFile(self, chan):	
		return self.getChanHeader(chan, "ats_data_file")

	# Purely for ease, just absolute these distance values
	# for Dx, Dy and Dz and then add them up
	def getChanDx(self, chan):		
		x1 = np.absolute(self.getChanHeader(chan, "pos_x1"))
		x2 = np.absolute(self.getChanHeader(chan, "pos_x2"))
		return x2 + x1

	def getChanDy(self, chan):		
		y1 = np.absolute(self.getChanHeader(chan, "pos_y1"))
		y2 = np.absolute(self.getChanHeader(chan, "pos_y2"))
		return y2 + y1

	def getChanDz(self, chan):	
		z1 = np.absolute(self.getChanHeader(chan, "pos_z1"))
		z2 = np.absolute(self.getChanHeader(chan, "pos_z2"))
		return z2 + z1

	def getChanSensor(self, chan):
		return self.getChanHeader(chan, "sensor_type")		

	def getChanSerial(self, chan):
		return self.getChanHeader(chan, "sensor_sernum")

	def getChanChopper(self, chan):
		echopper = self.getChanHeader(chan, "echopper")
		hchopper = self.getChanHeader(chan, "hchopper")
		# return true if the chopper amplifier was on
		if isElectric(chan) and echopper:
			return True
		if isMagnetic(chan) and hchopper:
			return True
		return False

	###################
	### GET INFO FOR MULTIPLE CHANNELS
	##################
	def getSensors(self, chans):
		sensors = {}
		for chan in chans:
			sensors[chan] = self.getChanSensor(chan)
		return sensors

	def getSerials(self, chans):
		serials = {}
		for chan in chans:
			serials[chan] = self.getChanSerial(chan)
		return serials

	def getChoppers(self, chans):
		choppers = {}
		for chan in chans:
			choppers[chan] = self.getChanChopper(chan)		
		return choppers				

	###################
	### GET INFO FOR MULTIPLE CHANNELS
	##################
	def setHeader(headerName, headerVal):
		self.headerName = headerVal

	def setChanHeader(chan, headerName, headerVal):
		chanIndex = self.chanMap[chan]
		self.chanHeaders[chanIndex][headerName] = headerVal

	###################
	### COUNT DATA
	##################
	# get data in raw count format
	# this function works for ATS and internal format
	# spam data should has its own implmentation
	def getUnscaledSamples(self, **kwargs):
		# initialise chans, startSample and endSample with the whole dataset
		options = self.parseGetDataKeywords(kwargs)
	
		# get samples - this is inclusive
		dSamples = options["endSample"] - options["startSample"] + 1

		# loop through chans and get data
		data = {}
		for chan in options["chans"]:		
			# check to make sure channel exists
			self.checkChan(chan)
			# get data file
			dFile = os.path.join(self.getDataPath(), self.getChanDataFile(chan))
			# get the data
			byteOff = self.dataByteOffset + options["startSample"] * self.dataByteSize
			# now check if lsb applied or not and read data as float32 or int32 accordingly
			if self.getChanLSBApplied(chan):
				data[chan] = np.memmap(dFile, dtype="float32", mode="r", offset=byteOff, shape=(dSamples))
			else:
				data[chan] = np.memmap(dFile, dtype="int32", mode="r", offset=byteOff, shape=(dSamples))
		# return data for the channels
		return data	
	
	def getUnscaledData(self, startTime, endTime, **kwargs):
		options = self.parseGetDataKeywords(kwargs)	
		startSample, endSample = self.time2sample(startTime, endTime)
		return self.getUnscaledSamples(chans=options["chans"], startSample=startSample, endSample=endSample)	

	###################
	### PHYSICAL DATA
	##################
	# get data in physical units
	# ATS provides unscaled data in counts (*lsb gives mV)
	# SPAM provides unscaled data in 
	# the count data/unscaled data * lsb is assumed to give data in mV
	# using field units
	# electric field in mV/km
	# magnetic fields is mV - need to calibrate to get magnetic field in nT
	def getPhysicalSamples(self, **kwargs):
		# initialise chans, startSample and endSample with the whole dataset
		options = self.parseGetDataKeywords(kwargs)		
		
		# get data
		data = self.getUnscaledSamples(chans=options["chans"], startSample=options["startSample"], endSample=options["endSample"])
		# multiply each chan by least significant bit of chan
		for chan in options["chans"]:	
			if not self.getChanLSBApplied(chan):
				# apply LSB
				data[chan] = data[chan]*self.getChanLSB(chan) # this gives data in mV
				# print chan, self.getChanLSB(chan)
				# divide by the distance - this should only be for the electric channels
				# again, this might already be applied
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
			# do this after all scaling and removing nans and zeros
			if options["remaverage"]:
				data[chan] = data[chan] - np.average(data[chan])					
		# if magnetic channel, just return
		return data

	def getPhysicalData(self, startTime, endTime, **kwargs):
		options = self.parseGetDataKeywords(kwargs)		
		startSample, endSample = self.time2sample(startTime, endTime)
		return self.getPhysicalSamples(chans=options["chans"], startSample=startSample, endSample=endSample)

	def parseGetDataKeywords(self, keywords):
		# set the defaults
		options = {}
		options["chans"] = self.getChannels()
		options["startSample"] = 0
		options["endSample"] = self.getNumSamples()-1
		options["startTime"] = self.getStartDatetime()
		options["endTime"] = self.getStopDatetime()
		options["remaverage"] = True
		options["remzeros"] = False
		options["remnans"] = False
		# now take the options from the keywords
		for w in options:
			if w in keywords:
				options[w] = keywords[w]
		# do some checks
		if options["endSample"] >= self.getNumSamples():
			options["endSample"] = self.getNumSamples()-1
			self.printWarning("End sample greater than number of samples. Adjusted to {:d}".format(endSample))	
		# return
		return options

	###################
	### TIME INFORMATION
	### THIS NEEDS A FEW MORE ERROR CHECKS
	##################
	# recall
	# the first sample is zero
	def time2sample(self, timeStart, timeEnd):
		# check to see times within range
		timeStart, timeEnd = self.getDataTimes(timeStart, timeEnd)
		# start sample
		deltaStart = timeStart - self.getStartDatetime()
		sampleStart = deltaStart.total_seconds()*self.getSampleFreq()	
		sampleStart = int(round(sampleStart)) # this will hopefully deal with fractional sampling
		# end sample	
		deltaEnd = timeEnd - timeStart
		deltaSamples = deltaEnd.total_seconds()/self.getSampleRate()
		deltaSamples = int(round(deltaSamples))
		sampleEnd = sampleStart + deltaSamples
		# return samples	
		return sampleStart, sampleEnd		

	def sample2time(self, sampleStart, sampleEnd):
		# convert samples to some data format
		deltaStart = timedelta(seconds=self.getSampleRate()*sampleStart)
		# delta end is inclusive
		deltaEnd = timedelta(seconds=self.getSampleRate()*(sampleEnd-sampleStart))
		timeStart = self.getStartDatetime() + deltaStart
		timeEnd = timeStart + deltaEnd
		return timeStart, timeEnd

	def getDataTimes(self, timeStart, timeEnd):
		deltaStart = timeStart - self.getStartDatetime()
		deltaEnd = self.getStopDatetime() - timeEnd
		if deltaStart.total_seconds() < 0:
			self.printText("Date {} before start of recording. Start date adjusted to {}".format(timeStart, self.getStartDatetime()))			
			timeStart = self.getStartDatetime() 	
		if deltaEnd.total_seconds() < 0:
			self.printText("Date {} after end of recording. Stop date adjusted to {}".format(timeEnd, self.getStopDatetime()))			
			timeEnd = self.getStopDatetime()
		return timeStart, timeEnd			

	###################
	### READ HEADER
	### each subclass should have a read header function
	### which reads the appropriate data format	
	##################
	# set the data types for each header
	def intHeaders(self):
		intGlobal = ["meas_channels"]
		intChan = ["gain_stage1", "gain_stage2", "hchopper", "echopper", "num_samples", "sensor_sernum"]
		return intGlobal, intChan

	def floatHeaders(self):
		floatGlobal = ["sample_freq"]
		floatChan = ["sample_freq", "ts_lsb", "pos_x1", "pos_x2", "pos_y1", "pos_y2", "pos_z1", "pos_z2"]
		return floatGlobal, floatChan

	def boolHeaders(self):
		boolGlobal = []
		boolChan = ["lsb_applied"]
		return boolGlobal, boolChan

	# read the header information
	def readHeader(self):
		# this is implemented in the child classes
		return
	
	# deal with the data formats
	def formatHeaderData(self):
		# do the int formatting
		intGlobal, intChan = self.intHeaders()
		floatGlobal, floatChan = self.floatHeaders()
		boolGlobal, boolChan = self.boolHeaders()
		# deal with the global headers
		for h in intGlobal:
			self.headers[h] = int(self.headers[h])
		for h in floatGlobal:
			self.headers[h] = float(self.headers[h])
		for h in boolGlobal:
			self.headers[h] = bool(self.headers[h])
		# deal with the channel headers
		numChans = len(self.chanHeaders)
		for iChan in xrange(0, numChans):
			for cH in self.chanHeaders[iChan]:
				if cH in intChan:
					self.chanHeaders[iChan][cH] = int(self.chanHeaders[iChan][cH])
				if cH in floatChan:
					self.chanHeaders[iChan][cH] = float(self.chanHeaders[iChan][cH])
				if cH in boolChan:
					if isinstance(self.chanHeaders[iChan][cH], basestring): # if string, parse true and false
						if self.chanHeaders[iChan][cH] == "True":
							self.chanHeaders[iChan][cH] = True
						else:
							self.chanHeaders[iChan][cH] = False
					else: # a bool or an integer
						self.chanHeaders[iChan][cH] = bool(self.chanHeaders[iChan][cH])

	# set some values
	def initialise(self):
		# create the type - index map
		self.chans = []
		self.chanMap = {}
		for iChan in xrange(0, self.getNumChannels()):
			chanType = self.chanHeaders[iChan]["channel_type"]		
			self.chanMap[chanType] = iChan 
			self.chans.append(chanType)

		# check the number of samples of each channel
		numSamples = []
		for c in self.chans:
			numSamples.append(self.getChanSamples(c))

		self.headers["num_samples"] = min(numSamples)
		#self.printText("Minimum number of samples = {:d}".format(self.getNumSamples()))
		for c, n in zip(self.chans, numSamples):
			if n != self.getNumSamples():
				self.printWarning("Not all channels have the same number of samples")
				self.printWarning("{} has {:d} samples more than the minimum".format(c, n - self.getNumSamples()))

		# create datetime objects
		datetimeStart = "{} {}".format(self.getHeaderVal("start_date"), self.getHeaderVal("start_time"))
		datetimeStop = "{} {}".format(self.getHeaderVal("stop_date"), self.getHeaderVal("stop_time"))
		self.datetimeStart = datetime.strptime(datetimeStart, "%Y-%m-%d %H:%M:%S.%f")
		self.datetimeStop = datetime.strptime(datetimeStop, "%Y-%m-%d %H:%M:%S.%f")

		# lastly, check the stop time - sometimes this appears to be one sample after the number of samples (ATS)
		# do a calculation and amend appropriately
		# this is the internal convention - start and end times should reflect the times of the first and last sample
		startTime, endTime = self.sample2time(0, self.getNumSamples()-1)
		if endTime != self.datetimeStop:		
			self.datetimeStop = endTime
			self.headers["stop_date"] = self.datetimeStop.strftime("%Y-%m-%d")
			self.headers["stop_time"] = self.datetimeStop.strftime("%H:%M:%S.%f")
			for idx, chan in enumerate(self.chanHeaders):
				self.chanHeaders[idx]["stop_date"] = self.datetimeStop.strftime("%Y-%m-%d")
				self.chanHeaders[idx]["stop_time"] = self.datetimeStop.strftime("%H:%M:%S.%f")

	###################
	### Error checking functions
	##################		
	def checkChan(self, chan):
		if chan not in self.chans:
			"Error - Channel does not exist"

	###################
	### DEBUG
	##################		
	# print headers
	def printInfo(self):		
		# print the headers
		self.printInfoBegin()
		self.printText("Data Path = {}".format(self.getDataPath()))
		self.printText("Global Headers")
		self.printText(self.headers)
		self.printText("Channels found:")
		self.printText(self.chans)
		self.printText("Channel Map")
		self.printText(self.chanMap)
		self.printText("Channel Headers")
		for c in self.chans:		
			self.printText(c)
			self.printText(self.chanHeaders[self.chanMap[c]])
		self.printText("Note: Field units used. Physical data has units mV/km for electric fields and mV for magnetic fields")
		self.printText("Note: To get magnetic field in nT, please calibrate")	
		self.printInfoEnd()

	def printInfoBegin(self):
		self.printText("####################")	
		self.printText("DATA READER INFO BEGIN")		
		self.printText("####################")	

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("DATA READER INFO END")		
		self.printText("####################")			
		
	def printText(self, infoStr):
		generalPrint("Data Reader Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Data Reader Warning", warnStr)



