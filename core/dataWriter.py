# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop

class DataWriter
Parent class for data writers
NOTE: all data formats written have the same headers
Writing of ats style xml headers is not supported
"""
import os
import glob
import struct
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import numpy as np
# utils
from utilsIO import *

class DataWriter(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self):
		self.outPath = ""
		# in subclasses, extension might change
		# i.e. .ats
		self.extension = ".dat"
		# data type - the format of the data being written out
		self.dtype = np.int32
		# information about data being written
		self.headers = None
		self.chans = None
		self.chanMap = None
		self.chanHeaders = None


	###################
	### GET FUNCTIONS
	##################
	def getOutPath(self):
		return self.outPath


	###################
	### SET FUNCTIONS
	##################
	def setOutPath(self, path):
		self.outPath = path

	def setGlobalHeadersFromKeywords(self, headers, keywords):
		globalHeaderwords = self.globalHeaderwords()
		for gH in globalHeaderwords:
			hdrVal = ""
			if gH in headers:
				hdrVal = headers[gH]
			if gH in keywords:
				hdrVal = keywords[gH]
			headers[gH] = hdrVal
		return headers	

	def setChanHeadersFromKeywords(self, chanHeaders, keywords):
		chanHeaderwords = self.chanHeaderwords()
		for iChan in xrange(0, len(chanHeaders)):
			for cH in chanHeaderwords:
				hdrVal = ""
				if cH in chanHeaders[iChan]:
					hdrVal = chanHeaders[iChan][cH]
				if cH in keywords:
					hdrVal = keywords[cH]
				chanHeaders[iChan][cH] = hdrVal
		return chanHeaders


	###################
	### CALCULATION METHODS
	##################
	def calcStopDateTime(self, fs, numSamples, datetimeStart):
		# calculate duration in seconds
		duration = 1.0*(numSamples-1)/fs # numSamples - 1 because have to remove the initial sample which is taken at start time
		datetimeStop = datetimeStart + timedelta(seconds=duration)
		# now get this in start and stop times
		stopDate = datetimeStop.strftime("%Y-%m-%d") 	
		stopTime = datetimeStop.strftime("%H:%M:%S.%f")
		return datetimeStop, stopDate, stopTime


	###################
	### HEADERWORDS SUPPORTED
	###################
	def globalHeaderwords(self):
		gHeaders = ["sample_freq", "num_samples", "start_time", "start_date", "stop_time", "stop_date", "meas_channels"]
		return gHeaders

	def chanHeaderwords(self):
		cHeaders = ["sample_freq", "num_samples", "start_time", "start_date", "stop_time", "stop_date", "ats_data_file", "sensor_type", "channel_type","ts_lsb", "lsb_applied",
			  "pos_x1", "pos_x2", "pos_y1", "pos_y2", "pos_z1", "pos_z2", "sensor_sernum", "gain_stage1", "gain_stage2", "hchopper", "echopper"]
		return cHeaders	


	###################
	### WRITE DATA
	###################
	def writeDataset(self, reader, **kwargs):
		if self.getOutPath() == "":
			self.printWarning("No output filepath given")
			return
		# make the directory
		checkAndMakeDir(self.getOutPath())
		# write using information from a reader file
		headers = reader.getHeaders()
		chanHeaders, chanMap = reader.getChanHeaders()
		chans = reader.getChannels()
		# now write depending on whether lsb_applied or not
		if "lsb_applied" in kwargs and kwargs["lsb_applied"]:
			self.write(headers, chanHeaders, chanMap, reader.getPhysicalSamples(), **kwargs)
			self.dtype = np.float32
		else:
			self.write(headers, chanHeaders, chanMap, reader.getUnscaledSamples(), **kwargs)

	def writeData(self, headers, chanHeaders, data, **kwargs):
		if self.getOutPath() == "":
			self.printWarning("No output filepath given")
			return
		# make the directory
		checkAndMakeDir(self.getOutPath())
		# calculate our own cMap
		chanMap = {}
		for iChan in xrange(0, len(chanHeaders)):
			chanType = chanHeaders[iChan]['channel_type']		
			chanMap[chanType] = iChan 
		# check the data type
		if "lsb_applied" in kwargs and kwargs["lsb_applied"]:
			self.dtype = np.float32		
		# write the data
		self.write(headers, chanHeaders, chanMap, data, **kwargs)

	# write out the dataset
	def write(self, headers, chanHeaders, chanMap, data, **kwargs):
		# set global headers for keyword arguments
		headers = self.setGlobalHeadersFromKeywords(headers, kwargs)
		# set channel headers for keyword arguments
		chanHeaders = self.setChanHeadersFromKeywords(chanHeaders, kwargs)

		# now overwrite the options by checking the actual data
		# let's check all the data sizes
		chans = sorted(list(data.keys()))
		dataSizes = []
		for c in chans:
			dataSizes.append(data[c].size)
		if min(dataSizes) != max(dataSizes):
			self.printWarning("Channels do not have the same number of samples: {} - {}".format(", ".join(chans), ", ".join(dataSizes)))
			self.printWarning("Only the smallest number of samples will be written out")
		numSamples = min(dataSizes)

		# set number of samples from actual data
		headers["num_samples"] = numSamples

		# limit data and set the chan header 
		for c in chans:
			data[c] = data[c][:numSamples]
			cIndex = chanMap[c]
			chanHeaders[cIndex]["num_samples"] = numSamples

		# deal with start and end time
		# the start time does not change on resampling, only the end time
		duration = numSamples/headers["sample_freq"]
		# create datetime objects
		startString = '{} {}'.format(headers["start_date"], headers["start_time"])
		stopString = '{} {}'.format(headers["stop_date"], headers["stop_time"])
		datetimeStart = datetime.strptime(startString, "%Y-%m-%d %H:%M:%S.%f")
		datetimeStop = datetime.strptime(stopString, "%Y-%m-%d %H:%M:%S.%f")
		datetimeRecalc, stopDate, stopTime = self.calcStopDateTime(headers["sample_freq"], numSamples, datetimeStart)
		# compare to datetimeStop already
		if datetimeRecalc != datetimeStop:
			self.printWarning("Note, discrepancy between stop time in given headers and those calculated from data")
			self.printWarning("Causes of this might be resampling or interpolation processes and the limiting of data")
			self.printWarning("Stop time calculated from data will be used")
			self.printWarning("If no resampling, interpolation or limiting of data has been performed, please check all times")
		headers["stop_date"] = stopDate
		headers["stop_time"] = stopTime
		# do the same for the chan headers
		for c in chans:
			cIndex = chanMap[c]
			chanHeaders[cIndex]["stop_date"] = stopDate
			chanHeaders[cIndex]["stop_time"] = stopTime

		# finally, check the number of measurement channels 
		headers["meas_channels"] = len(chans) 

		# now write out the headers and save to class variables
		self.writeHeaders(headers, chans, chanMap, chanHeaders)
		self.headers = headers
		self.chans = chans
		self.chanMap = chanMap
		self.chanHeaders = chanHeaders

		# write out the data files
		self.writeDataFiles(chans, data)

	# write out the header files
	def writeHeaders(self, headers, chans, chanMap, chanHeaders):
		# write out the global headers
		f = open(os.path.join(self.getOutPath(), "global.hdr"), "w")
		f.write("HEADER = GLOBAL\n")
		globalHeaderwords = self.globalHeaderwords()
		for gH in globalHeaderwords:
			f.write("{} = {}\n".format(gH, headers[gH]))
		f.close()

		# write out the channel headers
		chanHeaderwords = self.chanHeaderwords()
		for idx, c in enumerate(chans):
			cf = open(os.path.join(self.getOutPath(), "chan_{:02d}.hdr".format(idx)), "w")
			cf.write("HEADER = CHANNEL\n")
			# now need to use the cMap to get the index of the cHeaders array
			cIndex = chanMap[c]
			# change the data file
			chanHeaders[cIndex]["ats_data_file"] = "chan_{:02d}{}".format(idx, self.extension)
			for cH in chanHeaderwords:
				cf.write("{} = {}\n".format(cH, chanHeaders[cIndex][cH]))
			cf.close()
		return True

	def writeDataFiles(self, chans, data):
		# implement in the child class
		return

	###################
	### DEBUG
	##################		
	# print headers
	def printInfo(self):		
		# print the headers
		self.printInfoBegin()
		self.printText("Output file path for data = {}".format(self.getOutPath()))
		# if it exists, print out the headers
		if self.headers:
			self.printText("Global Headers")
			self.printText(self.headers)
		# if exists, print out a list of chans
		if self.chans:
			self.printText("Channels found:")
			self.printText(self.chans)
		# if exists, print out the chanMap
		if self.chanMap:
			self.printText("Channel Map")
			self.printText(self.chanMap)
		# if it exists, print out the chanHeaders
		if self.chanHeaders:
			self.printText("Channel Headers")
			for c in self.chans:		
				self.printText(c)
				self.printText(self.chanHeaders[self.chanMap[c]])
		self.printInfoEnd()

	def printInfoBegin(self):
		self.printText("####################")	
		self.printText("DATA WRITER INFO BEGIN")		
		self.printText("####################")	

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("DATA WRITER INFO END")		
		self.printText("####################")			
		
	def printText(self, infoStr):
		generalPrint("Data Writer Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Data Writer Warning", warnStr)