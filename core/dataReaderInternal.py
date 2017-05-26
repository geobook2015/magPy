# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
"""
import os
import glob
from datetime import datetime, timedelta
import numpy as np
# import dataReaderATS
# inherit from dataReaderATS because the binary format is the same
# the headers are kept in plain text format
from dataReader import DataReader
# utils
from utilsIO import *

class DataReaderInternal(DataReader):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, dataPath):
		self.dataPath = dataPath
		# get a list of the xml files in the folder
		self.headerF = glob.glob(os.path.join(self.dataPath,"*.hdr"))
		self.dataF = glob.glob(os.path.join(self.dataPath,"*.dat"))
		self.dataByteOffset = 0
		self.dataByteSize = 4
		self.readHeader()
		self.formatHeaderData()
		self.initialise()	

	###################
	### READ HEADER
	##################
	# the headers in the header file
	def dataHeaders():
		# note, in comparison to ats headers, this also has one called lsb_applied
		recordingHeaders = ["start_time", "start_date", "stop_time", "stop_date"]
		globalHeaders = ["meas_channels", "sample_freq"]
		channelHeadersInput = ["gain_stage1", "gain_stage2", "hchopper", "echopper"]
		channelHeadersOutput = ["start_time", "start_date", "sample_freq", "num_samples", "ats_data_file", "sensor_type", "channel_type", "ts_lsb", "lsb_applied", "pos_x1", "pos_x2", "pos_y1", "pos_y2", "pos_z1", "pos_z2", "sensor_sernum"]
		return recordingHeaders, globalHeaders, channelHeadersInput, channelHeadersOutput
	
	# read in the header file
	def readHeader(self):
		# first read the global hears
		# look in headerF for global.hdr
		if os.path.join(self.dataPath, "global.hdr") not in self.headerF:
			 self.printWarning("Global header not found. The global.hdr file is required")
		globalF = open(os.path.join(self.dataPath, "global.hdr"))
		lines = globalF.readlines()
		globalF.close()
		# ignore the first line
		lines.pop(0)
		# now go through and get the headers
		self.headers = {}
		for l in lines:
			if l == "":
				continue
			key, val = self.lineToKeyAndValue(l.strip())
			self.headers[key] = val
		
		# now want to deal with the chan headers
		numChans = int(self.headers["meas_channels"])
		self.chanHeaders = []
		for iChan in xrange(0, numChans):
			chanF = open(os.path.join(self.dataPath, "chan_{:02d}.hdr".format(iChan)))
			lines = chanF.readlines()
			chanF.close()
			# remove first line and read the headers for the channel
			lines.pop(0)
			cHeader = {}
			for l in lines:
				if l == "":
					continue
				key, val = self.lineToKeyAndValue(l.strip())
				cHeader[key] = val
			self.chanHeaders.append(cHeader)

	# a helper function to read the headers
	def lineToKeyAndValue(self, line):
		split = line.split("=")
		key = split[0].strip()
		val = split[1].strip()
		return key, val

	###################
	### DEBUG
	##################		
	def printInfoBegin(self):
		self.printText("####################")	
		self.printText("INTERNAL READER INFO BEGIN")		
		self.printText("####################")	

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("INTERNAL READER INFO END")		
		self.printText("####################")			
		
	def printText(self, infoStr):
		generalPrint("INTERNAL Reader Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("INTERNAL Reader Warning", warnStr)