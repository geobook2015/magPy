# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
"""
import os
import glob
import re, struct
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import numpy as np
# import dataReader
from dataReader import DataReader
# utils
from utilsIO import *
from utilsChecks import *
from utilsProcess import removeZeros, removeZerosSingle

class DataReaderPhoenix(DataReader):
	"""
	The Phoenix data and recording format is different. There are three frequencies recorded consecutively (e.g. 2400Hz, 150Hz, 15Hz)
	The lowest sampling frequency is continuous whilst the others record small bits at regular intervals
	There is no issue with the continous sampling frequency. However, with the others, this is going to lead to lots of small data folders
	But this is the only way to include it into the code and still make it interoperable with other data formats
	"""

	###################
	### PREPARE
	##################
	def prepare(self):
		# get a list of the header and data files in the folder
		self.headerF = glob.glob(os.path.join(self.dataPath,"*.TBL"))
		self.dataF = glob.glob(os.path.join(self.dataPath,"*.TS*"))
		print self.headerF
		print self.dataF
		# there will be multiple TS files in here
		# need to figure out
		self.numHeaderFiles = len(self.headerF)
		self.numDataFiles = len(self.dataF)

	def initialise():
		return

	###################
	### COUNT DATA
	##################
	# this only returns the data from the continuous channel
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

	def readTag(dataFile):
	    second = struct.unpack('b', dataFile.read(1))[0]
	    minute = struct.unpack('b', dataFile.read(1))[0]
	    hour = struct.unpack('b', dataFile.read(1))[0]
	    day = struct.unpack('b', dataFile.read(1))[0]
	    month = struct.unpack('b', dataFile.read(1))[0]
	    year = struct.unpack('b', dataFile.read(1))[0]
	    dayOfWeek = struct.unpack('b', dataFile.read(1))[0]
	    century = struct.unpack('b', dataFile.read(1))[0]
	    print "{:02d}:{:02d}:{:02d} {:02d}/{:02d}/{:02d}".format(hour, minute, second, day, month, year)
	    # serial number
	    serialNum = struct.unpack('h', dataFile.read(2))
	    # num scans
	    numScans = struct.unpack('h', dataFile.read(2))[0]
	    print "numScans = {}".format(numScans)
	    # channels per scan
	    numChans = struct.unpack('b', dataFile.read(1))[0]
	    print "numChans = {}".format(numChans)
	    # tag length
	    tagLength = struct.unpack('b', dataFile.read(1))
	    # status code
	    statusCode = struct.unpack('b', dataFile.read(1))
	    # bit-wise saturation flags
	    saturationFlag = struct.unpack('b', dataFile.read(1))
	    # reserved
	    reserved = struct.unpack('b', dataFile.read(1))
	    # sample length
	    sampleLength = struct.unpack('b', dataFile.read(1))
	    # sample rate
	    sampleRate = struct.unpack('h', dataFile.read(2))
	    # units of sample rate: 0 = Hz, 1 = minute, 2 = hour, 3 = day
	    sampleUnits = struct.unpack('b', dataFile.read(1))
	    # clock status
	    clockStatus = struct.unpack('b', dataFile.read(1))
	    # clock error in micro seconds
	    clockError = struct.unpack('i', dataFile.read(4))
	    # reserved
	    res1 = struct.unpack('b', dataFile.read(1))
	    res2 = struct.unpack('b', dataFile.read(1))
	    res3 = struct.unpack('b', dataFile.read(1))
	    res4 = struct.unpack('b', dataFile.read(1))
	    res5 = struct.unpack('b', dataFile.read(1))
	    res6 = struct.unpack('b', dataFile.read(1))

	    return numScans, numChans

	def readRecord(dataFile, numChans, numScans):
	    print "samples to read = {}".format(numChans*numScans)
	    data = np.zeros(shape=(numChans, numScans), dtype='int')
	    for scan in xrange(0, numScans):
	        for chan in xrange(0, numChans):
	            dataBytes = dataFile.read(3)
	            data[chan, scan] = twosComplement(dataBytes)
	    return data

	###################
	### PHYSICAL DATA
	##################
	# get data in physical units
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

	def twosComplement(dataBytes):
	    # two's complement 24-bit integer, little endian, unsigned and signed
	    unsigned = struct.unpack('<I', dataBytes + '\x00')[0]
	    signed = unsigned if not (unsigned & 0x800000) else unsigned - 0x1000000
	    return signed

	###################
	### READ HEADERS
	###################
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

	# extend a method in the base class
	# read in header (table) file
	# this is binary
	def readHeader(self):
		# first, find which ts files are available (2,3,4,5)
		# and the continuous recording frequency (the max)
		self.tsNums = []
		for tsfile in self.dataF:
			self.tsNums.append(int(tsfile[-1]))
		self.continuous = max(self.tsNums)
		self.continuousI = self.tsNums.index(self.continuous)
		print self.tsNums, self.continuous
		# read the table data
		tableData = self.readTable()
		# and then populate the headers
		headers, chanHeaders = self.headersFromTable(tableData)
		# finally, check the number of samples in each file
		self.checkSamples()

	def readTable(self):
		if len(self.headerF) > 1:
			self.printWarning("More table files than expected. Using: {}".format(self.headerF[0]))
		tableFile = open(self.headerF[0], "rb")
		tableData = {}
		# why 138 - because this seems to be how many header words there are
		for i in xrange(0, 138):
			# formats for reading in
		    ints1 = ["EGN", "EGNC", "HGN", "MTSR", "L2NS", "L3NS", "L4NS", "TXPR", "TBVO", "TBVI", "INIT", "RQST", "MODE", "AQST", "HSMP", "CCLS",
		        "TEMP", "TMAX", "CHEX", "CHEY", "CHHX", "CHHY", "CHHZ", "NREF", "CCLT", "PZLT", "NSAT", "OCTR", "CLST", "TALS", "TCMB", "TOTL"
		    ]
		    ints2 = ["SNUM", "MXSC", "SATR", "BADR", "BAT1","BAT2","BAT3", "EXR", "EYR", "ELEV", "SRL2", "SRL3", "SRL4", "SRL5", "LPFR", "LFRQ"]
		    ints4 = ["STDE", "STDH"]
		    ints8 = []
		    ints1_4 = []
		    ints1_8 = ["TDSP", "LFIX", "TSYN", "STIM", "ETIM", "HTIM", "ETMH", "NUTC", "FTIM", "LTIM"]
		    floats = []
		    doubles = ["EXAC", "EXDC", "EYAC", "EYDC", "HXAC", "HXDC", "HYAC", "HYDC", "HZAC", "HZDC", "EXNR", "EXPR", "EYNR", "EYPR", "GNDR",
		        "TSTV", "FSCV", "CCMN", "CCMX", "HATT", "HAMP", "LFIX", "EXLN", "EYLN", "TSTR", "INPR", "CFMN", "CFMX", "HNOM"
		    ]
		    # get the header word
		    header = struct.unpack('12s', tableFile.read(12))
		    header = self.removeControl(header[0])
		    if header in ints1:
		        value = struct.unpack('b', tableFile.read(1))[0]
		        tableFile.read(12)
		    elif header in ints2:
		        value = struct.unpack('h', tableFile.read(2))[0]
		        tableFile.read(11)
		    elif header in ints4:
		        value = struct.unpack('i', tableFile.read(4))[0]
		        tableFile.read(9)
		    elif header in ints8:
		        value = struct.unpack('q', tableFile.read(8))[0]
		        tableFile.read(5)
		    elif header in ints1_4:
		        value = struct.unpack('4b', tableFile.read(4))[0]
		        tableFile.read(9)
		    elif header in ints1_8:
		        value = struct.unpack('8b', tableFile.read(8))[0]
		        tableFile.read(5)
		    elif header in floats:
		        value = struct.unpack('f', tableFile.read(4))[0]
		        tableFile.read(9)
		    elif header in doubles:
		        value = struct.unpack('d', tableFile.read(8))[0]
		        tableFile.read(5)
		    else:
		        value = struct.unpack('13s', tableFile.read(13))
		        value = self.removeControl(value[0])
		    tableData[header] = value
		for h, v in tableData.items():
			print "{} = {}".format(h,v)
		tableFile.close()
		return tableData

	def removeControl(self, inString):
	    # return re.sub(r'[\x00-\x1f\x7f-\x9f]', '', inString)
	    return re.sub(r'[\x00-\x05\xff]', '', inString)

	# populate the headers from the table values
	def headersFromTable(self, tableData):
		# initialise storage
		headers = {}
		chanHeaders = []
		# get the sample freqs for each ts file
		self.tsSampleRates = []
		for tsNum in self.tsNums:
			self.tsSampleRates.append(tableData["SRL{}".format(tsNum)])
		print self.tsSampleRates
		# for sample frequency, use the continuous channel
		headers["sample_freq"]	= self.tsSampleRates[self.continuousI]
		# these are the unix time stamps
		# startDate = float(timeSplit[1] + "." + timeSplit[2])
		# datetimeStart = datetime.utcfromtimestamp(startDate)
		# stopDate = float(timeSplit[3] + "." + timeSplit[4])
		# datetimeStop = datetime.utcfromtimestamp(stopDate)
		# headers["start_date"] = datetimeStart.strftime("%Y-%m-%d")
		# headers["start_time"] = datetimeStart.strftime("%H:%M:%S.%f")
		# headers["stop_date"] = datetimeStop.strftime("%Y-%m-%d")
		# headers["stop_time"] = datetimeStop.strftime("%H:%M:%S.%f")
		# # here calculate number of samples
		# deltaSeconds = (datetimeStop - datetimeStart).total_seconds()
		# # calculate number of samples - have to add one because the time given in SPAM recording is the actual time of the last sample
		# numSamples = int(deltaSeconds*headers["sample_freq"]) + 1
		# put these in headers for ease of future calculations in merge headers
		numSamples = 1
		headers["num_samples"] = numSamples
		headers["ats_data_file"] = self.dataF[self.continuousI]
		# deal with the channel headers
		# now want to do this in the correct order
		# chan headers should reflect the order in the data
		chans = ["Ex", "Ey", "Hx", "Hy", "Hz"]
		chanOrder = []
		for chan in chans:
			chanOrder.append(tableData["CH{}".format(chan.upper())])
		print chans, chanOrder
		# sort the lists in the right order based on chanOrder
		chanOrder, chans = (list(x) for x in zip(*sorted(zip(chanOrder, chans), key=lambda pair: pair[0])))
		print chans, chanOrder
		for chan in chans:
			chanH = self.chanDefaults()
			# set the sample frequency from the main headers
			chanH["sample_freq"] = headers["sample_freq"]
			# channel input information (gain_stage1, gain_stage2, hchopper, echopper)
			chanH["gain_stage1"] = 1
			chanH["gain_stage2"] = 1
			# channel output information (sensor_type, channel_type, ts_lsb, pos_x1, pos_x2, pos_y1, pos_y2, pos_z1, pos_z2, sensor_sernum)
			chanH["ats_data_file"] = self.dataF[self.continuousI]
			chanH["num_samples"] = numSamples
			# channel information
			chanH["channel_type"] = consistentChans(chan) # consistent chan naming
			# serial numbers of coils
			if isMagnetic(chanH["channel_type"]):
				chanH["sensor_sernum"] = tableData["{}SN".format(chan.upper())]
				chanH["sensor_type"] = "Phoenix"

			# the distances
			if chan == "EX":
				chanH["pos_x1"] = float(tableData["EXLN"])/2.0
				chanH["pos_x2"] = chanH["pos_x1"]
			if chan == "EY":
				chanH["pos_y1"] = float(tableData["EYLN"])/2.0
				chanH["pos_y2"] = chanH["pos_y1"]

			# append chanHeaders to the list
			chanHeaders.append(chanH)

		# data information (meas_channels, sample_freq)
		headers["meas_channels"] = len(chans) # this gets reformatted to an int later

		# return the headers and chanHeaders from this file
		return headers, chanHeaders

		def checkSamples(self):
			# make sure the continuous number of samples are ok
			# and for the other sampling frequencies, just calculate the number of samples


	###################
	### REFORMAT DATA TO INTERNAL FORMAT
	###################
	# def reformat():
		# go through the data files and then reformat them to the internal data format



	###################
	### DEBUG
	##################
	def printInfoBegin(self):
		self.printText("####################")
		self.printText("PHOENIX READER INFO BEGIN")
		self.printText("####################")

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("PHOENIX READER INFO BEGIN")
		self.printText("####################")

	def printDataFileList(self):
		self.printText("####################")
		self.printText("PHOENIX READER DATA FILE LIST BEGIN")
		self.printText("####################")
		self.printText("Data File\t\tSample Ranges")
		for dFile, sRanges in zip(self.dataFileList, self.dataRanges):
			self.printText("{}\t\t{} - {}".format(dFile, sRanges[0], sRanges[1]))
		self.printText("Total samples = {}".format(self.getNumSamples()))
		self.printText("####################")
		self.printText("PHOENIX READER DATA FILE LIST END")
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Emerald Reader Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Emerald Reader Warning", warnStr)
