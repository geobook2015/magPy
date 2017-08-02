# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
"""
import os
import glob
import re, struct
import xml.etree.ElementTree as ET
import collections
import copy
from datetime import datetime, timedelta
import numpy as np
# import dataReader
from dataReader import DataReader
# import the dataWriter for Reformatting
from dataWriterInternal import DataWriterInternal
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
		# set the sample byte size
		self.sampleByteSize = 3 # two's complement
		self.tagByteSize = 32
		self.dtype = int
		# there will be multiple TS files in here
		# need to figure out
		self.numHeaderFiles = len(self.headerF)
		self.numDataFiles = len(self.dataF)


	###################
	### GET METHODS
	### Some additional get methods for phoenix data
	##################
	def getSamplesRatesTS(self):
		info = {}
		for num, sr in zip(self.tsNums, self.tsSampleFreqs):
			info[num] = sr
		return info

	def getNumberSamplesTS(self):
		info = {}
		for num, ns in zip(self.tsNums, self.tsNumSamples):
			info[num] = ns
		return info

	###################
	### COUNT DATA
	##################
	# this only returns the data from the continuous channel
	def getUnscaledSamples(self, **kwargs):
		# initialise chans, startSample and endSample with the whole dataset
		options = self.parseGetDataKeywords(kwargs)

		# get the files to read and the samples to take from them, in the correct order
		recordsToRead, samplesToRead = self.getRecordsForSamples(options["startSample"], options["endSample"])
		# set up the dictionary to hold the data
		data = {}
		for chan in options["chans"]:
			data[chan] = np.zeros(shape=(options["endSample"] - options["startSample"] + 1), dtype=self.dtype)

		# open the file
		dFile = open(self.continuousF, "rb")

		# loop through chans and get data
		sampleCounter = 0
		for record, sToRead in zip(recordsToRead, samplesToRead):
			# number of samples to read in record
			dSamples = sToRead[1] - sToRead[0] + 1
			# find the byte read start and byte read end
			recordByteStart = self.recordBytes[self.continuous][record]
			recordSampleStart = self.recordSampleStarts[self.continuous][record]
			# find the offset on the readFrom bytes
			# now recall, each sample is recorded as a scan (all channels recorded at the same time)
			# so multiply by number of channels to get the number of bytes to read
			byteReadStart = recordByteStart + (sToRead[0] - recordSampleStart)*self.sampleByteSize*self.getNumChannels()
			bytesToRead = dSamples*self.sampleByteSize*self.getNumChannels()
			# read the data - numpy does not support 24 bit two's complement (3 bytes) - hence use struct
			dFile.seek(byteReadStart, 0) # seek to start byte from start of file
			dataBytes = dFile.read(bytesToRead)
			#dataBytes = struct.unpack("{}s".format(bytesToRead), dFile.read(bytesToRead))[0]
			dataRead = self.twosComplement(dataBytes)
			# now need to unpack this
			for chan in options["chans"]:
				# check to make sure channel exists
				self.checkChan(chan)
				# get the channel index - the chanIndex should give the right order in the data file
				# as it is the same order as in the header file
				chanIndex = self.chanMap[chan]
				# now populate the channel data appropriately
				data[chan][sampleCounter : sampleCounter + dSamples] = dataRead[chanIndex:dSamples*self.getNumChannels():self.getNumChannels()]
			# increment sample counter
			sampleCounter = sampleCounter + dSamples # get ready for the next data read

		# close file and return the data
		dFile.close()
		return data

	def getRecordsForSamples(self, startSample, endSample):
		# this reads the data from the continuous recording only
		recordsToRead = []
		samplesToRead = []
		for record, timeStart in enumerate(self.recordStarts[self.continuous]):
			recordStartSamp = self.recordSampleStarts[self.continuous][record]
			recordEndSamp = self.recordSampleStops[self.continuous][record]
			if recordStartSamp > endSample or recordEndSamp < startSample:
				continue # nothing to read from this file
			# in this case, there is some overlap with the samples to read
			recordsToRead.append(record)
			readFrom = recordStartSamp # i.e. the first sample in the datafile
			readTo = recordEndSamp # this the last sample in the file
			if recordStartSamp < startSample:
				readFrom = startSample
			if recordEndSamp > endSample:
				readTo = endSample
			# this is an inclusive number readFrom to readTo including readTo
			samplesToRead.append([readFrom, readTo])
		return recordsToRead, samplesToRead

	def readTag(self, dataFile):
	    second = struct.unpack('b', dataFile.read(1))[0]
	    minute = struct.unpack('b', dataFile.read(1))[0]
	    hour = struct.unpack('b', dataFile.read(1))[0]
	    day = struct.unpack('b', dataFile.read(1))[0]
	    month = struct.unpack('b', dataFile.read(1))[0]
	    year = struct.unpack('b', dataFile.read(1))[0]
	    dayOfWeek = struct.unpack('b', dataFile.read(1))[0]
	    century = struct.unpack('b', dataFile.read(1))[0]
	    dateString = "{:02d}{:02d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}.000".format(century, year, month, day, hour, minute, second)
	    # serial number
	    serialNum = struct.unpack('h', dataFile.read(2))
	    # num scans
	    numScans = struct.unpack('h', dataFile.read(2))[0]
	    # channels per scan
	    numChans = struct.unpack('b', dataFile.read(1))[0]
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

	    return numScans, numChans, dateString
		# return numScans, numChans

	def readRecord(self, dataFile, numChans, numScans):
	    data = np.zeros(shape=(numChans, numScans), dtype='int')
	    for scan in xrange(0, numScans):
	        for chan in xrange(0, numChans):
	            dataBytes = dataFile.read(3)
	            data[chan, scan] = twosComplement(dataBytes)
	    return data

	def twosComplement(self, dataBytes):
	    # two's complement 24-bit integer, little endian, unsigned and signed
		# this is padding 3 bytes out with a null byte
		# and reading as unsigned integer with little endian (<)
		numSamples = len(dataBytes)/self.sampleByteSize # this should be exact
		dataRead = np.zeros(shape=(numSamples), dtype=self.dtype)
		for i in xrange(0, numSamples):
			sampleBytes = dataBytes[ i*self.sampleByteSize : (i+1)*self.sampleByteSize ]
			unsigned = struct.unpack("<I", sampleBytes + "\x00")[0]
			signed = unsigned if not (unsigned & 0x800000) else unsigned - 0x1000000
			dataRead[i] = signed
		return dataRead


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
		# need to remove the gain
		for chan in options["chans"]:
			# apply the lsb
			# remove the gain
			data[chan] = data[chan]/self.getChanGain1(chan)
			# divide by distance in km
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
		self.continuousF = self.dataF[self.continuousI]
		# read the table data
		self.tableData = self.readTable()
		# and then populate the headers
		self.headers, self.chanHeaders = self.headersFromTable(self.tableData)
		# finally, check the number of samples in each file
		self.checkSamples()

	def readTable(self):
		if len(self.headerF) > 1:
			self.printWarning("More table files than expected. Using: {}".format(self.headerF[0]))
		tableFile = open(self.headerF[0], "rb")
		tableData = collections.OrderedDict()
		# why 138 - because this seems to be how many header words there are
		for i in xrange(0, 138):
			# formats for reading in
		    ints1 = ["EGN", "EGNC", "HGN", "HGNC", "MTSR", "L2NS", "L3NS", "L4NS", "TXPR", "TBVO", "TBVI", "INIT", "RQST", "MODE", "AQST", "HSMP", "CCLS",
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
		        value = struct.unpack('8b', tableFile.read(8))
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
		self.tsSampleFreqs = []
		for tsNum in self.tsNums:
			self.tsSampleFreqs.append(tableData["SRL{}".format(tsNum)])
		# for sample frequency, use the continuous channel
		headers["sample_freq"]	= self.tsSampleFreqs[self.continuousI]
		# these are the unix time stamps
		firstDate, firstTime, lastDate, lastTime = self.getDates(tableData)
		# the start date is equal to the time of the first record
		headers["start_date"] = firstDate
		headers["start_time"] = firstTime
		datetimeStart = datetime.strptime("{} {}".format(firstDate, firstTime), "%Y-%m-%d %H:%M:%S.%f")
		# the stop date
		datetimeLast = datetime.strptime("{} {}".format(lastDate, lastTime), "%Y-%m-%d %H:%M:%S.%f")
		# records are usually equal to one second (beginning on 0 and ending on the last sample before the next 0)
		datetimeStop = datetimeLast + timedelta(seconds=(1.0-1.0/headers["sample_freq"]))
		# put the stop date and time in the headers
		headers["stop_date"] = datetimeStop.strftime("%Y-%m-%d")
		headers["stop_time"] = datetimeStop.strftime("%H:%M:%S.%f")
		# here calculate number of samples
		deltaSeconds = (datetimeStop - datetimeStart).total_seconds()
		# calculate number of samples - have to add one because the time given in SPAM recording is the actual time of the last sample
		numSamples = round(deltaSeconds*headers["sample_freq"]) + 1
		headers["num_samples"] = numSamples
		headers["ats_data_file"] = self.continuousF
		# deal with the channel headers
		# now want to do this in the correct order
		# chan headers should reflect the order in the data
		chans = ["Ex", "Ey", "Hx", "Hy", "Hz"]
		chanOrder = []
		for chan in chans:
			chanOrder.append(tableData["CH{}".format(chan.upper())])
		# sort the lists in the right order based on chanOrder
		chanOrder, chans = (list(x) for x in zip(*sorted(zip(chanOrder, chans), key=lambda pair: pair[0])))
		for chan in chans:
			chanH = self.chanDefaults()
			# set the sample frequency from the main headers
			chanH["sample_freq"] = headers["sample_freq"]
			# channel output information (sensor_type, channel_type, ts_lsb, pos_x1, pos_x2, pos_y1, pos_y2, pos_z1, pos_z2, sensor_sernum)
			chanH["ats_data_file"] = self.dataF[self.continuousI]
			chanH["num_samples"] = numSamples
			# channel information
			chanH["channel_type"] = consistentChans(chan) # consistent chan naming

			# magnetic channels only
			if isMagnetic(chanH["channel_type"]):
				chanH["sensor_sernum"] = tableData["{}SN".format(chan.upper())][-4:]
				chanH["sensor_type"] = "Phoenix"
				# channel input information (gain_stage1, gain_stage2, hchopper, echopper)
				chanH["gain_stage1"] = tableData["HGN"]
				chanH["gain_stage2"] = 1

			# electric channels only
			if isElectric(chanH["channel_type"]):
				# the distances
				if chan == "Ex":
					chanH["pos_x1"] = float(tableData["EXLN"])/2.0
					chanH["pos_x2"] = chanH["pos_x1"]
				if chan == "Ey":
					chanH["pos_y1"] = float(tableData["EYLN"])/2.0
					chanH["pos_y2"] = chanH["pos_y1"]
				# channel input information (gain_stage1, gain_stage2, hchopper, echopper)
				chanH["gain_stage1"] = tableData["EGN"]
				chanH["gain_stage2"] = 1

			# append chanHeaders to the list
			chanHeaders.append(chanH)

		# data information (meas_channels, sample_freq)
		headers["meas_channels"] = len(chans) # this gets reformatted to an int later
		# return the headers and chanHeaders from this file
		return headers, chanHeaders

	def getDates(self, tableData):
		# this is the start time of the first record
		firstSecond = tableData["FTIM"][0]
		firstMinute = tableData["FTIM"][1]
		firstHour = tableData["FTIM"][2]
		firstDay = tableData["FTIM"][3]
		firstMonth = tableData["FTIM"][4]
		firstYear = tableData["FTIM"][5]
		firstCentury = tableData["FTIM"][-1]
		firstDate = "{:02d}{:02d}-{:02d}-{:02d}".format(firstCentury, firstYear, firstMonth, firstDay)
		firstTime = "{:02d}:{:02d}:{:02d}.000".format(firstHour, firstMinute, firstSecond)
		# this is the start time of the last record
		lastSecond = tableData["LTIM"][0]
		lastMinute = tableData["LTIM"][1]
		lastHour = tableData["LTIM"][2]
		lastDay = tableData["LTIM"][3]
		lastMonth = tableData["LTIM"][4]
		lastYear = tableData["LTIM"][5]
		lastCentury = tableData["LTIM"][-1]
		lastDate = "{:02d}{:02d}-{:02d}-{:02d}".format(lastCentury, lastYear, lastMonth, lastDay)
		lastTime = "{:02d}:{:02d}:{:02d}.000".format(lastHour, lastMinute, lastSecond)
		return firstDate, firstTime, lastDate, lastTime

	# check the number of samples for all the ts files
	# recall, the format is 3 bytes two's complement per sample
	def checkSamples(self):
		self.recordStarts = {}
		self.recordScans = {}
		self.recordBytes = {}
		self.recordSampleStarts = {}
		self.recordSampleStops = {}
		# loop over the tsNums
		samplesDict = {}
		for dFileName in self.dataF:
			ts = int(dFileName[-1])
			self.recordStarts[ts] = []
			self.recordScans[ts] = []
			self.recordBytes[ts] = []
			self.recordSampleStarts[ts] = []
			self.recordSampleStops[ts] = []
			# logFile = open("log{}.txt".format(ts), "w")
			# start number of samples at 0
			samples = 0
			# get file size in samples
			numBytes = os.path.getsize(dFileName)
			bytesread = 0
			# now run through the file and figure out the number of samples
			dFile = open(dFileName, "rb")
			while bytesread < numBytes:
				numScans, numChans, dateString = self.readTag(dFile) # this is 32 bytes
				self.recordBytes[ts].append(bytesread + self.tagByteSize)
				dataBytes = numScans*numChans*self.sampleByteSize
				dFile.seek(dataBytes, 1)
				bytesread += self.tagByteSize + dataBytes
				# save the record start times and scan lengths
				self.recordStarts[ts].append(dateString)
				self.recordScans[ts].append(numScans)
				# save the sample starts
				self.recordSampleStarts[ts].append(samples)
				# increment the number of samples
				# recall, a scan is all channels recorded at one time
				# this is equivalent to one sample
				samples += numScans # this is the count
				# sample stop is samples -1 because inclusive of the current sample
				self.recordSampleStops[ts].append(samples-1)
				# logFile.write("{} : {} : {} - {}\n".format(dateString, numScans, self.recordSampleStarts[ts][-1], self.recordSampleStops[ts][-1]))
			dFile.close()
			# save number of samples in dict
			samplesDict[ts] = samples
			# logFile.close()

		self.tsNumSamples = []
		for tsNum in self.tsNums:
			self.tsNumSamples.append(samplesDict[tsNum])

		# check the samples of the continuous file
		if self.tsNumSamples[self.continuousI] != self.getNumSamples():
			self.printWarning("Number of samples calculated from times is different to that in file")
			self.printWarning("{} samples in file, {} calculated from time".format(self.tsNumSamples[self.continuousI], self.getNumSamples()))


	###################
	### REFORMAT DATA TO INTERNAL FORMAT
	###################
	def reformatHigh(self, path):
		writer = DataWriterInternal()
		for idx, ts in enumerate(self.tsNums):
			# let's get the headers
			headers = self.getHeaders()
			chanHeaders, chanMap = self.getChanHeaders()
			chans = self.getChannels()
			# now go through the different ts files to get ready to output
			if ts == self.continuous:
				continue
			sampleFreq = self.tsSampleFreqs[idx]
			# set sample frequency in headers
			headers["sample_freq"] = sampleFreq
			for cH in chanHeaders:
				cH["sample_freq"] = sampleFreq
			# now open the data file
			dFile = open(self.dataF[idx], "rb")
			# each record has to be read separately and then compare time to previous
			outStartTime = datetime.strptime(self.recordStarts[ts][0], "%Y-%m-%d %H:%M:%S.%f")
			# set up the data dictionary
			data = {}
			for record, startDate in enumerate(self.recordStarts[ts]):
				# start date is a string
				startByte = self.recordBytes[ts][record]
				startDateTime = datetime.strptime(startDate, "%Y-%m-%d %H:%M:%S.%f")
				# read the record - numpy does not support 24 bit two's complement (3 bytes) - hence use struct
				bytesToRead = self.recordScans[ts][record]*self.sampleByteSize*self.getNumChannels()
				dFile.seek(startByte, 0) # seek to start byte from start of file
				dataBytes = dFile.read(bytesToRead)
				dataRead = self.twosComplement(dataBytes)
				dataRecord = {}
				for chan in chans:
					# as it is the same order as in the header file
					chanIndex = self.chanMap[chan]
					dataRecord[chan] = dataRead[chanIndex:self.recordScans[ts][record]*self.getNumChannels():self.getNumChannels()]
				# need to compare to previous record
				if record != 0 and startDateTime != prevEndTime:
					# then need to write out the current data before saving the new data
					# write out current data
					outStopTime = prevEndTime - timedelta(seconds=1.0/sampleFreq) # because inclusive of first sample (previous end time for continuity comparison)
					# calculate number of samples
					numSamples = data[chans[0]].size
					headers["start_date"] = outStartTime.strftime("%Y-%m-%d")
					headers["start_time"] = outStartTime.strftime("%H:%M:%S.%f")
					headers["stop_date"] = outStopTime.strftime("%Y-%m-%d")
					headers["stop_time"] = outStopTime.strftime("%H:%M:%S.%f")
					headers["num_samples"] = numSamples
					for cH in chanHeaders:
						cH["start_date"] = headers["start_date"]
						cH["start_time"] = headers["start_time"]
						cH["stop_date"] = headers["stop_date"]
						cH["stop_time"] = headers["stop_time"]
						cH["num_samples"] = numSamples
					# get the outpath
					dataOutpath = os.path.join(path, "meas_ts{}_{}_{}".format(ts, outStartTime.strftime("%Y_%m_%d_%H_%M_%S"), outStopTime.strftime("%Y_%m_%d_%H_%M_%S")))
					# write out
					writer.setOutPath(dataOutpath)
					writer.writeData(headers, chanHeaders, data)
					# then save current data
					outStartTime = startDateTime
					data = copy.deepcopy(dataRecord)
					prevEndTime = startDateTime + timedelta(seconds=((1.0/sampleFreq)*self.recordScans[ts][record]))
				else:
					# then record == 0 or startDateTime == prevEndTime
					# update prevEndTime
					prevEndTime = startDateTime + timedelta(seconds=((1.0/sampleFreq)*self.recordScans[ts][record]))
					if record == 0:
						data = copy.deepcopy(dataRecord)
						continue
					# otherwise, want to concatenate the data
					for chan in chans:
						data[chan] = np.concatenate((data[chan], dataRecord[chan]))
			# close the data file
			dFile.close()

	def reformatContinuous(self, path):
		writer = DataWriterInternal()
		outpath = "meas_ts{}_{}_{}".format(self.continuous, self.getStartDatetime().strftime("%Y_%m_%d_%H_%M_%S"), self.getStopDatetime().strftime("%Y_%m_%d_%H_%M_%S"))
		outpath = os.path.join(path, outpath)
		writer.setOutPath(outpath)
		headers = self.getHeaders()
		chanHeaders, chanMap = self.getChanHeaders()
		writer.writeData(headers, chanHeaders, self.getUnscaledSamples())

	def reformat(self, prepend):
		self.reformatContinuous(prepend)
		self.reformatHigh(prepend)


	###################
	### DEBUG
	##################
	def printInfoBegin(self):
		self.printText("####################")
		self.printText("PHOENIX READER INFO BEGIN")
		self.printText("NOTE: Information given for continuous recording file")
		self.printText("####################")

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("PHOENIX READER INFO BEGIN")
		self.printText("####################")

	def printDataFileList(self):
		self.printText("####################")
		self.printText("PHOENIX READER DATA FILE LIST BEGIN")
		self.printText("####################")
		self.printText("TS File\t\tSampling frequency (Hz)\t\tNum Samples")
		for dF, tsF, tsN in zip(self.dataF, self.tsSampleFreqs, self.tsNumSamples):
			self.printText("{}\t\t{}\t\t{}".format(os.path.basename(dF), tsF, tsN))
		self.printText("####################")
		self.printText("PHOENIX READER DATA FILE LIST END")
		self.printText("####################")

	def printTableFile(self):
		self.printText("####################")
		self.printText("PHOENIX READER TABLE FILE BEGIN")
		self.printText("####################")
		for h, v in self.tableData.items():
			self.printText("{} = {}".format(h,v))
		self.printText("####################")
		self.printText("PHOENIX READER TABLE FILE END")
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Phoenix Reader Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Phoenix Reader Warning", warnStr)
