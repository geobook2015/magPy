#!/usr/bin/python
import os
import sys
sys.path.append(os.path.join("..", "core"))
sys.path.append(os.path.join("..", "utils"))
import numpy as np
import math
from datetime import datetime, timedelta
# import readers
from dataReaderSpam import DataReaderSPAM
from dataReaderATS import DataReaderATS
from dataReaderInternal import DataReaderInternal
# import writers
from dataWriterInternal import DataWriterInternal
import matplotlib.pyplot as plt
# import utils
from utilsProcess import *
from utilsIO import *

###################
### SOME GLOBAL PATHS
###################
testProjectPath = os.path.join("..","..","testProject")

###################
### TEST INTERPOLATION NUMBERS
###################
def testInterpValues():
	generalPrint("testsDataReader", "Running test function: testInterpValues")
	# set a dummy start time
	# startTime = datetime.strptime("2016-03-01 13:22:36.9334", "%Y-%m-%d %H:%M:%S.%f")
	# fs = 250.0
	startTime = datetime.strptime("2016-03-01 13:22:36.888", "%Y-%m-%d %H:%M:%S.%f")
	fs = 500.0	

	# create some dummy data
	chans = ["Ex", "Ey", "Hx", "Hy"]
	sinFreq = [1, 3, 2, 4]
	data = {}
	samples = 2001
	period = 1.0/fs
	for idx, chan in enumerate(chans):
		vals = np.arange(0, samples)*period*2*math.pi*sinFreq[idx]
		data[chan] =  np.sin(vals)
		
	startTimeInterp, samplesInterp, dataInterp = interpolateToSecond(fs, startTime, data)
	# create datetime for the initial data
	x = np.empty(shape=(samples), dtype=datetime)
	for i in xrange(0, samples):
		x[i] = startTime + timedelta(seconds=1.0*i/fs)

	generalPrint("testsDataReader", "Non interpolated: start time {}, end time {}".format(startTime.strftime("%Y-%m-%d %H:%M:%S.%f"), x[-1].strftime("%Y-%m-%d %H:%M:%S.%f")))
	# create datetime for the interpolated data
	xInterp = np.empty(shape=(samplesInterp), dtype=datetime)
	for i in xrange(0, samplesInterp):
		xInterp[i] = startTimeInterp + timedelta(seconds=1.0*i/fs)
	generalPrint("testsDataReader", "Non interpolated: start time {}, end time {}".format(startTimeInterp.strftime("%Y-%m-%d %H:%M:%S.%f"), xInterp[-1].strftime("%Y-%m-%d %H:%M:%S.%f")))

	# now plot
	fig = plt.figure()
	for idx, chan in enumerate(chans):
		plt.subplot(len(chans), 1, idx + 1)
		plt.plot(x, data[chan], "o-", markersize=8, label="Original")
		plt.plot(xInterp, dataInterp[chan], "*--", markersize=12, label="Interpolated")
		plt.grid()
		plt.legend()
		plt.title(chan)
	plt.show()

###################
### TEST ATS FILES AND RESAMPLING
###################
def testResampleAndWriter():
	generalPrint("testsDataReader", "Running test function: testResampleAndWriter")
	# ATS data file
	chans = ["Ex", "Ey", "Hx", "Hy"]
	dataPathATS = os.path.join(testProjectPath, "timeData", "M10", "meas_2016-03-02_12-10-21")
	dataPathInternal = os.path.join(testProjectPath, "timeData", "M10",  "testWriteUnscaled")
	atsReader = DataReaderATS(dataPathATS)
	atsReader.printInfo()
	# write out dataset as an example - write out Unscaled data
	writer = DataWriterInternal()
	writer.setOutPath(dataPathInternal)
	writer.printInfo()
	writer.writeDataset(atsReader)
	writer.printInfo()
	# now try reading in again
	internalReader = DataReaderInternal(dataPathInternal)
	internalReader.printInfo()

	# compare both datasets
	dataATSUnscaled = atsReader.getUnscaledSamples()
	dataInternalUnscaled = internalReader.getUnscaledSamples()

	# let's plot and compare
	ncols = 1
	nrows = 4
	fig = plt.figure(figsize=(20,nrows*4))
	for idx, c in enumerate(chans):
		plt.subplot(nrows, ncols, idx+1)
		plt.plot(dataATSUnscaled[c][:200], "o-", markersize=8, label="ATS")
		plt.plot(dataInternalUnscaled[c][:200], "^", markersize=8, label="Internal")
		plt.grid()
		plt.legend()
		plt.title(c)

	# now write a dataset out in LSB
	dataPathInternal = os.path.join(testProjectPath, "timeData", "M10",  "testWritePhysical")	
	writer.setOutPath(dataPathInternal)
	writer.printInfo()
	writer.writeDataset(atsReader, lsb_applied=True)
	writer.printInfo()

	# now compare the physical data
	dataATS = atsReader.getPhysicalSamples()
	dataInternal = internalReader.getPhysicalSamples()

	# let's plot and compare
	ncols = 1
	nrows = 4
	fig = plt.figure(figsize=(20,nrows*4))
	for idx, c in enumerate(chans):
		plt.subplot(nrows, ncols, idx+1)
		plt.plot(dataATS[c][:200], "o-", markersize=8, label="ATS")
		plt.plot(dataInternal[c][:200], "^", markersize=8, label="Internal")
		plt.grid()
		plt.legend()
		plt.title(c)

	# now let's resample dataATS
	# and then save that dataset
	dataPathResample = os.path.join(testProjectPath, "timeData", "M10", "testResample")
	fsOrig = atsReader.getSampleFreq()
	fsNew = 349
	dataResamp = resample(dataATS, fsOrig, fsNew)
	writer.setOutPath(dataPathResample)
	headers = atsReader.getHeaders()
	chanHeaders, chanMap = atsReader.getChanHeaders()
	writer.writeData(headers, chanHeaders, dataResamp, sample_freq=fsNew, lsb_applied=True)
	writer.printInfo()

	# read resampled data
	resampledReader = DataReaderInternal(dataPathResample)
	resampledReader.printInfo()
	dataResampLoad = resampledReader.getPhysicalSamples()
	numResampSamples = resampledReader.getNumSamples()
	x = np.arange(0, numResampSamples)
	x = x*fsOrig/fsNew

	# let's plot and compare
	ncols = 1
	nrows = 4
	fig = plt.figure(figsize=(20,nrows*4))
	for idx, c in enumerate(chans):
		plt.subplot(nrows, ncols, idx+1)
		plt.plot(dataInternal[c][:200], "o-", markersize=8, label="Internal")
		plt.plot(x[:200], dataResampLoad[c][:200], "^", markersize=8, label="Resampled load")
		plt.plot(x[:200], dataResamp[c][:200], ".--", markersize=8, label="Resample data")
		plt.grid()
		plt.legend()
		plt.title(c)
	plt.show()
 

###################
### TEST SPAM READERS AND INTERPOLATION
###################
def testSPAM():
	generalPrint("testsDataReader", "Running test function: testSPAM")
	# read in spam data
	dataPathSPAM = os.path.join(testProjectPath, "timeData", "Wittstock", "runComb")
	spamReader = DataReaderSPAM(dataPathSPAM)
	spamReader.printInfo()
	spamReader.printDataFileList()

	# lets try and get some data
	# data = spamReader.getUnscaledSamples(chans=["Hx", "Hy", "Ex", "Ey"], startSample=21599900, endSample=21605000)
	data = spamReader.getUnscaledSamples(chans=["Hx", "Hy", "Ex", "Ey"], startSample=2000, endSample=10000)	
	# data = spamReader.getUnscaledSamples(chans=["Hx", "Hy", "Ex", "Ey"])
	x = np.arange(0, spamReader.getNumSamples())	

	# let's try and interpolate on to the second
	startTimeInterp, numSamplesInterp, dataInterp = interpolateToSecond(spamReader.getSampleFreq(), spamReader.getStartDatetime(), data)
	generalPrint("testSPAM", "start time of interpolated data = {}".format(startTimeInterp))
	generalPrint("testSPAM", "number of samples of interpolated data = {}".format(numSamplesInterp))
	# calculate shifted array for plotting
	sampleShift = startTimeInterp - spamReader.getStartDatetime()
	shift = sampleShift.total_seconds()
	plotShift = shift/spamReader.getSampleRate()
	xInterp = np.arange(0, numSamplesInterp) + plotShift

	# now let's do a test write of this dataset
	dataPathInterp = os.path.join(testProjectPath,"timeData", "Wittstock", "runInterp")
	writer = DataWriterInternal()
	writer.setOutPath(dataPathInterp)
	headers = spamReader.getHeaders()
	chanHeaders, chanMap = spamReader.getChanHeaders()
	writer.writeData(headers, chanHeaders, dataInterp, start_time=startTimeInterp.strftime("%H:%M:%S.%f"), start_date=startTimeInterp.strftime("%Y-%m-%d"), num_samples=numSamplesInterp, lsb_applied=True)
	writer.printInfo()

	# now try and resample this dataset
	fsOrig = spamReader.getSampleFreq()
	fsResamp = 128.0
	dataResamp = resample(dataInterp, fsOrig, fsResamp)
	dataPathResamp = os.path.join(testProjectPath,"timeData", "Wittstock", "runResamp")	
	writer.setOutPath(dataPathResamp)
	writer.writeData(headers, chanHeaders, dataResamp, start_time=startTimeInterp.strftime("%H:%M:%S.%f"), start_date=startTimeInterp.strftime("%Y-%m-%d"), sample_freq=fsResamp, lsb_applied=True)
	writer.printInfo()
	# calculate the x axis for the resamp data
	xResamp = plotShift + np.arange(0, numSamplesInterp)*(fsOrig/fsResamp)

	# get filtered data for comparison to resampled data
	dataFilt = lpFilter(dataInterp, 256, 55)
	
	# get chans
	chans = sorted(list(data.keys()))
	# create figure	
	ncols = 1
	nrows = 4
	samplesToPlot = 500
	fig = plt.figure(figsize=(20,nrows*4))
	# calculate sampling frequency for each one
	for idx, c in enumerate(chans):
		plt.subplot(nrows, ncols, idx+1)
		# plt.plot(x[0:samplesToPlot], data[c][0:samplesToPlot], "o-", label="Original")
		# plt.plot(xInterp[0:samplesToPlot], dataInterp[c][0:samplesToPlot], "x-", label="Interpolated")
		plt.plot(xInterp[0:samplesToPlot], dataFilt[c][0:samplesToPlot], "x-", label="Filtered")
		plt.plot(xResamp[0:samplesToPlot/2], dataResamp[c][0:samplesToPlot/2], "s-", label="Resampled")	
		plt.legend()
		plt.grid()
		plt.title(c)
	plt.show()


###################
### RUN THE TESTING METHODS
###################
testInterpValues()
# testResampleAndWriter()
# testSPAM()