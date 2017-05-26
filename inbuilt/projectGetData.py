#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('pythonMT_dev', 'core'))
sys.path.append(os.path.join('pythonMT_dev', 'utils'))
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import copy
# import my classes
from project import Project
from calibrator import Calibrator
from spectrumCalculator import SpectrumCalculator
from spectrumReader import SpectrumReader
# utilities
from utilsPlotter import *
from utilsProcess import *
from utilsIO import *

def projectGetTime(proj, site, meas, **kwargs):
	# print project information
	proj.printInfo()

	# print site info
	proj.printSiteInfo(site)

	# print measurement info
	proj.printMeasInfo(site, meas)

	# get measurement sample frequency
	fs = proj.getMeasSampleFreq(site, meas)
	# get measurement start and end times
	datetimeStart = proj.getMeasStartTime(site, meas)
	datetimeEnd = proj.getMeasEndTime(site, meas)
	# get the measurement reader and chans
	reader = proj.getMeasDataReader(site, meas)
	chans = reader.getChannels()

	# set the default parameters
	decimation = [1]
	calibrate = False
	notch = []
	detrend = True

	# apply options
	if "chans" in kwargs:
		chans = kwargs["chans"]
	if "start" in kwargs:
		datetimeStart = datetime.strptime(kwargs["start"], '%Y-%m-%d %H:%M:%S')
	if "stop" in kwargs:
		datetimeEnd = datetime.strptime(kwargs["stop"], '%Y-%m-%d %H:%M:%S')
	if "decimation" in kwargs:
		if type(kwargs["decimation"]) != type(list()):
			decimation = [kwargs["decimation"]]
		else:
			decimation = kwargs["decimation"]
	if "calibrate" in kwargs:
		calibrate = kwargs["calibrate"]	
	if "detrend" in kwargs:
		calibrate = kwargs["detrend"]			
	if "notch" in kwargs:
		notch = kwargs["notch"]			

	# get correct data start and end times
	dataStartTime, dataEndTime = reader.getDataTimes(datetimeStart, datetimeEnd)

	# get data, sensor info, serial info, chopper info for calibration
	data = reader.getPhysicalData(chans, dataStartTime, dataEndTime)
	sensors = reader.getSensors(chans)
	serials = reader.getSerials(chans)
	choppers = reader.getChoppers(chans)

	if calibrate:
		cal = Calibrator(proj.getCalDataPath())
		cal.printInfo()	
		# calibrate		
		data = cal.calibrate(data, fs, sensors, serials, choppers)

	# notch filter if required
	for n in notch:
		for c in data:
			data[c] = notchFilter(data[c], fs, n, n/5.0)		

	# lists for saving the data
	timeData = []
	timeX = []
	specData = []
	specX = []
	sampleFreqs = []

	# set fsDec
	fsDec = fs
	# now for each decimation, get data
	for d in decimation:
		# decimate
		data = downsampleTime(data, d, 51)
		fsDec = fsDec/d
		size = data[chans[0]].size
		# create spectrum calculator
		specCalc = SpectrumCalculator(fsDec, size)			

		# now calculate the spectrum
		f, fData = specCalc.calcFourierCoeff(data)

		# save data
		if detrend:
			for c in list(data.keys()):
				data[c] = signal.detrend(data[c], type="linear")	
		timeData.append(copy.deepcopy(data))
		timeX.append([dataStartTime, dataEndTime])
		specData.append(fData)
		specX.append(f)
		sampleFreqs.append(fsDec)	

	return timeX, timeData, specX, specData


# stack spectra and save images 
# good for QC
def projectGetSpec(proj, site, meas, iDec, windows, **kwargs):
	# print project information
	proj.printInfo()

	# print site info
	proj.printSiteInfo(site)

	# print measurement info
	proj.printMeasInfo(site, meas)

	# open the spec reader
	check = specReader.openBinaryForReading('spectra', iDec)
	if not check:
		generalPrint("ProjectGetSpec", "spectra file does not exist")
		return # probably because this decimation level not calculated
	specReader.printInfo()
	fsDec = specReader.getSampleFreq()
	dataChans = specReader.getChannels()	
	# calculate frequency array
	f = np.linspace(0, fsDec/2.0, specReader.getDataSize())				
	winData = specReader.readBinaryWindowLocal(iW)

	return f, specData



# decimation = [1,2,2,4,4,4]

# chans = ["Ex"]
# # fig = dataView(proj, site, meas, start=startTime, stop=stopTime, decimation=decimation, chans=chans, calibrate=True)
# fig = dataView(proj, site, meas, start=startTime, stop=stopTime, decimation=decimation, chans=chans)
# fig.tight_layout()
# fig.savefig("dataEx")



