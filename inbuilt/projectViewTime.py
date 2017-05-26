#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('pythonMT_dev', 'core'))
sys.path.append(os.path.join('pythonMT_dev', 'utils'))
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
# import my classes
from project import Project
from calibrator import Calibrator
# utilities
from utilsProcess import *
from utilsIO import generalPrint
from utilsChecks import *
from utilsPlotter import getPlotFonts

def projectViewTime(proj, startDate, endDate, **kwargs):
	generalPrint("ProjectViewTime", "Showing time data with options: {}".format(kwargs))
	# get options
	options = parseKeywords(getDefaultOptions(proj), kwargs)
	# format startDate and endDate
	start = datetime.strptime("{}.000".format(startDate), "%Y-%m-%d %H:%M:%S.%f")
	end = datetime.strptime("{}.000".format(endDate), "%Y-%m-%d %H:%M:%S.%f")

	# create a calibrator in case
	cal = Calibrator(proj.getCalDataPath())
	if options["calibrate"]:
		cal.printInfo()			

	dataAll = {}
	xAll = {}
	fsAll = {}
	# first need to collect the relevant data
	# loop over sites
	for s in options["sites"]:
		# print site info
		# proj.printSiteInfo(s)

		# get measurement directories for site s
		timeMeas = proj.getSiteTimeFiles(s)

		# get the dictionary ready
		dataAll[s] = {}
		xAll[s] = {}
		fsAll[s] = {}
		# loop over measurement folders and calculate spectra for each one
		for meas in timeMeas:
			# get measurement sample frequency
			fs = proj.getMeasSampleFreq(s, meas)
			# check to see if in given frequency list
			if int(fs) not in options["freqs"]:
				continue

			# now find recordings that are within the recording time
			siteStart = proj.getMeasStartTime(s, meas)
			siteEnd = proj.getMeasEndTime(s, meas)
			if siteEnd < start or siteStart > end:
				continue

			# now get the data
			reader = proj.getMeasDataReader(s, meas)
			reader.printInfo()
			# get the samples of the datetimes
			sampleStart, sampleEnd = reader.time2sample(start, end)
			# and then go back and get the times - this protects against non second sampling
			# as the samples returned from time2sample are rounded 
			# using sample2time, we get the appropriate start and end times for those samples
			readStart, readEnd = reader.sample2time(sampleStart, sampleEnd)
			# get the data
			# data = reader.getPhysicalData(readStart, readEnd, chans=options["chans"])
			data = reader.getPhysicalData(readStart, readEnd, chans=options["chans"])
			# if calibration is on, calibrate the data
			if options["calibrate"]:
				generalPrint("ProjectViewTime", "Calibrating time data: site {}, time data {}".format(s, meas))
				sensors = reader.getSensors(reader.getChannels())
				serials = reader.getSerials(reader.getChannels())
				choppers = reader.getChoppers(reader.getChannels())	
				data = cal.calibrate(data, fs, sensors, serials, choppers)		

			# have the data, now want to calculate the x array
			samples = data[options["chans"][0]].size
			fsDelta = timedelta(seconds=1.0/fs)
			x = np.empty(shape=(samples), dtype=datetime)
			for i in xrange(0, samples):
				x[i] = readStart + timedelta(seconds=1.0*i/fs)
			# save the data
			dataAll[s][meas] = data
			xAll[s][meas] = x
			fsAll[s][meas] = fs

	# once all the data has been collected, plot it all
	fig = plt.figure(figsize=options["figsize"])
	plotFonts = options["plotfonts"]
	# suptitle
	st = fig.suptitle("Time data from {} to {}".format(start.strftime("%Y-%m-%d %H-%M-%S"), end.strftime("%Y-%m-%d %H-%M-%S")), fontsize=plotFonts["suptitle"])			
	st.set_y(0.98)	
	for s in dataAll:
		for meas in dataAll[s]:
			# apply the filter options
			if options["lpfilt"]:
				dataAll[s][meas] = lpFilter(dataAll[s][meas], fsAll[s][meas], options["lpfilt"])
			if options["hpfilt"]:
				dataAll[s][meas] = hpFilter(dataAll[s][meas], fsAll[s][meas], options["hpfilt"])
			if options["bpfilt"]:
				dataAll[s][meas] = hpFilter(dataAll[s][meas], fsAll[s][meas], options["bpfilt"][0], options["bpfilt"][1])

			generalPrint("ProjectViewTime", "Plotting {} - {} from {}  to {}".format(s, meas, xAll[s][meas][0], xAll[s][meas][-1]))	
			# now plot the data
			for idx, chan in enumerate(options["chans"]):
				plt.subplot(4, 1, idx+1)
				plotData = dataAll[s][meas][chan]
				if options["normalise"]: # then normalise the data
					plotData = plotData/np.linalg.norm(plotData)
				plt.plot(xAll[s][meas], plotData, label="{} - {}".format(s, meas))
	
	for idx, chan in enumerate(options["chans"]):
		ax = plt.subplot(4, 1, idx+1)	
		plt.title("Channel {}".format(chan), fontsize=plotFonts["title"])
		plt.grid()
		if idx == len(options["chans"])-1:
			plt.xlabel("Time", fontsize=plotFonts["axisLabel"])
		# limit the x-axis
		plt.xlim([start, end])
		# do the yaxis
		if isElectric(chan):
			plt.ylabel("mV/km", fontsize=plotFonts["axisLabel"])
			if len(options["Eylim"]) > 0:
				plt.ylim(options["Eylim"])
		else:
			if options["calibrate"]:
				plt.ylabel("nT", fontsize=plotFonts["axisLabel"]) 
			else:
				plt.ylabel("mV", fontsize=plotFonts["axisLabel"])	
			if len(options["Hylim"]) > 0:
				plt.ylim(options["Hylim"])			
		plt.legend(fontsize=plotFonts["legend"])	
		# set tick sizes
		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		    label.set_fontsize(plotFonts["axisTicks"])		
	
	plt.tight_layout()
	# shift subplots down, make room for suptitle				
	fig.subplots_adjust(top=0.92)	
	if options["save"]:
		fig.savefig(os.path.join(proj.getImageDataPath(), "timeData_{}_{}".format(start.strftime("%Y-%m-%d_%H-%M-%S_"), end.strftime("%Y-%m-%d_%H-%M-%S"))))	
	if options["show"]:
		plt.show()
	plt.close("all")
	return fig

def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()	
	default["chans"] = ["Ex", "Ey", "Hx", "Hy"]
	default["calibrate"] = False		
	default["normalise"] = False
	default["lpfilt"] = 0
	default["hpfilt"] = 0
	default["bpfilt"] = 0
	default["figsize"] = (30,15)
	default["show"] = True
	default["save"] = False
	default["Eylim"] = []
	default["Hylim"] = []
	default["plotfonts"] = getPlotFonts()
	# need to think about how to do masks well
	return default

def parseKeywords(default, keywords):
	# check user options
	for w in default:
		if w in keywords:
			default[w] = keywords[w]
	return default	 