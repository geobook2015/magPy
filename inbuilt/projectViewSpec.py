#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('..', 'core'))
sys.path.append(os.path.join('..', 'utils'))
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
# import my classes
from project import Project
from decimationParameters import DecimationParams
from windowParameters import WindowParams
from spectrumReader import SpectrumReader
# utils
from utilsWindow import *
from utilsIO import *
from utilsChecks import *
from utilsPlotter import *

# view a single spectra
def projectViewSpec(proj, site, meas, **kwargs):
	generalPrint("projectViewSpecSingle", "Plotting spectra for measumrent {} and site {} with options: {}".format(meas, site, kwargs))
	# default options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	

	if meas not in proj.getSiteTimeFiles(site):
		warningPrint("projectViewSpecSingle", "Measurement directory {} for site {} not found".format(meas, site))
		exit()

	# print measurement info
	proj.printMeasInfo(site, meas)	
	fs = proj.getMeasSampleFreq(site, meas)

	# define decimation parameters
	decParams = DecimationParams(fs)
	decParams.setDecimationParams(options["declevels"], options["freqlevel"])
	#decParams.printInfo()
	numLevels = decParams.getNumLevels()

	# create the spectrum reader
	specReader = SpectrumReader(proj.getSpecDataPathMeas(site, meas))
	specReader.printInfo()

	# open the spectra file for the current decimation level
	check = specReader.openBinaryForReading("spectra", options["plotlevel"])
	if not check:
		# probably because this decimation level not calculated
		warningPrint("projectViewSpecSingle", "Spectra file does not exist at level {}".format(options["plotlevel"])) 
	specReader.printInfo()

	# get the channels
	dataChans = specReader.getChannels()
	if len(options["chans"]) > 0:
		dataChans = options["chans"]

	# get windows
	numWindows = specReader.getNumWindows()
	fsDec = specReader.getSampleFreq()			
		
	# get the window data
	windows = options["plotwindow"]
	if isinstance(windows, basestring) and windows == "all":
		if numWindows > 400:
			windows = list(np.linspace(0, numWindows, 400, endpoint=False, dtype=np.int32))
		else:
			windows = list(np.arange(0,numWindows))
	elif isinstance(windows, int):
		windows = [windows] # if an integer, make it into a list
	elif isinstance(windows, dict):
		windows = list(np.arange(dict["start"], dict["stop"]+1))

	# collect the data
	winSpec = np.empty(shape=(len(windows), len(dataChans), specReader.getDataSize()), dtype="complex")
	winMap = {}
	# get the window data
	count = 0
	for iW in windows:
		winMap[count] = iW
		for idx, chan in enumerate(dataChans):
			winSpec[count, idx] = specReader.readBinaryWindowLocal(iW)[chan]
		count += 1

	if options["section"]:
		fig = projectViewSpecSection(proj, site, meas, specReader, dataChans, winSpec, winMap, **kwargs)
	else:
		fig = projectViewSpecSingle(proj, site, meas, specReader, dataChans, winSpec, winMap, **kwargs)

	if options["save"]:
		fig.savefig(os.path.join(proj.getImageDataPath(), "{}_{}_dec{}".format(options["prepend"], site, options["plotlevel"])))
	if options["show"]:
		plt.show()
	plt.close("all")	


# view individual plots
def projectViewSpecSingle(proj, site, meas, specReader, dataChans, specData, winMap, **kwargs):
	# default options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	
	# freq array
	f = specReader.getFrequencyArray()	
	# create figure	
	plotFonts = options["plotfonts"]					
	fig = plt.figure(figsize=(20,3*len(dataChans)))
	st = fig.suptitle("Spectra plot, site = {}, meas = {}, fs = {:.2f} [Hz], decimation level = {:2d}".format(
			site, meas, specReader.getSampleFreq(), options["plotlevel"]
		), 
		fontsize=plotFonts["suptitle"]
	)			
	st.set_y(0.98)
	nrows = len(dataChans)
	for iWin, winData in enumerate(specData): 
		for idx, chan in enumerate(dataChans):
			# calculate start time of spec section
			globalI = winMap[iWin] + specReader.getGlobalOffset()
			startStart, startEnd = gIndex2datetime(globalI, proj.getRefTime(), specReader.getSampleFreq(), 
				specReader.getWindowSize(), specReader.getWindowOverlap()
			)			
			ax = plt.subplot(nrows, 1, idx+1)
			ax.semilogy(f, np.absolute(winData[idx]), label="global window {}, time {}".format(globalI, startStart))			

	# put on axis labels etc
	for idx, chan in enumerate(dataChans):
		ax = plt.subplot(nrows, 1, idx+1)		
		plt.title("Amplitude {}".format(chan), fontsize=plotFonts["title"])
		ax.set_ylim(options["amplim"])
		ax.set_xlim(0, specReader.getSampleFreq()/2.0)
		if isMagnetic(chan):
			ax.set_ylabel("Amplitude [nT]", fontsize=plotFonts["axisLabel"])
		else:
			ax.set_ylabel("Amplitude [mV/km]", fontsize=plotFonts["axisLabel"])
		if idx == len(dataChans)-1:
			ax.set_xlabel("Frequency [Hz]", fontsize=plotFonts["axisLabel"])
		# set tick sizes
		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		    label.set_fontsize(plotFonts["axisTicks"])											
		plt.grid(True)
		plt.legend()
	
	fig.tight_layout()
	# shift subplots down:					
	fig.subplots_adjust(top=0.92)
	return fig

# view a spectra section
def projectViewSpecSection(proj, site, meas, specReader, dataChans, specData, winMap, **kwargs):
	# default options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	
	# freq array
	f = specReader.getFrequencyArray()	
	# create figure	
	plotFonts = options["plotfonts"]					
	fig = plt.figure(figsize=(9*len(dataChans),15))
	st = fig.suptitle("Spectra sections, site = {}, meas = {}, fs = {:.2f} [Hz], decimation level = {:2d}, windows = {:d}, {} to {}".format(
			site, meas, specReader.getSampleFreq(), options["plotlevel"], len(specData), winMap[0], winMap[len(specData)-1]
		), 
		fontsize=plotFonts["suptitle"]
	)			
	st.set_y(0.98)
	# calculate out the dates
	dates = []
	for iWin, winData in enumerate(specData):
		globalI = winMap[iWin] + specReader.getGlobalOffset()
		startStart, startEnd = gIndex2datetime(globalI, proj.getRefTime(), specReader.getSampleFreq(), 
			specReader.getWindowSize(), specReader.getWindowOverlap()
		)	
		dates.append(startStart)
			
	for idx, chan in enumerate(dataChans):
		ax = plt.subplot(1, len(dataChans), idx+1)
		plotData = np.transpose(np.absolute(np.squeeze(specData[:,idx,:])))
		plt.pcolor(dates, f, plotData, norm=LogNorm(vmin=plotData.min(), vmax=plotData.max()), cmap=cm.magma)
		plt.colorbar()
		# amplim
		ax.set_ylim(0, specReader.getSampleFreq()/2.0)
		ax.set_xlim([dates[0], dates[-1]])
		if isMagnetic(chan):
			plt.title("Amplitude {} [nT]".format(chan), fontsize=plotFonts["title"])
		else:
			plt.title("Amplitude {} [mV/km]".format(chan), fontsize=plotFonts["title"])
		ax.set_ylabel("Frequency [Hz]", fontsize=plotFonts["axisLabel"])
		ax.set_xlabel("Time", fontsize=plotFonts["axisLabel"])
		# set tick sizes
		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		    label.set_fontsize(plotFonts["axisTicks"])								
		plt.grid(True)

	fig.autofmt_xdate()
	fig.tight_layout()
	# shift subplots down:					
	fig.subplots_adjust(top=0.92)
	return fig


def getDefaultOptions(proj):
	# default options
	default = {}
	default["freqs"] = proj.getAllSampleFreq()	
	default["chans"] = []
	default["declevels"] = 7
	default["freqlevel"] = 6
	default["plotlevel"] = 0
	default["plotwindow"] = [0]
	default["plotfonts"] = getPlotFonts()	
	default["section"] = False
	default["amplim"] = [0.001, 10000]	
	default["show"] = True
	default["save"] = False
	default["prepend"] = "viewSpec"		
	return default

def parseKeywords(default, keywords):
	# check keyword arguments
	for w in default:
		if w in keywords:
			default[w] = keywords[w]	
	# do the except for numstacks
	if "numstacks" in keywords:
		default["numstacks"] = np.array(default["numstacks"], dtype=int)
	return default	