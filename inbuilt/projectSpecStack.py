#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('..', 'core'))
sys.path.append(os.path.join('..', 'utils'))
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
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


# stack spectra and save images 
# good for QC
def projectSpecStack(proj, **kwargs):
	generalPrint("ProjectSpecStack", "Stacking spectra with options: {}".format(kwargs))
	# default options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	

	# calculate spectra stacks
	# get the reference time
	datetimeRef = proj.getRefTime()	

	# loop over sites
	for s in options["sites"]:
		# print site info
		proj.printSiteInfo(s)

		# get measurement directories for site s
		timeMeas = proj.getSiteTimeFiles(s)

		# loop over measurement folders and calculate spectra for each one
		for meas in timeMeas:
			# get measurement sample frequency
			fs = proj.getMeasSampleFreq(s, meas)			
			# check to see if in given frequency list
			if int(fs) not in options["freqs"]:
				continue

			# print measurement info
			proj.printMeasInfo(s, meas)		

			# define decimation parameters
			decParams = DecimationParams(fs)
			decParams.setDecimationParams(options["declevels"], options["freqlevel"])
			#decParams.printInfo()
			numLevels = decParams.getNumLevels()

			# create the spectrum reader
			specReader = SpectrumReader(proj.getSpecDataPathMeas(s, meas))
			specReader.printInfo()

			# loop through decimation levels
			for iDec in xrange(0, numLevels):	
				# open the spectra file for the current decimation level
				check = specReader.openBinaryForReading("spectra", iDec)
				if not check:
					continue # probably because this decimation level not calculated
				specReader.printInfo()

				# get the channels
				dataChans = specReader.getChannels()
				if len(options["chans"]) > 0:
					dataChans = options["chans"]

				# get windows
				numWindows = specReader.getNumWindows()
				fsDec = specReader.getSampleFreq()
				
				# calculate frequency array
				f = np.linspace(0, fsDec/2.0, specReader.getDataSize())				

				# calculate num of windows to stack per plot
				nStacks = options["numstacks"][iDec] 
				# round up and then have a smaller stack at the end
				stackSize = int(np.floor(1.0*numWindows/nStacks))

				for iP in xrange(0, nStacks):
					stackStart = iP*stackSize
					stackStop = min(stackStart + stackSize, numWindows)
					# array to hold the stacked data
					stackedData = {}
					ampData = {}
					phaseData = {}
					powerData = {}
					# assign
					for c in dataChans:
						stackedData[c] = np.zeros(shape=(specReader.getDataSize()), dtype="complex")
						ampData[c] = np.zeros(shape=(specReader.getDataSize()), dtype="complex")
						phaseData[c] = np.zeros(shape=(specReader.getDataSize()), dtype="complex")												
						for c2 in dataChans:						
							powerData[c+c2] = np.zeros(shape=(specReader.getDataSize()), dtype="complex")
					
					# now stack the data and create nice plots					
					for iW in xrange(stackStart, stackStop):
						winSpec = specReader.readBinaryWindowLocal(iW)
						for c in dataChans:
							stackedData[c] += winSpec[c]
							ampData[c] += np.absolute(winSpec[c])
							phaseData[c] += np.angle(winSpec[c])*(180.0/np.pi)
							# get coherency data
							for c2 in dataChans:
								powerData[c+c2] += winSpec[c]*np.conjugate(winSpec[c2])
							
					# scale powers and stacks
					for c in dataChans:
						stackedData[c] = stackedData[c]/(stackStop-stackStart)
						ampData[c] = ampData[c]/(stackStop-stackStart)
						phaseData[c] = phaseData[c]/(stackStop-stackStart)						
						for c2 in dataChans:
							powerData[c+c2] = 2*powerData[c+c2]/(stackStop-stackStart) # some normalisation		
							powerData[c+c2][[0,-1]] = powerData[c+c2][[0,-1]]/2	# some normalisation					

					# get date time information for windows
					startG = stackStart + specReader.getGlobalOffset()
					stopG = stackStop + specReader.getGlobalOffset()
					startStart, startEnd = gIndex2datetime(startG, datetimeRef, fsDec, specReader.getWindowSize(), specReader.getWindowOverlap())
					stopStart, stopEnd = gIndex2datetime(stopG, datetimeRef, fsDec, specReader.getWindowSize(), specReader.getWindowOverlap())
					
					# now plot and save
					# calculate number of rows
					# calculate number of rows - in case interested in coherencies too
					nrows = 2
					if len(options["coherencies"]) > 0:
						nrows = nrows + np.ceil(1.0*len(options["coherencies"])/len(dataChans))

					# create figure	
					plotFonts = options["plotfonts"]					
					fig = plt.figure(figsize=(30,5*nrows))
					st = fig.suptitle("Spectra stack, fs = {:.6f} [Hz], decimation level = {:2d}, windows = {:d}, {} to {}".format(
						fsDec, iDec, stackStop-stackStart, startStart.strftime("%d %B %Y, %H:%M:%S.%f"), stopEnd.strftime("%d %B %Y, %H:%M:%S.%f")
						), 
						fontsize=plotFonts["suptitle"]
					)			
					st.set_y(0.98)
					for idx, c in enumerate(dataChans):
						ax1 = plt.subplot(nrows, len(dataChans), idx+1)
						plt.title("Amplitude {}".format(c), fontsize=plotFonts["title"])
						ax1.semilogy(f, ampData[c])
						ax1.set_ylim(0.001, 100000)
						ax1.set_xlim(0, fsDec/2)
						if isMagnetic(c):
							ax1.set_ylabel("Amplitude [nT]", fontsize=plotFonts["axisLabel"])
						else:
							ax1.set_ylabel("Amplitude [mV/km]", fontsize=plotFonts["axisLabel"])
						ax1.set_xlabel("Frequency [Hz]", fontsize=plotFonts["axisLabel"])							
						plt.grid(True)
						# plot phase
						ax2 = plt.subplot(nrows, len(dataChans), len(dataChans)+idx+1)
						plt.title("Phase {}".format(c), fontsize=plotFonts["title"])						
						ax2.plot(f, phaseData[c])
						ax2.set_ylim(-190, 190)
						ax2.set_xlim(0, fsDec/2)
						ax2.set_ylabel("Phase [degrees]", fontsize=plotFonts["axisLabel"])
						ax2.set_xlabel("Frequency [Hz]", fontsize=plotFonts["axisLabel"])
						plt.grid(True)
					
					# plot coherencies
					for idx, coh in enumerate(options["coherencies"]):
						c = coh[0]
						c2 = coh[1]
						cohNom = np.power(np.absolute(powerData[c+c2]), 2)
						cohDenom = powerData[c+c] * powerData[c2+c2]
						coherence = cohNom/cohDenom 
						ax = plt.subplot(nrows, len(dataChans), 2*len(dataChans)+idx+1)
						plt.title("Coherence {} - {}".format(c, c2), fontsize=plotFonts["title"])
						ax.plot(f, coherence)
						ax.set_ylim(0, 1.1)
						ax.set_xlim(0, fsDec/2)
						ax.set_ylabel("Coherence", fontsize=plotFonts["axisLabel"])
						ax.set_xlabel("Frequency [Hz]", fontsize=plotFonts["axisLabel"])						
						plt.grid(True)

					# save plot
					filename = os.path.join(proj.getSpecDataPathMeas(s, meas), "{}_dec{:02d}_stackNum{:03d}".format(options["prepend"], iDec, iP))
					fig.tight_layout()
					# shift subplots down:					
					fig.subplots_adjust(top=0.92)
					fig.savefig(filename)
					plt.close("all")


def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()	
	default["chans"] = []
	default["declevels"] = 7
	default["freqlevel"] = 6
	default["numstacks"] = np.ones(shape=(default["declevels"]), dtype=int)	
	default["prepend"] = "specStack"		
	default["plotfonts"] = getPlotFonts()
	default["coherencies"] = []
	return default

def parseKeywords(default, keywords):
	# check keyword arguments
	words = ["sites", "freqs", "chans", "declevels", "freqlevel", "numstacks", "prepend", "plotfonts", "coherencies"]
	for w in words:
		if w in keywords:
			default[w] = keywords[w]	
	# do the except for numstacks
	if "numstacks" in keywords:
		default["numstacks"] = np.array(default["numstacks"], dtype=int)
	return default	