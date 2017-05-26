#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('..', 'core'))
sys.path.append(os.path.join('..', 'stats'))
sys.path.append(os.path.join('..', 'utils'))
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
# import my classes
from project import Project
from decimationParameters import DecimationParams
from windowParameters import WindowParams
from windowSelector import WindowSelector
from statisticFrequency import StatisticFrequency
from plotter import PlotterStatistic
# utils
from utilsWindow import *
from utilsIO import *
from utilsStats import *
from utilsPlotter import *

# plot the various statistics by measurement file
def projectStatPlotMeas(proj, **kwargs):
	generalPrint("ProjectStatPlotMeas", "Plotting statistics for individual measurement files with options: {}".format(kwargs))
	# get options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	
	# get the reference time
	datetimeRef = proj.getRefTime()	
	# creat the statistic plotter
	plotter = PlotterStatistic()

	# loop over sites
	for s in options["sites"]:
		# print site info
		proj.printSiteInfo(s)

		# get measurement directories for site s
		timeMeas = proj.getSiteTimeFiles(s)

		# loop over and plot stats for each
		for meas in timeMeas:
			# get measurement sample frequency
			fs = proj.getMeasSampleFreq(s, meas)			
			# check to see if in given frequency list
			if fs not in options["freqs"]:
				continue	

			# output some text
			generalPrint("Statistic Plot", "Plotting stats for site {}, measurement {}".format(s, meas))

			# define decimation parameters
			decParams = DecimationParams(fs)
			if len(options["evalfreq"]) == 0:
				decParams.setDecimationParams(options["declevels"], options["freqlevel"])
			else:
				decParams.setFrequencyParams(options["evalfreq"], options["declevels"], options["freqlevel"])
			numLevels = decParams.getNumLevels()

			# now do window parameters - this helps calculate times for data
			winParams = WindowParams(decParams)	

			# now create the window selector - in case any masks are applied		
			winSelector = WindowSelector(proj, fs, decParams, winParams)
			winSelectorSites = [s]
			if options["remotesite"]:
				winSelectorSites.append(options["remotesite"])
			winSelector.setSites(winSelectorSites)
			# do the postpend for masks
			maskpend = ""			
			# add window masks - this might be different for each site
			if len(options["masks"].keys()) > 0:
				for maskSite in options["masks"]:
					if isinstance(options["masks"][maskSite], basestring): # single mask
						winSelector.addWindowMask(maskSite, options["masks"][maskSite])
						continue
					if all(isinstance(item, basestring) for item in options["masks"][maskSite]): # list of masks for the site
						for m in options["masks"][maskSite]:
							winSelector.addWindowMask(maskSite, m)	
			winSelector.calcSharedWindows()	

			# loop through decimation levels
			for iDec in xrange(0, numLevels):

				# sample frequency of the level
				fsDec = decParams.getSampleFreqLevel(iDec)

				# get the evaluation frequencies
				evalFreq = decParams.getEvalFrequenciesForLevel(iDec)				

				# create statistic readers
				statHandlers = {}
				statData = {}
				statWindows = [] 
				# create the statistic handlers
				for idx, stat in enumerate(options["stats"]):
					statName, statElements = statNameMap(stat)
					statHandlers[stat] = StatisticFrequency(statName)
					statData[stat] = []									
					# open the data file
					check = statHandlers[stat].readStatFile(proj.getStatDataPathMeas(s, meas), iDec)
					if check: # then exists, add data
						statData[stat].append(statHandlers[stat].getStats())	
						if idx == 0: # save window information for the statistics - only need to do this once	
							statWindows.append(statHandlers[stat].getGlobalIndices())

				# check for data, if no data, then continue
				check = True
				for sH in statData:
					check = check and len(statData[sH]) != 0
				if not check:
					continue

				# initialise properties of the plotter
				plotter.setDataInfo(datetimeRef, iDec, fsDec, winParams, statWindows)
				plotter.addWinSelector(winSelector)

				# now can plot the data
				for idx, eFreq in enumerate(evalFreq):
					# plots
					for sH in statHandlers:
						if sH in options["plots"]:
							plotname = os.path.join(proj.getStatDataPathMeas(s, meas), statHandlers[sH].getStatName(), statHandlers[sH].getStatName() + "_{:d}_{:d}{}".format(iDec, idx, maskpend))
							plotSingleStat(plotter, plotname, statHandlers[sH], statData[sH], idx, eFreq)
	
					if "transFuncCoherence" in options["plots"]:
						transFuncPlot_Coh = os.path.join(proj.getStatDataPathMeas(s, meas), statHandlers["transFunc"].getStatName(), "coh_" + statHandlers["transFunc"].getStatName() + "_{:d}_{:d}{}".format(iDec, idx, maskpend))
						plotTFCoherenceStat(plotter, transFuncPlot_Coh, statHandlers["transFunc"], statData["transFunc"], statData["coherence"], idx, eFreq)

					if "resPhaseCoherence" in options["plots"]:					
						resPhasePlot_Coh = os.path.join(proj.getStatDataPathMeas(s, meas), statHandlers["resPhase"].getStatName(), "coh_" + statHandlers["resPhase"].getStatName() + "_{:d}_{:d}{}".format(iDec, idx, maskpend))
						plotRPCoherenceStat(plotter, resPhasePlot_Coh, statHandlers["resPhase"], statData["resPhase"], statData["coherence"], idx, eFreq)						


# plot the various statistics by sampling frequency
def projectStatPlotFs(proj, **kwargs):
	generalPrint("ProjectStatPlotFs", "Plotting statistics by sample frequency with options: {}".format(kwargs))
	# get options
	options = parseKeywords(getDefaultOptions(proj), kwargs)
	# get the reference time
	datetimeRef = proj.getRefTime()	
	# creat the statistic plotter
	plotter = PlotterStatistic()	

	# loop over sites
	for s in options["sites"]:
		# print site info
		proj.printSiteInfo(s)

		# loop over sampling frequencies
		siteFs = proj.getSiteSampleFreqs(s)

		for fs in siteFs:
			# skip if not included
			if int(fs) not in options["freqs"]:
				continue

			# output some text
			generalPrint("Statistic Plot", "Plotting stats for site {}, sample frequency {:f}".format(s, int(fs)))

			# get measurement directories for site s
			timeMeas = proj.getSiteTimeFilesFs(s, fs)

			# define the plotting path
			plotPath = os.path.join(proj.getStatDataPathSite(s), "{:d}".format(int(fs)))
			checkAndMakeDir(plotPath)			

			# define decimation parameters
			decParams = DecimationParams(fs)
			if len(options["evalfreq"]) == 0:
				decParams.setDecimationParams(options["declevels"], options["freqlevel"])
			else:
				decParams.setFrequencyParams(options["evalfreq"], options["declevels"], options["freqlevel"])
			numLevels = decParams.getNumLevels()

			# now do window parameters - this helps calculate times for data
			winParams = WindowParams(decParams)

			# now create the window selector - in case any masks are applied		
			winSelector = WindowSelector(proj, fs, decParams, winParams)
			winSelectorSites = [s]
			if options["remotesite"]:
				winSelectorSites.append(options["remotesite"])
			winSelector.setSites(winSelectorSites)
			# do the postpend for masks
			maskpend = ""
			# add window masks
			if len(options["masks"].keys()) > 0:
				maskpend = "_withMask"
				for maskSite in options["masks"]:
					if isinstance(options["masks"][maskSite], basestring): # single mask
						winSelector.addWindowMask(maskSite, options["masks"][maskSite])
						continue
					if all(isinstance(item, basestring) for item in options["masks"][maskSite]): # list of masks for the site
						for m in options["masks"][maskSite]:
							winSelector.addWindowMask(maskSite, m)	
			
			# # apply datetime constraint m1
			# winSelector.addLevelTimeConstraint("19:00:00", "09:00:00", 0)
			# winSelector.addLevelTimeConstraint("19:00:00", "09:00:00", 1)
			winSelector.addDateConstraint("2016-03-31")			

			# # apply datetime constraint m13
			# winSelector.addLevelTimeConstraint("20:00:00", "08:00:00", 0)
			# winSelector.addLevelTimeConstraint("20:00:00", "08:00:00", 1)
			# winSelector.addLevelTimeConstraint("20:00:00", "08:00:00", 2)			

			winSelector.calcSharedWindows()														

			# loop through decimation levels
			for iDec in xrange(0, numLevels):
				# sample frequency of the level
				fsDec = decParams.getSampleFreqLevel(iDec)

				# get the evaluation frequencies
				evalFreq = decParams.getEvalFrequenciesForLevel(iDec)

				# create statistic readers
				statHandlers = {}
				statData = {}
				statWindows = [] # this records the number of windows				
				# create the statistic handlers
				for idx, stat in enumerate(options["stats"]):
					statName, statElements = statNameMap(stat)
					statHandlers[stat] = StatisticFrequency(statName)
					statData[stat] = []	
					# loop over time files			
					for meas in timeMeas:								
						# open the data file
						check = statHandlers[stat].readStatFile(proj.getStatDataPathMeas(s, meas), iDec)
						if check: # then exists, add data
							statData[stat].append(statHandlers[stat].getStats())	
							if idx == 0: # save window information for the statistics - only need to do this once
								statWindows.append(statHandlers[stat].getGlobalIndices())						
				
				# check for data, if no data, then quit
				check = True
				for sH in statData:
					check = check and len(statData[sH]) != 0
				if not check:
					continue

				# initialise properties of the plotter
				plotter.setDataInfo(datetimeRef, iDec, fsDec, winParams, statWindows)
				plotter.addWinSelector(winSelector)

				# now can plot the data
				for idx, eFreq in enumerate(evalFreq):
					# plots
					for sH in statHandlers:
						if sH in options["plots"]:
							plotname = os.path.join(plotPath, statHandlers[sH].getStatName() + "_{:d}_{:d}{}".format(iDec, idx, maskpend))
							plotSingleStat(plotter, plotname, statHandlers[sH], statData[sH], idx, eFreq)
	
					if "transFuncCoherence" in options["plots"]:
						transFuncPlot_Coh = os.path.join(plotPath, "coh_" + statHandlers["transFunc"].getStatName() + "_{:d}_{:d}{}".format(iDec, idx, maskpend))
						plotTFCoherenceStat(plotter, transFuncPlot_Coh, statHandlers["transFunc"], statData["transFunc"], statData["coherence"], idx, eFreq)

					if "resPhaseCoherence" in options["plots"]:					
						resPhasePlot_Coh = os.path.join(plotPath, "coh_" + statHandlers["resPhase"].getStatName() + "_{:d}_{:d}{}".format(iDec, idx, maskpend))
						plotRPCoherenceStat(plotter, resPhasePlot_Coh, statHandlers["resPhase"], statData["resPhase"], statData["coherence"], idx, eFreq)					

# helper plot functions
def plotSingleStat(plotter, plotname, statistic, data, eIdx, eFreq):
	name = statistic.getStatName()
	if "psdStat" in name:
		plotPSDStat(plotter, plotname, statistic, data, eIdx, eFreq)		
	if "cohStat" in name:
		plotCohStat(plotter, plotname, statistic, data, eIdx, eFreq)
	if "polStat" in name:
		plotPolStat(plotter, plotname, statistic, data, eIdx, eFreq)
	if "pcohStat" in name:
		plotPCohStat(plotter, plotname, statistic, data, eIdx, eFreq)		
	if "tfStat" in name:
		plotTFStat(plotter, plotname, statistic, data, eIdx, eFreq)
	if "resPhaseStat" in name:
		plotRPStat(plotter, plotname, statistic, data, eIdx, eFreq)	
	if "absvalEqn" in name:
		plotCrossPlotStat(plotter, plotname, statistic, data, eIdx, eFreq)

	plotter.close()

# individual plot functions
def getNoTime(data):
	# return False
	notime = False
	if len(data) > 1:
		notime = True
	return notime	

def plotPSDStat(plotter, plotname, psdStatistic, psdData, eIdx, eFreq):
	numWindows = calcNumWindows(psdData)
	fig = plotter.plotData1d(psdStatistic.getWindowStats(), psdData, 
		eIndex=eIdx, title="Power Spectral Density, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		xlabel="Time", ylabel="PSD", logy=True, ylim=[0.0000001, 1000], colsPerRow=2, notime=getNoTime(psdData)
	)
	fig.savefig(plotname)

def plotCohStat(plotter, plotname, cohStatistic, cohData, eIdx, eFreq):
	numWindows = calcNumWindows(cohData)
	fig = plotter.plotData1d(cohStatistic.getWindowStats(), cohData, 
		eIndex=eIdx, title="Coherence, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		# xlabel="Time", ylabel="Coherency", colsPerRow=4, notime=getNoTime(cohData)
		ylim=[-0.1, 1.1], xlabel="Time", ylabel="Coherency", colsPerRow=4, notime=getNoTime(cohData), hist=True
	)
	fig.savefig(plotname)

def plotPolStat(plotter, plotname, polStatistic, polData, eIdx, eFreq):
	numWindows = calcNumWindows(polData)
	fig = plotter.plotData1d(polStatistic.getWindowStats(), polData, 
		eIndex=eIdx, title="Polarisation Direction, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows),
		ylim=[-100, 100], xlabel="Time", ylabel="Polarisation Direction", colsPerRow=2, notime=getNoTime(polData)
	)
	fig.savefig(plotname)

def plotPCohStat(plotter, plotname, pcohStatistic, pcohData, eIdx, eFreq):
	numWindows = calcNumWindows(pcohData)
	fig = plotter.plotData1d(pcohStatistic.getWindowStats(), pcohData, 
		eIndex=eIdx, title="Partial and Bivariate Coherence, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		ylim=[-0.1, 1.1], xlabel="Time", ylabel="Coherence", colsPerRow=2, notime=getNoTime(pcohData)
	)
	fig.savefig(plotname)

def plotCrossPlotStat(plotter, plotname, statistic, statData, eIdx, eFreq):
	numWindows = calcNumWindows(statData)
	fig = plotter.plotData2d(statistic.getWindowStats(), statData, 
		eIndex=eIdx, title="Cross plot, input and output {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		xlabel="Input", ylabel="Output", colsPerRow=4, notime=getNoTime(statData), gradients=True
	)
	fig.savefig(plotname)	

def plotTFStat(plotter, plotname, transFuncStatistic, transFuncData, eIdx, eFreq):
	numWindows = calcNumWindows(transFuncData)
	maxVal = np.sqrt(1000*eFreq*5)	
	fig = plotter.plotData2d(transFuncStatistic.getWindowStats(), transFuncData, 
		eIndex=eIdx, title="Impedance tensor components, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		colsPerRow=2, xlim=[-maxVal, maxVal], ylim=[-maxVal, maxVal], xlabel="Impedance tensor real", ylabel="Impedance tensor imaginary"
	)
	fig.savefig(plotname)
	# also plot the vs time for each individual component	
	fig = plotter.plotData1d_2comp(transFuncStatistic.getWindowStats(), transFuncData, 
		eIndex=eIdx, title="Impedance tensor components, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		colsPerRow=4, ylim=[-maxVal, maxVal], ylimalt=[-maxVal,maxVal], xlabel="Time", ylabel="Value", notime=getNoTime(transFuncData), hist=True
	)
	fig.savefig(plotname + "_2comp")	

def plotTFCoherenceStat(plotter, plotname, transFuncStatistic, transFuncData, altData, eIdx, eFreq):
	numWindows = calcNumWindows(transFuncData)
	maxVal = np.sqrt(1000*eFreq*5)				
	fig = plotter.plotData2d_1d(transFuncStatistic.getWindowStats(), transFuncData, altData, 
		eIndex=eIdx, title="Impedance tensor components, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		colsPerRow=2, xlim=[-maxVal, maxVal], ylim=[-maxVal, maxVal], xlabel="Impedance tensor real", ylabel="Impedance tensor imaginary", cTitle="Coherence"
	)
	fig.savefig(plotname)

def plotRPStat(plotter, plotname, resPhaseStatistic, resPhaseData, eIdx, eFreq):
	numWindows = calcNumWindows(resPhaseData)		
	fig = plotter.plotData2d(resPhaseStatistic.getWindowStats(), resPhaseData, 
		eIndex=eIdx, title="App. Resistivity and Phase, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		colsPerRow=2, logx=True, xlim=[0.01, 1000000], ylim=[-180,180], xlabel="Resistivity [Ohm m]", ylabel="Phase [degrees]"
	)
	fig.savefig(plotname)
	# also plot the vs time for each individual component	
	fig = plotter.plotData1d_2comp(resPhaseStatistic.getWindowStats(), resPhaseData, 
		eIndex=eIdx, title="App. Resistivity and Phase, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		colsPerRow=4, logy=True, ylim=[0.0001,10000], logyalt=False, ylimalt=[-270, 270], xlabel="Time", ylabel="Value", notime=getNoTime(resPhaseData), hist=True
	)
	fig.savefig(plotname + "_2comp")				

def plotRPCoherenceStat(plotter, plotname, resPhaseStatistic, resPhaseData, altData, eIdx, eFreq):
	numWindows = calcNumWindows(resPhaseData)		
 	fig = plotter.plotData2d_1d(resPhaseStatistic.getWindowStats(), resPhaseData, altData,
		eIndex=eIdx, title="App. Resistivity and Phase, period {:6f} [s], frequency {:6f} [Hz], {:d} windows".format(1.0/eFreq, eFreq, numWindows), 
		colsPerRow=2, logx=True, xlim=[0.01, 1000000], ylim=[-180,180], xlabel="Resistivity [Ohm m]", ylabel="Phase [degrees]", cTitle="Coherence"
	)
	fig.savefig(plotname)

def calcNumWindows(data):
	numWindows = 0
	for d in data:
		# first index is global index, second is the data
		numWindows = numWindows + d.shape[0]
	return numWindows

def getStatName(stat, mode):
	statNames = {"coherence": "cohStat", "transFunc": "tfStat", "resPhase": "resPhaseStat"}
	if stat not in statNames:
		warningPrint("Statistic Plot", "Statistic not found. Exiting...")
		exit()
	name = statNames[stat]
	if mode == "remote":
		name = "RR_{}".format(name)
	return name

def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()
	default["evalfreq"] = []	
	default["declevels"] = 7
	default["freqlevel"] = 6
	default["remotesite"] = ""
	default["masks"] = {}	
	default["stats"] = ["absvalEqn", "coherence", "psd", "poldir", "transFunc", "resPhase", "partialcoh"]	
	default["plots"] = ["absvalEqn", "coherence", "psd", "poldir", "partialcoh", "transFunc", "transFuncCoherence", "resPhase", "resPhaseCoherence"]
	return default

def parseKeywords(default, keywords):
	# check keyword arguments
	for w in default:
		if w in keywords:
			default[w] = keywords[w]			
	return default	