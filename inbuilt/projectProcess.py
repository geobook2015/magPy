#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('pythonMT_dev', 'core'))
sys.path.append(os.path.join('pythonMT_dev', 'utils'))
import numpy as np
from datetime import datetime
# import my classes
from project import Project
from decimationParameters import DecimationParams
from windowParameters import WindowParams
from windowSelector import WindowSelector
from processorSingleSite import ProcessorSingleSite
from processorRemoteReference import ProcessorRemoteReference
# utils
from utilsIO import generalPrint

def projectProcess(proj, **kwargs):
	generalPrint("ProjectProcess", "Processing project data")
	# get options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	
	
	# get the reference time
	datetimeRef = proj.getRefTime()	

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
			projectProcessSingle(proj, fs, s, **kwargs)

# process a single sampling frequency for a single site
def projectProcessSingle(proj, fs, site, **kwargs):
	generalPrint("ProjectProcessSingle", "Processing site {}, sample frequency {}".format(site, fs))
	# get options
	options = parseKeywords(getDefaultOptionsSingle(proj, site), kwargs)	

	# define decimation parameters
	decParams = DecimationParams(fs)
	if len(options["evalfreq"]) == 0:
		decParams.setDecimationParams(options["declevels"], options["freqlevel"])
	else:
		decParams.setFrequencyParams(options["evalfreq"], options["declevels"], options["freqlevel"])

	# now do window parameters
	winParams = WindowParams(decParams)			

	# now do the window selector
	winSelector = WindowSelector(proj, fs, decParams, winParams)
	
	# if two sites are the same, winSelector only uses distinct sites
	# hence using site and inputSite is no problem even if they are the same
	if options["remotesite"]: # only remote site
		winSelector.setSites([site, options["inputsite"], options["remotesite"]])
	else: # single site
		winSelector.setSites([site, options["inputsite"]])

	# NEED a way to pass time information for datetime constraints
	# here, can set various constraints (datetime, statistic)
	# if int(fs) == 128 or int(fs) == 4096:
	# 	winSelector.addTimeConstraint("20:00:00", "06:00:00")				
	# print datetime constraints	
	
	# add window masks
	if len(options["masks"].keys()) > 0:
		for maskSite in options["masks"]:
			if isinstance(options["masks"][maskSite], basestring): # single mask
				winSelector.addWindowMask(maskSite, options["masks"][maskSite])
				continue
			if all(isinstance(item, basestring) for item in options["masks"][maskSite]): # list of masks for the site
				for m in options["masks"][maskSite]:
					winSelector.addWindowMask(maskSite, m)

	# # apply datetime constraint m1 128Hz
	# winSelector.addLevelTimeConstraint("19:00:00", "09:00:00", 0)
	# winSelector.addLevelTimeConstraint("19:00:00", "09:00:00", 1)
	# # apply for m1 65536 
	# winSelector.addDateConstraint("2016-03-31")

	# # apply datetime constraint m13 128Hz
	# winSelector.addLevelTimeConstraint("20:00:00", "08:00:00", 0)
	# winSelector.addLevelTimeConstraint("20:00:00", "08:00:00", 1)
	# winSelector.addLevelTimeConstraint("20:00:00", "08:00:00", 2)				
			
	# calculate the shared windows and print info
	winSelector.calcSharedWindows()
	winSelector.printInfo()
	winSelector.printDatetimeConstraints()
	winSelector.printWindowMasks()
	winSelector.printSharedIndices()
	winSelector.printWindowsForFrequency()
	
	# now have the windows, pass the winSelector to singleSiteProcessor
	if options["remotesite"]:
		generalPrint("projectProcessSingle", "Remote reference processing with sites: in = {}, out = {}, reference = {}".format(options["inputsite"], site, options["remotesite"]))
		processor = ProcessorRemoteReference(proj, winSelector)
		# add the remote site
		processor.setRemote(options["remotesite"], options["remotechans"])
	else:
		generalPrint("projectProcessSingle", "Single site processing with sites: in = {}, out = {}".format(options["inputsite"], site))		
		processor = ProcessorSingleSite(proj, winSelector)
	# add the input and output site	
	processor.setInput(options["inputsite"], options["inchans"])
	processor.setOutput(site, options["outchans"])
	processor.setPrepend(options["prepend"])	
	processor.printInfo()
	processor.process()	

def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()
	default["evalfreq"] = []		
	default["declevels"] = 7
	default["freqlevel"] = 6
	default["prepend"] = ""		
	default["inchans"] = ["Hx", "Hy"]
	default["outchans"] = ["Ex", "Ey"]
	default["remotechans"] = default["inchans"]
	default["remotesite"] = ""
	# need to think about how to do masks well
	return default

def getDefaultOptionsSingle(proj, site):
	default = getDefaultOptions(proj)
	default["inputsite"] = site
	default["masks"] = {}
	return default

def parseKeywords(default, keywords):
	# check user options
	for w in default:
		if w in keywords:
			default[w] = keywords[w]	
	return default	