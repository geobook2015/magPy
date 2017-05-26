#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('pythonMT_dev', 'core'))
sys.path.append(os.path.join('pythonMT_dev', 'stats'))
sys.path.append(os.path.join('pythonMT_dev', 'utils'))
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
# import my classes
from project import Project
from transferFunctionReader import TransferFunctionReader
# utilities
from utilsPlotter import *
from utilsIO import generalPrint
from utilsPlotter import getPlotFonts

def projectViewTF(proj, **kwargs):
	generalPrint("ProjectViewTF", "Showing project transfer functions with options: {}".format(kwargs))
	# get options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	
	
	# loop over sites
	for s in options["sites"]:
		# now plot TF for this site
		if options["oneplot"]:
			fig = projectPlotSingle(proj, s, **kwargs)
		else:
			fig = projectPlotMulti(proj, s, **kwargs)

	if options["show"]:
		plt.show()

# plot the polarisations on the same plot
def projectPlotSingle(proj, site, **kwargs):
	# get options
	options = parseKeywords(getDefaultOptions(proj), kwargs)	
	# print site info
	proj.printSiteInfo(site)

	# get the sample frequencies for the site
	sampleFreqs = set(proj.getSiteSampleFreqs(site))	
	# find the intersection with the options["freqs"]
	sampleFreqs = sampleFreqs.intersection(options["freqs"])
	# after the intersection, make this a sorted list again
	sampleFreqs = sorted(list(sampleFreqs))

	# now loop over the prepend options
	# if prepend is a string, then make it a list
	if isinstance(options["prepend"], basestring):
		options["prepend"] = [options["basestring"]]

	plotFonts = options["plotfonts"]	
	for prepend in options["prepend"]:
		# create figure		
		fig = plt.figure(figsize=options["figsize"])
		# sup title
		st = fig.suptitle("Site {}: {}: fs = {}".format(site, prepend, arrayToStringInt(sampleFreqs)), fontsize=plotFonts["suptitle"])			
		st.set_y(0.98)
		# plot			
		mks = ["o", "*", "d", "^", "h"]
		colours = {"ExHx":"orange", "EyHy":"green", "ExHy":"red", "EyHx":"blue"}
		# some plot options
		nrows = 2 # amplitude and phase
		ncols = 1
		dataFound = 0
		# loop over sampling frequencies
		for fs, mk in zip(sampleFreqs, mks):
			path = os.path.join(proj.getTransDataPathSite(site), "{:d}".format(int(fs)), "{}_fs{:d}_{}".format(site, int(fs), prepend))
			# check path - if does not exist, continue
			if not checkFilepath(path):
				continue
			dataFound += 1 # increment by 1 if something found	
			tfReader = TransferFunctionReader(path)
			periods = tfReader.getPeriods()
			res = {}
			phase = {}
			for pol in options["polarisations"]:
				res[pol], phase[pol] = tfReader.getResAndPhase(pol)

			# can do this in a loop
			for idx, pol in enumerate(options["polarisations"]):
				# amplitude
				plt.subplot(nrows, ncols, 1)
				plt.loglog(periods, res[pol], marker=mks[idx], markersize=9, markerfacecolor='none', markeredgecolor=colours[pol], mew=1.5,
					c=colours[pol], ls="dashed", label="{} - {}".format(fs, pol))
				# phase
				plt.subplot(nrows, ncols, 2)
				plt.semilogx(periods, phase[pol], marker=mks[idx], markersize=9, markerfacecolor='none', markeredgecolor=colours[pol], mew=1.5, 
					c=colours[pol], ls="dashed", label="{} - {}".format(fs, pol))

		# check if any files found
		if dataFound == 0:
			continue

		# now put the titles and such on
		# amplitude
		ax = plt.subplot(2, 1, 1)	
		plt.ylim(options["res_ylim"])
		plt.xlim(options["res_xlim"])	
		plt.grid(True)		
		plt.title("Apparent Resistivity", fontsize=plotFonts["title"])
		plt.xlabel("Period [s]", fontsize=plotFonts["axisLabel"])
		plt.ylabel("Apparent Resistivity [Ohm m]", fontsize=plotFonts["axisLabel"])
		plt.legend(loc="best", fontsize=plotFonts["legend"], framealpha=0.5)
		# make square
		ax.set_aspect('equal', adjustable='box')	
		# set tick sizes
		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		    label.set_fontsize(plotFonts["axisTicks"])							
		# phase
		ax = plt.subplot(2, 1, 2)
		plt.ylim(options["phase_ylim"])
		plt.xlim(options["phase_xlim"])
		plt.grid(True)		
		plt.title("Phase", fontsize=plotFonts["title"])
		plt.xlabel("Period [s]", fontsize=plotFonts["axisLabel"])
		plt.ylabel("Phase [degrees]", fontsize=plotFonts["axisLabel"])		
		plt.legend(loc="best", fontsize=plotFonts["legend"], framealpha=0.5)
		# set tick sizes
		for label in (ax.get_xticklabels() + ax.get_yticklabels()):
		    label.set_fontsize(plotFonts["axisTicks"])			

		fig.tight_layout()
		# shift subplots down, make room for suptitle				
		fig.subplots_adjust(top=0.92)		
		if options["save"]:
			fig.savefig(os.path.join(proj.getTransDataPathSite(site), "{}_{}.png".format(site, prepend)))	

	if not options["show"]:
		plt.close("all")	
	return fig

# plot the polarisations on multiple plots
def projectPlotMulti(proj, site, **kwargs):
	# get options
	options = parseKeywords(getDefaultOptions(proj), kwargs)		
	# print site info
	proj.printSiteInfo(site)

	# get the sample frequencies for the site
	sampleFreqs = set(proj.getSiteSampleFreqs(site))	
	# find the intersection with the options["freqs"]
	sampleFreqs = sampleFreqs.intersection(options["freqs"])
	# after the intersection, make this a sorted list again
	sampleFreqs = sorted(list(sampleFreqs))

	# now loop over the prepend options
	# if prepend is a string, then make it a list
	if isinstance(options["prepend"], basestring):
		options["prepend"] = [options["basestring"]]

	plotFonts = options["plotfonts"]	
	for prepend in options["prepend"]:
		# create figure
		fig = plt.figure(figsize=options["figsize"])
		# sup title
		st = fig.suptitle("Site {}: {}: fs = {}".format(site, prepend, arrayToStringInt(sampleFreqs)), fontsize=plotFonts["suptitle"])		
		st.set_y(0.98)
		# plot		
		mks = ["o", "*", "d", "^", "h"]
		# some plot info
		nrows = 2 # amplitude and phase
		ncols = len(options["polarisations"])
		dataFound = 0
		# loop over the sampling frequencies and plot	
		for fs, mk in zip(sampleFreqs, mks):
			path = os.path.join(proj.getTransDataPathSite(site), "{:d}".format(int(fs)), "{}_fs{:d}_{}".format(site, int(fs), prepend))
			# check path - if does not exist, continue			
			if not checkFilepath(path):
				continue
			dataFound += 1 # increment by 1 if something found				
			# plot the data
			tfReader = TransferFunctionReader(path)
			periods = tfReader.getPeriods()
			res = {}
			phase = {}
			for pol in options["polarisations"]:
				res[pol], phase[pol] = tfReader.getResAndPhase(pol)

			# can do this in a loop
			for idx, pol in enumerate(options["polarisations"]):
				# amplitude
				plt.subplot(nrows, ncols, idx + 1)
				plt.loglog(periods, res[pol], marker=mk, markersize=8, ls="dashed", label="{} - {}".format(fs, pol))
				# phase
				plt.subplot(nrows, ncols, ncols + idx + 1)
				plt.semilogx(periods, phase[pol], marker=mk, markersize=8, ls="dashed", label="{} - {}".format(fs, pol))

		# check if any files found
		if dataFound == 0:
			continue				

		# now put the titles and such on
		for idx, pol in enumerate(options["polarisations"]):	
			# amplitude
			ax = plt.subplot(nrows, ncols, idx + 1)	
			plt.ylim(options["res_ylim"])
			plt.xlim(options["res_xlim"])		
			plt.grid(True)		
			plt.title("{} Resistivity".format(pol), fontsize=plotFonts["title"])
			plt.xlabel("Period [s]", fontsize=plotFonts["axisLabel"])
			plt.ylabel("Apparent Resistivity [Ohm m]", fontsize=plotFonts["axisLabel"])
			plt.tick_params(axis='both', which='major', length=8, width=0.8, labelsize=plotFonts["axisTicks"])
			plt.legend(loc="lower right", fontsize=plotFonts["legend"], framealpha=0.5)					

			# phase
			ax = plt.subplot(nrows, ncols, ncols + idx + 1)
			plt.ylim(options["phase_ylim"])
			plt.xlim(options["phase_xlim"])
			plt.grid(True)		
			plt.title("{} Phase".format(pol), fontsize=plotFonts["title"])
			plt.xlabel("Period [s]", fontsize=plotFonts["axisLabel"])
			plt.ylabel("Phase [degrees]", fontsize=plotFonts["axisLabel"])	
			plt.tick_params(axis='both', which='major', length=8, width=0.8, labelsize=plotFonts["axisTicks"])	
			plt.legend(loc="lower right", fontsize=plotFonts["legend"], framealpha=0.5)

		fig.tight_layout()
		# shift subplots down, make room for suptitle				
		fig.subplots_adjust(top=0.92)		
		if options["save"]:
			fig.savefig(os.path.join(proj.getTransDataPathSite(site), "{}_{}.png".format(site, prepend)))	

	# return fig, which deals with showing
	if not options["show"]:
		plt.close("all")	
	return fig

def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()	
	default["polarisations"] = ["ExHx", "ExHy", "EyHx", "EyHy"]
	default["save"] = False
	default["show"] = True
	default["figsize"] = (8, 14)
	default["prepend"] = "transFunc"
	default["oneplot"] = True
	default["res_ylim"] = [0.000000001, 100000]
	default["res_xlim"] = [0.00001, 10000]
	default["phase_ylim"] = [-20,360]
	default["phase_xlim"] = [0.00001, 10000]
	default["plotfonts"] = getPlotFonts()
	return default

def parseKeywords(default, keywords):
	# parse the user supplied keywords
	for w in default:
		if w in keywords:
			default[w] = keywords[w]
	return default	 