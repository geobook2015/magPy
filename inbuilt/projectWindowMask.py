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
from windowMasker import WindowMasker
from windowSelector import WindowSelector
# from windowStatistic import WindowStatistic
# utils
from utilsWindow import *
from utilsIO import *
from utilsStats import *

# plot the various statistics by measurement file
def projectWindowMask(proj, **kwargs):
	generalPrint("ProjectWindowMask", "Creating window masks for project with options: {}".format(kwargs))
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

			# define decimation parameters
			decParams = DecimationParams(fs)
			if len(options["evalfreq"]) == 0:
				decParams.setDecimationParams(options["declevels"], options["freqlevel"])
			else:
				decParams.setFrequencyParams(options["evalfreq"], options["declevels"], options["freqlevel"])

			# now do window parameters - this helps calculate times for data
			winParams = WindowParams(decParams)

			# now create the window masker
			winMasker = WindowMasker(proj, s, fs, decParams, winParams)
			winMasker.addStats(options["stats"])
			winMasker.printInfo()
			# add constraints		


			### M1
			# #################
			# # 65536Hz masking
			# ################
			# before 2,0 looks good
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.3, 0.7], "cohEyHx": [0.3, 0.7]}, 1, 0)
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.3, 0.7], "cohEyHx": [0.3, 0.7]}, 1, 1)
			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [1500000, 2500000]}, 1, 2)		
			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [5e6, 2e7]}, 1, 3)		
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.6, 1.0], "cohEyHx": [0.6, 1.0]}, 2, 1)

			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [1000, 200000], "absEyEx": [1000, 200000], "absExEy": [1000, 200000], "absEyEy": [1000, 200000]}, 2, 2)	
			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [1000, 200000], "absEyEx": [1000, 200000], "absExEy": [1000, 200000], "absEyEy": [1000, 200000]}, 2, 3)
			
			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [50, 2500], "absEyEx": [50, 2500], "absExEy": [50, 2500], "absEyEy": [50, 2500]}, 3, 0)
			# winMasker.addConstraintFreq("absvalEqn", {"absEyEy": [200, 2500]}, 3, 1)
			# winMasker.addConstraintFreq("absvalEqn", {"absHyEy": [0, 40], "absExEy": [200, 2500]}, 4, 2)

			#################
			# 65536Hz masking
			################
			winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.20, 1.00], "cohEyHx": [0.20, 1.00]}, 0)			
			winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.20, 1.00], "cohEyHx": [0.20, 1.00]}, 1)			
			winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.20, 1.00], "cohEyHx": [0.20, 1.00]}, 2)			
			winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.20, 1.00], "cohEyHx": [0.20, 1.00]}, 3)
			winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.20, 1.00], "cohEyHx": [0.20, 1.00]}, 4)
			winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.20, 1.00], "cohEyHx": [0.20, 1.00]}, 5)
			winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.20, 1.00], "cohEyHx": [0.20, 1.00]}, 6)
			winMasker.addConstraintFreq("absvalEqn", {"absExEx": [0.0, 50000], "absExEy": [0.0, 60000], "absEyEx": [0.0, 100000], "absEyEy": [0.0, 200000]}, 0, 3)
			winMasker.addConstraintFreq("absvalEqn", {"absExEx": [0.0, 80000], "absExEy": [0.0, 60000], "absEyEx": [0.0, 100000], "absEyEy": [0.0, 200000]}, 0, 3)

			#################
			# 128Hz masking
			################
			## final constraints - along with using night data between 19:00 and 9:00 for decimation levels 0 and 1
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.60, 1.00], "cohEyHx": [0.70, 1.00]}, 2)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.50, 1.00], "cohEyHx": [0.70, 1.00]}, 3)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.40, 1.00], "cohEyHx": [0.60, 1.00]}, 4)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.30, 1.00], "cohEyHx": [0.30, 1.00]}, 5)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 5], "absExHyR": [0.0, 5], "absEyHxR": [0.0, 5], "absEyHyR": [0.0, 5]}, 1, 4)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 9], "absExHyR": [0.0, 10], "absEyHxR": [0.0, 10], "absEyHyR": [0.0, 10]}, 1, 5)		
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 0.3], "absExHyR": [0.0, 0.3], "absEyHxR": [0.0, 0.2], "absEyHyR": [0.0, 0.2]}, 2, 0)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 1], "absExHyR": [0.0, 0.5], "absEyHxR": [0.0, 0.4], "absEyHyR": [0.0, 0.5], "absHyHxR": [0, 0.02], "absHyHyR": [0, 0.015], "absHxHxR": [0, 0.03], "absHxHyR": [0, 0.02]}, 2, 1)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 5], "absExHyR": [0.0, 2], "absEyHxR": [0.0, 2], "absEyHyR": [0.0, 2], "absHyHxR": [0, 0.1], "absHyHyR": [0, 0.04], "absHxHxR": [0, 0.05], "absHxHyR": [0, 0.1]}, 2, 2)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 1], "absExHyR": [0.0, 1.3], "absEyHxR": [0.0, 2], "absEyHyR": [0.0, 2], "absHyHxR": [0, 0.04], "absHyHyR": [0, 0.07], "absHxHxR": [0, 0.1], "absHxHyR": [0, 0.1]}, 2, 3)					
			#### END OF M1


			### M13
			# #################
			# # 65536Hz masking
			# ################
			# before 2,0 looks good
			# winMasker.addConstraintFreq("absvalEqn", {"absEyEy": [1000, 20000]}, 2, 0)			
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.6, 1.0], "cohEyHx": [0.6, 1.0]}, 2, 1)
			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [1000, 200000], "absEyEx": [1000, 200000], "absExEy": [1000, 200000], "absEyEy": [1000, 200000]}, 2, 2)	
			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [1000, 200000], "absEyEx": [1000, 200000], "absExEy": [1000, 200000], "absEyEy": [1000, 200000]}, 2, 3)
			# winMasker.addConstraintFreq("absvalEqn", {"absExEx": [50, 2500], "absEyEx": [50, 2500], "absExEy": [50, 2500], "absEyEy": [50, 2500]}, 3, 0)
			# winMasker.addConstraintFreq("absvalEqn", {"absEyEy": [200, 2500]}, 3, 1)
			# winMasker.addConstraintFreq("absvalEqn", {"absHyEy": [0, 40], "absExEy": [200, 2500]}, 4, 2)

			#################
			# 128Hz masking
			################
			## final - along with using night data between 20:00 and 8:00 for decimation levels 0,1,2
			# winMasker.addConstraintFreq("cohStat", {"cohExHx": [0.00, 0.3], "cohEyHy": [0.00, 0.3]}, 0, 2)	
			# winMasker.addConstraintLevel("pcohStat", {"bivarEx": [0.60, 1.00], "bivarEy": [0.60, 1.00]}, 0)
			# winMasker.addConstraintLevel("pcohStat", {"bivarEx": [0.40, 1.00], "bivarEy": [0.60, 1.00]}, 1)
			# winMasker.addConstraintLevel("pcohStat", {"bivarEx": [0.40, 1.00], "bivarEy": [0.40, 1.00]}, 2)					
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.40, 1.00], "cohEyHx": [0.00, 1.00]}, 3,0)	
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.40, 1.00], "cohEyHx": [0.00, 1.00]}, 3,1)	
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.40, 1.00], "cohEyHx": [0.00, 1.00]}, 3,2)	
			# winMasker.addConstraintFreq("cohStat", {"cohExHy": [0.40, 1.00], "cohEyHx": [0.00, 1.00]}, 3,3)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 0.08], "absExHyR": [0.0, 0.08], "absEyHxR": [0.0, 0.1], "absEyHyR": [0.0, 0.1]}, 2, 0)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 0.08], "absExHyR": [0.0, 0.08], "absEyHxR": [0.0, 0.1], "absEyHyR": [0.0, 0.1]}, 2, 1)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 0.18], "absExHyR": [0.0, 0.22], "absEyHxR": [0.0, 0.17], "absEyHyR": [0.0, 0.25]}, 2, 2)			
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 0.2], "absExHyR": [0.0, 0.5], "absEyHxR": [0.0, 0.2], "absEyHyR": [0.0, 0.7]}, 2, 3)
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 0.5], "absExHyR": [0.0, 0.8], "absEyHxR": [0.0, 0.4], "absEyHyR": [0.0, 0.7]}, 2, 4)			
			# winMasker.addConstraintFreq("absvalEqn_RR", {"absExHxR": [0.0, 10], "absExHyR": [0.0, 20], "absEyHxR": [0.0, 20], "absEyHyR": [0.0, 50]}, 2, 5)			
			#### END OF M13


			# Single site coherence masking
			# cohStat
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.90, 1.00], "cohEyHx": [0.90, 1.00]}, 0)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.80, 1.00], "cohEyHx": [0.80, 1.00]}, 1)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.60, 1.00], "cohEyHx": [0.60, 1.00]}, 2)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.60, 1.00], "cohEyHx": [0.60, 1.00]}, 3)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.60, 1.00], "cohEyHx": [0.60, 1.00]}, 4)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.30, 1.00], "cohEyHx": [0.30, 1.00]}, 5)

			# cohStat2
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.95, 1.00], "cohEyHx": [0.95, 1.00]}, 0)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.90, 1.00], "cohEyHx": [0.90, 1.00]}, 1)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.80, 1.00], "cohEyHx": [0.80, 1.00]}, 2)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.70, 1.00], "cohEyHx": [0.70, 1.00]}, 3)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.60, 1.00], "cohEyHx": [0.60, 1.00]}, 4)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.30, 1.00], "cohEyHx": [0.30, 1.00]}, 5)

			# cohStat3
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.70, 1.00], "cohEyHx": [0.70, 1.00]}, 0)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.65, 1.00], "cohEyHx": [0.65, 1.00]}, 1)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.60, 1.00], "cohEyHx": [0.60, 1.00]}, 2)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.50, 1.00], "cohEyHx": [0.50, 1.00]}, 3)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.40, 1.00], "cohEyHx": [0.40, 1.00]}, 4)
			# winMasker.addConstraintLevel("cohStat", {"cohExHy": [0.30, 1.00], "cohEyHx": [0.30, 1.00]}, 5)

			# print information
			winMasker.printConstraints()
			# at the end
			winMasker.applyConstraints()
			winMasker.printMaskedWindows()
			winMasker.writeWindowFile(options["maskname"])
			# test with the window selector
			winSelector = WindowSelector(proj, fs, decParams, winParams)
			winSelector.setSites([s])
			winSelector.addWindowMask(s, options["maskname"])
			winSelector.calcSharedWindows()
			winSelector.printInfo()
			winSelector.printDatetimeConstraints()
			winSelector.printWindowMasks()
			winSelector.printSharedIndices()
			winSelector.printWindowsForFrequency()


def getDefaultOptions(proj):
	# default options
	default = {}
	default["sites"] = proj.getAllSites()
	default["freqs"] = proj.getAllSampleFreq()
	default["evalfreq"] = []
	default["declevels"] = 7
	default["freqlevel"] = 6
	default["stats"] = ["cohStat", "resPhaseStat"]
	default["maskname"] = "singleSite"
	return default

def parseKeywords(default, keywords):
	# check the user options
	for w in default:
		if w in keywords:
			default[w] = keywords[w]
	return default
