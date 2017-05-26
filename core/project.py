# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop

The project holds all the information about the project
And imposes the project structure
Time Data - any time data
Spec Data - this is the fourier data
Stat Data - statistic data
TF Data - transfer functions
Cal Data - the path to the calibration files
"""
import os
from datetime import datetime, timedelta
# utils
from utilsIO import *
from utilsReader import *

class Project(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self):
        # project file path
		self.projectFile = ""
		# paths to various parts of the project
		self.timeDataPath = ""
		self.specDataPath = ""
		self.statDataPath = ""
		self.transDataPath = ""
		self.calDataPath = ""
		self.imageDataPath = ""
		self.refTime = None
		self.projStart = datetime.now()
		self.projEnd = datetime.now()
		# site information
		self.sites = []
		self.sitesStarts = {}
		self.sitesEnds = {}
		self.sitesFs = {}
		self.sitesOcc = {}
		self.sitesTime = {}	
		self.sitesSpec = {}
		self.sitesStats = {}
		self.sitesTransFuncs = {}
		# cal files list
		self.calFiles = []
		# save readers - depending on speed
		self.sitesReaders = {}

	###################
	### GET VARIOUS PROJECT PROPERTIES
	##################
	def getTimeDataPath(self):
		return self.timeDataPath

	def getSpecDataPath(self):
		return self.specDataPath

	def getStatDataPath(self):
		return self.statDataPath

	def getTransDataPath(self):
		return self.transDataPath
	
	def getCalDataPath(self):
		return self.calDataPath

	def getImageDataPath(self):
		# data path to store images
		return self.imageDataPath
	
	def getRefTime(self):
		return self.refTime

	def getProjectStart(self):
		return self.projStart

	def getProjectStop(self):
		return self.projEnd
	
	def getAllSites(self):
		return self.sites

	def getNumSites(self):
		return len(self.sites)

	def getAllSampleFreq(self):
		sampleFreq = set()
		for s in self.getAllSites():
			sampleFreq.update(self.getSiteSampleFreqs(s))
		return sorted(list(sampleFreq))

 	###################
	### GET INFORMATION ABOUT 
    ### A SINGLE SITE
	##################      
	def getTimeDataPathSite(self, site):
		return os.path.join(self.getTimeDataPath(), site)

	def getSpecDataPathSite(self, site):
		return os.path.join(self.getSpecDataPath(), site)

	def getStatDataPathSite(self, site):
		return os.path.join(self.getStatDataPath(), site)

	def getTransDataPathSite(self, site):
		return os.path.join(self.getTransDataPath(), site)

	def getSiteTimeFiles(self, site):		
		return self.sitesTime[site]

	# get raw files of certain sample frequency
	def getSiteTimeFilesFs(self, site, fs):	
		dataFiles = []
		for measdir in self.getSiteTimeFiles(site):
			if self.getMeasSampleFreq(site, measdir) == fs:
				dataFiles.append(measdir)
		return dataFiles

	def getSiteSampleFreqs(self, site):
		# find the distinct samples frequencies
		fs = []
		for measdir in self.getSiteTimeFiles(site):
			measFs = self.getMeasSampleFreq(site, measdir)
			if measFs not in fs:
				fs.append(measFs)
		return fs

	def getSiteProcFiles(self, site):
		return self.sitesProc[site]
		
	def getSiteSpectraFiles(self, site):
		return self.sitesSpec[site]
		
	def getSiteStatsFiles(self, site):
		return self.sitesStats[site]
		
	def getSiteTransFuncFiles(self, site):
		return self.sitesTransFuncs[site]
		
	def getSiteStart(self, site):
		return self.getSiteOcc(site)[0]

	def getSiteStop(self, site):
		return self.getSiteOcc(site)[1]

	def getSiteOcc(self, site):
		return self.sitesOcc[site]

	# find sites which were occupied at the same time
	# might be useful to multi-site processing
	def getSiteConcurrent(self, site):
		occTimes = self.getSiteOcc(site)
		concSites = []
		for s in self.getAllSites():
			if s == site:
				continue
			occTimesSite = self.getSiteOcc(s)
			# check that start time is before end time of site and that end time is after start time of site
			if (occTimesSite[0] < occTimes[1]) and (occTimesSite[1] > occTimes[0]):
				concSites.append(s)
		return concSites		

 	###################
	### GET INFORMATION ABOUT 
    ### A SINGLE MEASUREMENT
	##################
	def getTimeDataPathMeas(self, site, measdir):
		return os.path.join(self.getTimeDataPath(), site, measdir)

	def getProcDataPathMeas(self, site, measdir):
		return os.path.join(self.getProcDataPath(), site, measdir)

	def getSpecDataPathMeas(self, site, measdir):
		return os.path.join(self.getSpecDataPath(), site, measdir)

	def getStatDataPathMeas(self, site, measdir):
		return os.path.join(self.getStatDataPath(), site, measdir)

	def getTransDataPathMeas(self, site, measdir):
		return os.path.join(self.getTransDataPath(), site, measdir)

	def getMeasStartTime(self, site, measdir):
		return self.sitesStarts[site][measdir]

	def getMeasEndTime(self, site, measdir):
		return self.sitesEnds[site][measdir]

	def getMeasSampleFreq(self, site, measdir):
		return self.sitesFs[site][measdir]

	def getMeasDataReader(self, site, measdir):
		if (self.checkSite(site) and self.checkMeasdir(site, measdir)):
			return self.sitesReaders[site][measdir]
		self.printText("Data directory {} for site {} not found".format(site, measdir))
		return False

	###################
	### SET VARIOUS PROJECT PROPERTIES
	##################
	def setTimeDataPath(self, path):
		self.timeDataPath = path

	def setReformatDataPath(self, path):
		self.reformatDataPath = path

	def setProcDataPath(self, path):
		self.procDataPath = path

	def setSpecDataPath(self, path):
		self.specDataPath = path

	def setStatDataPath(self, path):
		self.statDataPath = path	

	def setTransDataPath(self, path):
		self.transDataPath = path

	def setCalDataPath(self, path):
		self.calDataPath = path

	def setImageDataPath(self, path):
		self.imageDataPath = path
		
	def setRefTime(self, refTime):
		self.refTime = datetime.strptime(refTime, '%Y-%m-%d %H:%M:%S')

	def setProjectFile(self, filepath):
		self.projectFile = filepath

	###################
	### SAVE AND LOAD PROJECT FILES
	##################		
	def loadProjectFile(self, filepath):
		f = open(filepath, "r")
		lines = f.readlines()
		f.close()
		for l in lines:
			split = l.strip().split("=")
			if "Time data path" in l:
				self.setTimeDataPath(split[1].strip())
			if "Spectra data path" in l:
				self.setSpecDataPath(split[1].strip())
			if "Statistics data path" in l:
				self.setStatDataPath(split[1].strip())
			if "TransFunc data path" in l:
				self.setTransDataPath(split[1].strip())
			if "Calibration data path" in l:
				self.setCalDataPath(split[1].strip())	
			if "Image data path" in l:
				self.setImageDataPath(split[1].strip())			
			if "Reference time" in l:
				self.setRefTime(split[1].strip())
		self.setProjectFile(filepath)				
		self.initialiseProject()
	
	def saveProjectFile(self, filepath):
		self.setProjectFile(filepath)
		f = open(filepath, "w")
		f.write("Time data path = {}\n".format(self.getTimeDataPath()))
		if self.getSpecDataPath() != '':			
			f.write("Spectra data path = {}\n".format(self.getSpecDataPath()))
		if self.getStatDataPath() != '':			
			f.write("Statistics data path = {}\n".format(self.getStatDataPath()))
		if self.getTransDataPath() != '':			
			f.write("TransFunc data path = {}\n".format(self.getTransDataPath()))
		if self.getCalDataPath() != '':			
			f.write("Calibration data path = {}\n".format(self.getCalDataPath()))
		if self.getImageDataPath() != '':			
			f.write("Image data path = {}\n".format(self.getImageDataPath()))						
		f.write("Reference time = {}\n".format(self.getRefTime().strftime('%Y-%m-%d %H:%M:%S')))
		f.close()

	###################
	### GET INFORMATION ABOUT PROJECT
	##################	
	def initialiseProject(self):
		# sites are given by directors under the raw data path
		if self.getTimeDataPath() == '':
			self.printText("No time data path")
			return
		self.getTimeInfo()
		self.getSpecInfo()
		self.getStatInfo()
		self.getTransInfo()
		self.getImageInfo()
		self.getCalInfo()

	# suppose new time directories have been added (saved)
	# this will refresh the project so that they are included
	def refreshProject(self):
		self.getTimeInfo()
		self.getSpecInfo()
		self.getStatInfo()
		self.getTransInfo()
		self.getCalInfo()		
	
	# get the data for the individual bits
	def getTimeInfo(self):
		checkAndMakeDir(self.getTimeDataPath())
		self.sites = getDirsInDirectory(self.getTimeDataPath())
		self.sites = sorted(self.sites)
		self.sitesTime = self.getSiteFiles(self.getTimeDataPath())
		self.getTimeReaders()

	def getTimeReaders(self):
		for s in self.getAllSites():
			self.sitesReaders[s] = {}
			self.sitesStarts[s] = {}
			self.sitesEnds[s] = {}
			self.sitesFs[s] = {}

			# get time files
			timeFiles = self.getSiteTimeFiles(s)
			if len(timeFiles) == 0:
				# no data in directory
				self.sitesOcc[s] = [self.getRefTime(), self.getRefTime()]
				continue				

			for d in timeFiles:
				# initialise datapath and the various readers
				datapath = os.path.join(self.getTimeDataPath(), s, d)
				reader = getDataReader(datapath)
				if not reader: # then no data found in this directory
					self.printWarning("No data found in measurement directory {} for site {}. Please remove this directory".format(d, s))
					continue
				# otherwise get the reader and information	
				self.sitesReaders[s][d] = reader
				self.sitesStarts[s][d] = reader.getStartDatetime()
				self.sitesEnds[s][d] = reader.getStopDatetime()
				self.sitesFs[s][d] = reader.getSampleFreq()
			occStart = min(self.sitesStarts[s].values())
			occEnd = max(self.sitesEnds[s].values())
			self.sitesOcc[s] = [occStart, occEnd]
			occTimes = self.sitesOcc.values()
			starts = []
			ends = []
			for times in occTimes:
				starts.append(times[0])
				ends.append(times[1])
			self.projStart = min(self.projStart, min(starts))
			self.projEnd = max(self.projEnd, max(ends))
	
	def getSpecInfo(self):
		if self.getSpecDataPath() == '':
			return		
		checkAndMakeDir(self.getSpecDataPath())
		self.sitesSpec = self.getSiteFiles(self.getSpecDataPath())			
	
	def getStatInfo(self):
		if self.getStatDataPath() == '':
			return		
		checkAndMakeDir(self.getStatDataPath())		
		self.sitesStat = self.getSiteFiles(self.getStatDataPath())				
	
	def getTransInfo(self):
		if self.getStatDataPath() == '':
			return			
		checkAndMakeDir(self.getTransDataPath())		
		self.sitesTrans = self.getSiteFiles(self.getTransDataPath())				
		
	def getCalInfo(self):
		if self.getCalDataPath() == '':
			return			
		checkAndMakeDir(self.getCalDataPath())
		self.printText("Note: Calibration data needs to be at the top level in directory {}".format(self.getCalDataPath()))	
		self.calFiles = getFilesInDirectory(self.getCalDataPath())

	def getImageInfo(self):
		if self.getImageDataPath() == '':
			return		
		checkAndMakeDir(self.getImageDataPath())

	def getSiteFiles(self, path):
		siteDic = {}
		for s in self.getAllSites():
			path2 = os.path.join(path,s)
			siteDic[s] = getDataDirsInDirectory(path2)
		return siteDic

	###################
	### CHECK FNCS
	##################
	def checkSite(self, site):
		if site in self.getAllSites():
			return True
		else:
			return False

	def checkMeasdir(self, site, measdir):
		if measdir in self.getSiteTimeFiles(site):
			return True
		else:
			return False			

	###################
	### DEBUG
	##################		
	# print headers
	def printInfo(self):		
		# print the headers
		self.printText("####################")
		self.printText("PROJECT INFO BEGIN")		
		self.printText("####################")		
		self.printText("Time data path = {}".format(self.getTimeDataPath()))
		self.printText("Spectra data path = {}".format(self.getSpecDataPath()))
		self.printText("Statistics data path = {}".format(self.getStatDataPath()))
		self.printText("TransFunc data path = {}".format(self.getTransDataPath()))
		self.printText("Calibration data path = {}".format(self.getCalDataPath()))	
		self.printText("Images data path = {}".format(self.getImageDataPath()))	
		self.printText("Reference time = {}".format(self.getRefTime().strftime('%Y-%m-%d %H:%M:%S')))
		self.printText("Project start time = {}".format(self.getProjectStart().strftime('%Y-%m-%d %H:%M:%S.%f')))
		self.printText("Project stop time = {}".format(self.getProjectStop().strftime('%Y-%m-%d %H:%M:%S.%f')))	
		self.printSiteOccTimes()			
		self.printText("####################")
		self.printText("PROJECT INFO END")		
		self.printText("####################")	

	def printSiteOccTimes(self):
		self.printText("Project found {} sites:".format(self.getNumSites()))
		for s in self.getAllSites():
			self.printText("{}\t\tstart: {}\tend: {}".format(s, self.getSiteStart(s), self.getSiteStop(s)))

	def printSiteInfo(self, site):
		self.printText("####################")
		self.printText("SITE INFO BEGIN")		
		self.printText("####################")		
		self.printText("Site = {}".format(site))
		self.printText("Time data path = {}".format(self.getTimeDataPathSite(site)))			
		self.printText("Spectra data path = {}".format(self.getSpecDataPathSite(site)))
		self.printText("Statistics data path = {}".format(self.getStatDataPathSite(site)))
		self.printText("TransFunc data path = {}".format(self.getTransDataPathSite(site)))	
		self.printText("Site start time = {}".format(self.getSiteStart(site).strftime("%Y-%m-%d %H:%M:%S.%f")))
		self.printText("Site stop time = {}".format(self.getSiteStop(site).strftime("%Y-%m-%d %H:%M:%S.%f")))	
		self.printText("Sampling frequencies recorded = {}".format(arrayToString(self.getSiteSampleFreqs(site))))
		# print any measurement directors
		timeFiles = self.getSiteTimeFiles(site)
		if len(timeFiles) > 0:
			self.printText("Time data files (measurement files) = {}".format(", ".join(timeFiles)))	
		else:
			self.printText("No time data files (measurement files) found")
		# list any concurrent sites
		if len(self.getSiteConcurrent(site)) == 0:
			self.printText("Concurrent sites = None")	
		else:
			self.printText("Concurrent sites = {}".format(", ".join(self.getSiteConcurrent(site))))	
		self.printText("####################")
		self.printText("SITE INFO END")		
		self.printText("####################")

	def printMeasInfo(self, site, measdir):
		# simply print the ats information
		self.printText("####################")
		self.printText("MEASUREMENT INFO BEGIN")		
		self.printText("####################")	
		self.printText("Measurement = {}".format(self.getTimeDataPathMeas(site, measdir)))			
		self.getMeasDataReader(site, measdir).printInfo()		
		self.printText("####################")
		self.printText(" MEASUREMENT INFO END")		
		self.printText("####################")			

	def printText(self, infoStr):
		generalPrint("Project Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Project Warning", warnStr)




