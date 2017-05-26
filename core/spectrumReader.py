"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
spectrum reader reads fourier coefficients
Have two files
Info files have information about spectra
.dat files are ascii formatted data
.bin files are binary formatted data
"""
import os
import numpy as np
from datetime import datetime, timedelta
# utils
from utilsIO import *

class SpectrumReader(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, dataRoot):
		# read and write spectra
		self.dataRoot = dataRoot
		if not os.path.exists(self.dataRoot):
			self.printWarning("Directory {} does not exist".format(self.dataRoot))
		self.headerKeys = ["Reference time", "Sample frequency", "Window size", "Window overlap", 
			"Global offset", "Number of windows", "Number of channels", "Data size", "Channels"]		
		self.headers = {}
		self.dataType = np.dtype('complex64')
		self.dataByteSize = self.dataType.itemsize 
		self.filepath = ''
		self.file = False

	###################
	### GET GENERAL INFO
	##################
	def getDataRoot(self):
		return self.dataRoot

	def getHeaders(self):
		return self.headerKeys

	def getHeaderDict(self):
		return self.headers

	def getReferenceTime(self):
		return self.headers["Reference time"]

	def getSampleFreq(self):
		return float(self.headers["Sample frequency"])

	def getWindowSize(self):
		return int(self.headers["Window size"])

	def getWindowOverlap(self):
		return int(self.headers["Window overlap"])

	def getGlobalOffset(self):
		return int(self.headers["Global offset"])

	def getNumWindows(self):
		return int(self.headers["Number of windows"])

	def getGlobalRange(self):
		# the first global index is self.getGlobalOffset()
		# as the counting starts from zero 
		# the last is self.getGlobalOffset() + self.getNumWindows() - 1
		# because inclusive
		return [self.getGlobalOffset(), self.getGlobalOffset() + self.getNumWindows() - 1]

	def getNumChannels(self):
		return int(self.headers["Number of channels"])

	def getDataSize(self):
		return int(self.headers["Data size"])

	def getChannels(self):
		return self.headers["Channels"]	

	# return the frequency array
	def getFrequencyArray(self):
		return np.linspace(0, self.getSampleFreq()/2.0, self.getDataSize())

	###################
	### SET GENERAL INFO
	##################
	def setDataRoot(self, dataRoot):
		self.dataRoot = dataRoot	

	###################
	### READ FILES
	##################	
	def openBinaryForReading(self, filename, fileInc):
		filebase = self.getFileBase(filename, fileInc)
		filepathInfo = os.path.join(self.dataRoot, filebase + '.info')	
		self.filepath = os.path.join(self.dataRoot, filebase + '.bin')			
		# check files exist
		if not os.path.exists(filepathInfo) or not os.path.exists(self.filepath):
			self.printWarning("No data found in either {} or {}".format(filepathInfo, self.filepath))			
			return False
		# read info file	
		self.readInfoFile(filepathInfo)	
		self.channelByteSize = self.dataByteSize*self.getDataSize()
		self.windowByteSize = self.channelByteSize*self.getNumChannels()
		# set file to filepath - this is because binary files do not require opening using memap
		# self.file = self.filepath		
		return True		

	def readBinaryWindowLocal(self, localIndex):
		if localIndex >= self.getNumWindows():
			self.printWarning("Local index {:d} out of bounds".format(localIndex))
			self.printWarning("Min index = {:d}, Max index = {:d}".format(0, self.getNumWindows()-1))
		# with binary files, want the correct bytes
		byteOff = localIndex*self.windowByteSize
		specData = {}
		for cI, c in enumerate(self.getChannels()):
			chanOff = cI*self.channelByteSize
			specData[c] = np.memmap(self.filepath, dtype=self.dataType, mode='r', offset=byteOff + chanOff, shape=(self.getDataSize()))	
		return specData			
								

	def readBinaryWindowGlobal(self, globalIndex):
		if globalIndex >= self.getNumWindows() + self.getGlobalOffset() or globalIndex < self.getGlobalOffset():
			self.printWarning("Global index {:d} out of bounds".format(globalIndex))
			self.printWarning("Min index = {:d}, Max index = {:d}".format(self.getGlobalOffset(), self.getGlobalOffset() + self.getNumWindows()-1))	
		# convert global index to local index and return readAsciiWindowLocal						
		localIndex = globalIndex - self.getGlobalOffset()
		return self.readBinaryWindowLocal(localIndex)

	def openAsciiForReading(self, filename, fileInc):
		filebase = self.getFileBase(filename, fileInc)
		filepathInfo = os.path.join(self.dataRoot, filebase + '.info')			
		self.filepath = os.path.join(self.dataRoot, filebase + '.dat')	
		# check files exist
		if not os.path.exists(filepathInfo) or not os.path.exists(self.filepath):
			self.printWarning("No data found in either {} or {}".format(filepathInfo, self.filepath))
			return False
		# read info file	
		self.readInfoFile(filepathInfo)	
		# open file for reading
		self.file = open(self.filepath, "rb")
		# run through and find line endings
		self.lineOffset = []
		offset = 0
		for line in self.file:
		    self.lineOffset.append(offset)
		    offset += len(line)
		self.file.seek(0)		
		return True				

	def readAsciiWindowLocal(self, localIndex):
		if localIndex >= self.getNumWindows():
			self.printWarning("Local index {:d} out of bounds".format(localIndex))
			self.printWarning("Min index = {:d}, Max index = {:d}".format(0, self.getNumWindows()-1))
		# with ascii files, want the correct lines
		# find line where local index starts
		windowStartLine = localIndex*self.getNumChannels()
		specData = {}
		for cI, c in enumerate(self.getChannels()):
			indexC = windowStartLine + cI
			self.file.seek(self.lineOffset[indexC])
			line = self.file.readline()
			specData[c] = np.loadtxt(line.strip().split(","), dtype=complex)
		return specData	

	def readAsciiWindowGlobal(self, globalIndex):
		if globalIndex >= self.getNumWindows() + self.getGlobalOffset() or globalIndex < self.getGlobalOffset():
			self.printWarning("Global index {:d} out of bounds".format(globalIndex))
			self.printWarning("Min index = {:d}, Max index = {:d}".format(self.getGlobalOffset(), self.getGlobalOffset() + self.getNumWindows()-1))	
		# convert global index to local index and return readAsciiWindowLocal
		localIndex = globalIndex - self.getGlobalOffset()
		return self.readAsciiWindowLocal(localIndex)						
	
	def readInfoFile(self, filepath):
		# open file
		infoFile = open(filepath, 'r')	
		lines = infoFile.readlines()
		infoFile.close()
		# loop through all headers and get values
		for h in self.getHeaders():
			for l in lines:
				if h in l:
					self.headers[h] = self.getInfoValue(h, l)
					break 	

	###################
	### UTILS
	##################
	def getFileBase(self, filename, fileInc):
		return filename + "{:02d}".format(fileInc)	

	def getInfoValue(self, header, line):
		split = line.split("=")
		split[1] = split[1].strip()
		if header == "Channels":
			return split[1].split()
		elif header == "Reference time":
			return datetime.strptime(split[1], '%Y-%m-%d %H:%M:%S')
		elif header == "Sample frequency":
			return float(split[1])
		else:
			return int(float(split[1]))

	def closeFile(self):
		if self.filepath != '' and self.file:
			self.printText("Closing file {}".format(self.filepath))
			self.file.close()
			self.filepath = ''
		else:
			print 'No file open'	

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("SPECTRUM READER INFO BEGIN")		
		self.printText("####################")		
		self.printText("Data root = {}".format(self.getDataRoot()))
		if len(self.getHeaderDict()) > 0:
			self.printText("Filepath = {}".format(self.filepath))
			for h in self.getHeaders():
				self.printText("{} = {}".format(h, self.getHeaderDict()[h]))
		self.printText("####################")
		self.printText("SPECTRUM READER INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Spectrum Reader Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Spectrum Reader Warning", warnStr)
