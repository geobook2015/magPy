"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
spectrum writer writes fourier coefficients
Have two files
Info files have information about spectra
.dat files are ascii formatted data
.bin files are binary formatted data
"""
import os
import numpy as np
# utils
from utilsIO import *

class SpectrumWriter(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, dataRoot, refTime):
		# read and write spectra
		self.dataRoot = dataRoot
		self.filepath = ''
		self.refTime = refTime
		self.file = None

	###################
	### GET GENERAL INFO
	##################
	def getDataRoot(self):
		return self.dataRoot

	def getFilePath(self):
		return self.filepath

	def getFile(self):
		return self.file

	def getRefTime(self):
		return self.refTime

	###################
	### SET GENERAL INFO
	##################
	def setDataRoot(self, dataRoot):
		self.dataRoot = dataRoot		

	###################
	### WRITE FILES
	##################	
	def openBinaryForWriting(self, filename, fileInc, fs, winSize, winOverlap, globalOffset, numWindows, channels):
		# check to see if directory exists
		checkAndMakeDir(self.dataRoot)	
		# create filepath without extension
		filebase = self.getFileBase(filename, fileInc)
		# create info file
		filepathInfo = os.path.join(self.dataRoot, filebase + '.info')		
		self.writeInfoFile(filepathInfo, fs, winSize, winOverlap, globalOffset, numWindows, channels)
		# create file for data	
		self.filepath = os.path.join(self.dataRoot, filebase + '.bin')
		self.printText("Opening file {}".format(self.filepath))
		self.file = open(self.filepath, 'wb')		

	def writeBinary(self, data, dataInc):
		chans = list(data.keys())
		# sort channels alphebetically and save fourier coefficients		
		chans = sorted(chans)
		for c in chans:
			# this needs to be better - there are spaces in this file which are not required
			# this needs to be updated for newer version of numpy and scipy
			#self.file.write(data[c].tostring())
			#stack = np.hstack((data[c].real, data[c].imag))
			# write the real array first and then the imaginary array		
			#self.file.write(stack.tobytes())
			#float_array = array('d', stack)
			#float_array.tofile(self.file)
			# save as complex64 instead of 128 - otherwise too big
			self.file.write(data[c].astype('complex64').tobytes())

	def openAsciiForWriting(self, filename, fileInc, fs, winSize, winOverlap, globalOffset, numWindows, channels):
		# check to see if directory exists
		checkAndMakeDir(self.dataRoot)	
		# create filepath without extension
		filebase = filename + '{:02d}'.format(fileInc)
		# create info file
		filepathInfo = os.path.join(self.dataRoot, filebase + '.info')			
		self.writeInfoFile(filepathInfo, fs, winSize, winOverlap, globalOffset, numWindows, channels)			
		# create file for data	
		self.filepath = os.path.join(self.dataRoot, filebase + '.dat')
		self.printText("Opening file {}".format(self.filepath))		
		self.file = open(self.filepath, 'w')			

	def writeAscii(self, data, dataInc):
		chans = list(data.keys())
		# sort channels alphebetically and save fourier coefficients		
		chans = sorted(chans)
		for c in chans:
			outStr = arrayToStringSci(data[c])
			outStr = outStr + "\n"
			self.file.write(outStr)

	def writeInfoFile(self, filepath, fs, winSize, winOverlap, globalOffset, numWindows, channels):
		# open file
		infoFile = open(filepath, 'w')	
		# sort channels alphabetically - matching the order in the data files
		channels = sorted(channels)		
		# write out header information
		numChannels = len(channels)
		tmp = winSize + 1 # if winSize is odd, this will go down
		if winSize%2 == 0:
			tmp = tmp + 1
		dataSize = tmp/2
		infoFile.write("Reference time = {}\nSample frequency = {:.8f}\nWindow size = {:d}\nWindow overlap = {:d}\nGlobal offset = {:d}\n".format(
			self.getRefTime(), fs, winSize, winOverlap, globalOffset))	
		infoFile.write("Number of windows = {:d}\nData size = {:d}\n".format(numWindows, dataSize))
		infoFile.write("Number of channels = {:d}\n".format(numChannels))
		infoFile.write("Channels = " + " ".join(channels))
		infoFile.close()	

	###################
	### UTILS
	##################
	def getFileBase(self, filename, fileInc):
		return filename + '{:02d}'.format(fileInc)

	###################
	### CLOSE FILE
	##################
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
		self.printText("SPECTRUM WRITER INFO BEGIN")		
		self.printText("####################")		
		self.printText("Data root = {}".format(self.getDataRoot()))
		if self.getFile():
			self.printText("Current file open: {}".format(self.getFilePath()))
		self.printText("####################")
		self.printText("SPECTRUM WRITER INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Spectrum Writer Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Spectrum Writer Warning", warnStr)
