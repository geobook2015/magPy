# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop

class DataWriterInternal
Inherits from DataWriter
This is simply head files and binary data files
## In the future 
	The internal data format is in fact ATS
	With the 1024 header set at the beginning of the file for each channel
However, the header file saved is relevant only to this software
The header file means less processing to read the header information
"""
import os
import glob
import struct
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
import numpy as np
# parents
from dataWriter import DataWriter
# utils
from utilsIO import *

class DataWriterInternal(DataWriter):

	def writeDataFiles(self, chans, chanData):
		for idx, c in enumerate(chans):
			writePath = os.path.join(self.getOutPath(), "chan_{:02d}{}".format(idx, self.extension))
			dataF = open(writePath, "wb")
			# could put the ATS header here
			# not required right now
			# self.writeDataFileHeader()
			dataF.write(chanData[c].astype(self.dtype).tobytes())
			dataF.close()

	## write out the ATS header
	#def writeDataFileHeader(self, dataF, header):
	#	# header length, int16 in 0-1

	#	# header version, int16 in 2-3

	#	# number samples int32 in 4-7

	#	# float sampleFreq in Hz 

	#	# unix TIMESTAMP for starttime unsigned int

	###################
	### DEBUG
	###################		
	def printInfoBegin(self):
		self.printText("####################")	
		self.printText("DATA WRITER INTERNAL INFO BEGIN")		
		self.printText("####################")	

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("DATA WRITER INTERNAL INFO END")		
		self.printText("####################")			
		
	def printText(self, infoStr):
		generalPrint("Data Writer Internal Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Data Writer Internal Warning", warnStr)