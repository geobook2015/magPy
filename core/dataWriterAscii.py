# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop

class DataWriterAscii
Inherits from DataWriter
This is simply head files and ascii data files
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

class DataWriterAscii(DataWriter):

	def writeDataFiles(self, chans, chanData):
		self.extension = ".ascii"
		for idx, c in enumerate(chans):
			writePath = os.path.join(self.getOutPath(), "chan_{:02d}{}".format(idx, self.extension))
			# this could probably be made quicker - numpy savetxt maybe
			dataF = open(writePath, "w")
			size = chanData[c].size
			for i in xrange(0, size):
				dataF.write("{:9f}\n".format(chanData[c][i]))
			dataF.close()

	###################
	### DEBUG
	###################		
	def printInfoBegin(self):
		self.printText("####################")	
		self.printText("DATA WRITER ASCII INFO BEGIN")		
		self.printText("####################")	

	def printInfoEnd(self):
		self.printText("####################")
		self.printText("DATA WRITER ASCII INFO END")		
		self.printText("####################")			
		
	def printText(self, infoStr):
		generalPrint("Data Writer ASCII Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Data Writer ASCII Warning", warnStr)