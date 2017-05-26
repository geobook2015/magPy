"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
The calibrator takes a directory of calibration files 
A data object 
And does the calibration
"""
import os
import glob
import math
import scipy.interpolate as interp
import timeit
# utils
from utilsCalibrate import *
from utilsFreq import *
from utilsIO import *

class Calibrator(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, calDirectory):
		# data reader
		self.calDir = calDirectory
		self.findCalFiles()
		# set whether to use theoretical for magnetic channels
		self.useTheoretical = False

	###################
	### GET GENERAL INFO
	##################
	def getCalDirectory(self):
		return self.calDir

	def getCalFiles(self):
		return self.calFiles

	def getCalFormats(self):
		return self.calTypes	

	def getNumCalFiles(self):
		return self.numCalFiles	

	def getTheoretical(self):
		return self.useTheoretical

	###################
	### SET FUNCTION
	##################
	def setCalDir(self, calDirectory):
		self.calDir = calDirectory
		self.findCalFiles()

	def setTheoretical(self, theoretical):
		self.useTheoretical = theoretical

	###################
	### CALIBRATE
	##################	
	# calibrate data in the time domain
	# need to do a FFT
	# apply sensor TF
	# and then do an inverse FFT
	def calibrate(self, data, fs, sensor, serial, chopper):
		# iterate over data
		for chan in data:
			# output some info
			self.printText("Calibrating channel {}".format(chan))
			# try and find calibration file
			calIndex = self.getCalFile(sensor[chan], serial[chan], chopper[chan])
			if calIndex < 0:
				# no file found
				if self.getTheoretical() and (chan in magChannelsList()):
					# use theoretical
					calData = self.getTheoreticalCalData(sensor[chan])
					data[chan] = self.calibrateChan(data[chan], fs, calData)
					continue
				else:
					self.printText('No Calibration data found - channel will not be calibrated')	
					continue # nothing to do
			# else file found		
			# NOTE: no need to separately apply static gain, already included in cal data. Returned so that it can be output to provide info to user		
			calData, staticGain = readCalFile(self.calTypes[calIndex], self.calFiles[calIndex], sensor[chan], serial[chan], chopper[chan])	
			self.printText("Calibration file found for sensor {}, serial number {}, chopper {}: {}".format(sensor[chan], serial[chan], chopper[chan], self.calFiles[calIndex]))
			self.printText("Static gain correction of {} applied to calibration data".format(staticGain))
			self.printText("Format: {}".format(self.calTypes[calIndex]))									
			data[chan] = self.calibrateChan(data[chan], fs, calData)	
		# return data at the end
		return data

	# calibrate data in frequency domain
	# need to interpolate calibration info to frequencies in data
	# and then apply TF
	def calibrateChan(self, data, fs, calData):
		# do the forward transform
		dataSize = data.size
		# pad end of array
		data = np.pad(data, (0, padNextPower2(dataSize)), 'constant')
		fftData = forwardFFT(data)
		f = getFrequencyArray(fs, fftData.size) 	
		# now do the calibration in frequency domain
		transferFunc = self.interpolateCalData(calData, f)	
		# recall, zero element excluded, so start from 1
		# fft zero element should really be zero because average is removed from data
		fftData[1:] = fftData[1:]/transferFunc		
		# return the non padded part of the array		
		return inverseFFT(fftData, data.size)[:dataSize]

	def interpolateCalData(self, calData, f):
		# interpolate phase and magnitude
		# interpFuncMag = interp.interp1d(calData[:,0], calData[:,1])
		# interpFuncPhase = interp.interp1d(calData[:,0], calData[:,2])
		# try spline interpolation: InterpolatedUnivariateSpline
		interpFuncMag = interp.InterpolatedUnivariateSpline(calData[:,0], calData[:,1])
		interpFuncPhase = interp.InterpolatedUnivariateSpline(calData[:,0], calData[:,2])
		return interpFuncMag(f[1:])*np.exp(1j*interpFuncPhase(f[1:]))

	###################
	### GET CAL DATA
	##################
	def getCalFile(self, sensor, serial, chopper):
		index = -1
		for calE, calF in zip(self.calExt, self.calFormats):
			# get the name for this format
			calName = getCalName(calF, calE, sensor, serial, chopper)			
			# search to find a calibration file with that name
			# take the first encountered
			for iF in xrange(0, self.getNumCalFiles()):
				if calName in self.calFiles[iF]:
					index = iF
			# break out of calfile found
			if index != -1:
				break
		# return the index to the calibration file
		return index

	def getTheoreticalCalData(self, sensor):
		# should use the sensor to figure out what calibration func
	    if "mfs06" in sensor or "MFS06" in sensor:
	    	return unitCalibration()


	###################
	### GET THE CALIBRATION FILES IN THE DIRECTORY
	##################
	def findCalFiles(self):
		self.calFiles = []
		self.calTypes = []
		self.calExt, self.calFormats = getKnownCalibrationFormats()
		for cE, cF in zip(self.calExt, self.calFormats):
			# get all files of format cE
			tmp1 = glob.glob(os.path.join(self.calDir,"*.{}".format(cE)))
			tmp2 = [cF]*len(tmp1)
			self.calFiles = self.calFiles + tmp1
			self.calTypes = self.calTypes + tmp2	
		self.numCalFiles = len(self.calFiles)

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("CALIBRATOR INFO BEGIN")		
		self.printText("####################")		
		self.printText("Known extensions and calibration formats")
		self.printText("Extensions:")
		self.printText(", ".join(self.calExt))
		self.printText("Formats")
		self.printText(", ".join(self.calFormats))
		# now print the actual calibration files
		self.printText("Calibration files found:")
		for cE, cT in zip(self.calFiles, self.calTypes):
			self.printText("{}\t\t{}".format(cT, cE))
		self.printText("####################")
		self.printText("CALIBRATOR INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Calibrator Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Calibrator Warning", warnStr)
