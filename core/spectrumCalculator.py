"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
The spectrumCalculator calculates different types of spectra
Fourier transform
"""
import numpy as np
import scipy.signal as signal
import pyfftw
# utils
from utilsFreq import *
from utilsIO import *

class SpectrumCalculator(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, fs, winSamples):
		# data reader
		self.fs = fs
		self.numSamples = winSamples
		self.window = True	
		self.windowFunc = signal.hamming(winSamples)
		# create an pyfftw plan
		self.dataArray = pyfftw.empty_aligned(self.numSamples, dtype="float64")		
		self.fftObj = pyfftw.builders.rfft(self.dataArray)

	###################
	### GET GENERAL INFO
	##################
	def getSampleFreq(self):
		return self.fs

	def getNumSamples(self):
		return self.numSamples

	def getWindow(self):
		return self.window

	def getWindowFunc(self):
		return self.windowFunc

	###################
	### SET WINDOW PERCENT
	##################	
	def setWindow(self, window):
		self.window = window	

	def setWindowFunc(self, windowFunc):
		self.windowFunc = windowFunc

	###################
	### FOURIER TRANSFORM
	##################
	def calcFourierCoeff(self, data):
		fftData = {}
		chans = list(data.keys())
		for c in chans:		
			# first copy data into dataArray
			self.dataArray[:] = data[c][:]		
			# no need to pad, these are usually multiples of two
			# remove mean and window if set
			# detrend
			self.dataArray[:] = signal.detrend(self.dataArray, type="linear")			
			if self.getWindow():
				self.dataArray[:] = self.dataArray*self.getWindowFunc()
			# use pytfftw here
			fftData[c] = np.array(self.fftObj()) 
		# calculate frequency array
		f = getFrequencyArray(self.getSampleFreq(), len(fftData[chans[0]]))
		return f, fftData	

	#def calcRobustFourierCoeff(self, data):
	# use a robust method to calculate fourier coefficients

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("SPECTRUM CALCULATOR INFO BEGIN")		
		self.printText("####################")		
		self.printText("Sampling frequency = {:f} [Hz]".format(self.fs))				
		self.printText("####################")
		self.printText("SPECTRUM CALCULATOR INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Spectrum Calculator Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Spectrum Calculator Warning", warnStr)
