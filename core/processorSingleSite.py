"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
ProcessorSingleSite
Inherits from Processor
And implements single site processing
"""
import os
import random
import numpy as np
import scipy.signal as signal
import scipy.interpolate as interp
# utils
from utilsFreq import *
from utilsIO import *
from utilsRobust import *
from utilsProcess import *

class ProcessorSingleSite(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self, project, winSelector):
		# data reader
		self.proj = project
		self.winSelector = winSelector
		self.decParams = winSelector.getDecParams()
		self.winParams = winSelector.getWinParams()
		self.setDefaults()

	def setDefaults(self):
		# inputs
		self.inSite = ""
		self.inChannels = ["Hx", "Hy"]
		self.outSite = ""
		self.outChannels = ["Ex", "Ey"]
		self.allChannels = self.inChannels + self.outChannels
		# set cross channels - these are the actual channels to use in the solution
		self.crossChannels = self.allChannels
		# evaluation frequency data
		self.evalFreq = []
		self.evalFreqEqns = []
		# whether to include an intercept term
		self.intercept = False
		# smoothing options
		self.win = "hanning"
		self.winSmooth = -1
		# output filename
		self.prepend = ""


	###################
	### GET GENERAL INFO
	##################
	def getWinSelector(self):
		return self.winSelector

	def getDecParams(self):
		return self.decParams

	def getWinParams(self):
		return self.winParams

	def getInChannels(self):
		return self.inChannels

	def getInSite(self):
		return self.inSite

	def getInSize(self):
		return self.inSize

	def getOutChannels(self):
		return self.outChannels

	def getOutSite(self):
		return self.outSite

	def getOutSize(self):
		return self.outSize

	def getAllChannels(self):
		return self.allChannels

	def getCrossChannels(self):
		return self.crossChannels

	def getPrepend(self):
		return self.prepend

	def getIntercept(self):
		return self.intercept

	def getWindow(self):
		return self.win

	def getWindowSmooth(self, **kwargs):
		# first check if window size specified by user
		if self.winSmooth != -1 and self.winSmooth > 1:
			return self.winSmooth
		# if not, then calculate based on datasize
		if "datasize" in kwargs:
			winSmooth = kwargs["datasize"]*1.0/16.0
			if winSmooth < 3:
				return 3 # minimum smoothing
			# otherwise round to nearest odd number
			winSmooth = np.ceil(winSmooth) // 2 # this is floor division
			return int(winSmooth*2 + 1)
		# otherwise, return a default value
		return 15

	def getWinSelector(self):
		return self.winSelector

	###################
	### SET INPUTS AND OUTPUTS
	##################
	def setInput(self, inSite, inChannels):
		self.inSite = inSite
		self.inChannels = inChannels
		self.inSize = len(inChannels)
		# set all and cross channels
		self.allChannels = self.inChannels + self.outChannels
		self.crossChannels = self.allChannels

	def setOutput(self, outSite, outChannels):
		self.outSite = outSite
		self.outChannels = outChannels
		self.outSize = len(outChannels)
		# set all and cross channels
		self.allChannels = self.inChannels + self.outChannels
		self.crossChannels = self.allChannels

	def setCrossChannels(self, crossChannels):
		# this is for the case where we want to limit the cross channels used
		self.crossChannels = crossChannels 		

	def setPrepend(self, prepend):
		self.prepend = prepend

	def setSmooth(self, window, windowSize):
		self.win = window
		self.winSmooth = windowSize

	def setIntercept(self, intercept):
		self.intercept = intercept

	###################
	### PROCESS
	###################
	def process(self):
		# checks
		# for testing
		evalFreqEqnsTest = []
		evalFreqEqnsTest2 = []
		evalFreqEqnsTest3 = []
		evalFreqEqnsTest4 = []
		evalFreqVarsTest4 = []		
		evalFreqEqnsTest5 = []
		evalFreqEqnsTest6 = []
		self.crossChannels = ["Ex", "Ey"]
		# for each decimation level
		# read in the shared windows from all sites
		# for each evaluation frequency, store the data from each window
		# and then at the end, perform robust processing
		numLevels = self.getDecParams().getNumLevels()
		inChans = self.getInChannels()
		outChans = self.getOutChannels()
		# get set of shared windows for all decimation levels
		sharedWindows = self.getWinSelector().getSharedWindows()

		for iDec in xrange(0, numLevels):
			# print out some info
			self.printText("Processing decimation level {}".format(iDec))

			fs = self.getDecParams().getSampleFreqLevel(iDec)
			# get the number of all shared windows and the number of unmasked windows
			# unmasked windows are ones that will actually be used in the calculation
			numWindows = self.getWinSelector().getNumSharedWindows(iDec)
			unmaskedWindows = self.getWinSelector().getUnmaskedWindowsLevel(iDec)
			numUnmasked = len(unmaskedWindows)
			self.printText("Total shared windows for decimation level = {}".format(numWindows))
			self.printText("Total unmasked windows for decimation level = {}".format(numUnmasked))
			if numUnmasked == 0:
				self.printText("No unmasked windows found at this decimation level, continuing to next level".format(iDec))
				continue # continue to next decimation level
			self.printText("{} windows will be processed".format(numUnmasked))

			# get the evaluation frequencies
			evalFreq = self.getDecParams().getEvalFrequenciesForLevel(iDec)
			# set some variables
			totalSize = self.getInSize() + self.getOutSize()
			numEvalFreq = len(evalFreq)
			dataSize = self.getWinSelector().getDataSize(iDec)
			freq = np.linspace(0, fs/2, dataSize)
			# get the window smoothing params
			smoothLen = self.getWindowSmooth(datasize=dataSize)

			# create the data array
			# for each evaluation frequency
			# keep the spectral power information for all windows
			evalFreqData = np.empty(shape=(numEvalFreq, numWindows, totalSize, totalSize), dtype="complex")

			# an array for the window data
			winSpectraMatrix = np.empty(shape=(totalSize, totalSize, dataSize), dtype="complex")
			winDataArray = np.empty(shape=(totalSize, dataSize), dtype="complex")

			# loop over shared windows
			localWin = 0
			global2local = {}
			# randomise the windows
			# sharedWin = random.shuffle(sharedWindows[iDec])
			for iWin in unmaskedWindows:
				# do the local to global map
				global2local[iWin] = localWin

				# get the window for the input site
				inSF, inReader = self.getWinSelector().getSpecReaderForWindow(self.getInSite(), iDec, iWin)
				inData = inReader.readBinaryWindowGlobal(iWin)
				# get the window and channels for the output site
				if self.getOutSite() != self.getInSite():
					outSF, outReader = self.getWinSelector().getSpecReaderForWindow(self.getOutSite(), iDec, iWin)
					outData = outReader.readBinaryWindowGlobal(iWin)
				else:
					outData = inData

				# get data into the right part of the arrays
				for i in xrange(0, self.getInSize()):
					winDataArray[i] = inData[inChans[i]]
				for i in xrange(0, self.getOutSize()):
					winDataArray[self.getInSize() + i] = outData[outChans[i]]

				# and now can fill the parts of the matrix
				# recall, smooth the power spectra
				for i in xrange(0, totalSize):
					for j in xrange(i, totalSize):
						# winSpectraMatrix[i,j] = winDataArray[i] * np.conjugate(winDataArray[j])
						winSpectraMatrix[i,j] = smooth1d(winDataArray[i] * np.conjugate(winDataArray[j]), smoothLen, self.getWindow())
						if i != j:
							winSpectraMatrix[j,i] = np.conjugate(winSpectraMatrix[i,j]) # due to complex symmetry

				# after running through all windows, calculate evaluation frequencies
				# calculate frequency array
				evalFreqData[:, localWin] = self.calcEvalFrequencyData(freq, evalFreq, winSpectraMatrix)

				# increment local window
				localWin = localWin + 1

			# now all the data has been collected
			# for each evaluation frequency, do the robust processing
			# and get the evaluation frequency data
			for eIdx in xrange(0, numEvalFreq):
				self.printText("Processing evaluation frequency = {:.6f} [Hz], period = {:.6f} [s]".format(evalFreq[eIdx], 1/evalFreq[eIdx]))
				# get the constrained windows for the evaluation frequency
				evalFreqWindows = self.getWinSelector().getWindowsForFreq(iDec, eIdx)
				if len(evalFreqWindows) == 0: # no windows meet constraints
					self.printText("No windows found - possibly due to masking")
					continue
				localWinIndices = []
				for iW in evalFreqWindows:
					localWinIndices.append(global2local[iW])
				self.printText("{:d} windows will be solved for".format(len(localWinIndices)))
				# restrict processing to data that meets constraints for this evaluation frequency
				# add to class vars
				self.evalFreq.append(evalFreq[eIdx])
				# try a random smoothing of the spectral density estimates across windows
				# spectralEstimates = self.smoothSpectralEstimates(evalFreqData[eIdx, localWinIndices])
				# self.evalFreqEqns.append(self.robustProcess(spectralEstimates))
				# solution using all components
				numSolveWindows, obs, reg = self.prepareLinearEqn2(evalFreqData[eIdx, localWinIndices])
				self.evalFreqEqns.append(self.robustProcess(numSolveWindows, obs, reg))
				# evalFreqEqnsTest.append(self.stackedProcess(evalFreqData[eIdx, localWinIndices]))
				# evalFreqEqnsTest2.append(self.robustProcessReduced(evalFreqData[eIdx, localWinIndices]))
				# evalFreqEqnsTest3.append(self.robustProcessOLS(numSolveWindows, obs, reg))
				tmp1, tmp2 = self.robustProcessCM(numSolveWindows, obs, reg)
				evalFreqEqnsTest4.append(tmp1)
				evalFreqVarsTest4.append(tmp2)
				# evalFreqEqnsTest4.append(self.robustProcessCM(numSolveWindows, obs, reg))
				# evalFreqEqnsTest5.append(self.robustProcessCMMod(numSolveWindows, obs, reg))
				# evalFreqEqnsTest6.append(self.robustProcessStack(numSolveWindows, obs, reg))

		# write out all the data
		self.writeTF(self.getPrepend() + "_mmest", self.evalFreq, self.evalFreqEqns)
		# self.writeTF(self.getPrepend() + "_stacked", self.evalFreq, evalFreqEqnsTest)
		# self.writeTF(self.getPrepend() + "_reduced", self.evalFreq, evalFreqEqnsTest2)
		# self.writeTF(self.getPrepend() + "_ols", self.evalFreq, evalFreqEqnsTest3)
		self.writeTF(self.getPrepend() + "_cm", self.evalFreq, evalFreqEqnsTest4, variances=evalFreqVarsTest4)
		# self.writeTF(self.getPrepend() + "_cmMod", self.evalFreq, evalFreqEqnsTest5)
		# self.writeTF(self.getPrepend() + "_mmestStacked", self.evalFreq, evalFreqEqnsTest6)		

	###################
	### SOLVER ROUTINES
	###################
	def calcEvalFrequencyData(self, freq, evalFreq, winDataMatrix):
		# interpolate data to the evaluation frequencies
		inShape = winDataMatrix.shape
		data = np.empty(shape=(evalFreq.size, inShape[0], inShape[1]), dtype="complex")
		# get data from winDataMatrix
		for i in xrange(0, inShape[0]):
			for j in xrange(0, inShape[1]):
				interpFunc = interp.interp1d(freq, winDataMatrix[i,j])
				interpVals = interpFunc(evalFreq)
				for eIdx, eFreq in enumerate(evalFreq):
					data[eIdx,i,j] = interpVals[eIdx]
		return data

	def smoothSpectralEstimates(self, data):
		# takes the evaluation frequency data, which is indexed
		# windows, matrix of spectral components
		winSmooth = 9
		totalChans = self.getInSize() + self.getOutSize()
		for i in xrange(0, totalChans):
			for j in xrange(0, totalChans):
				data[:,i,j] = smooth1d(data[:,i,j], winSmooth, self.getWindow())
		return data

	def checkForBadValues(self, numWindows, data):
		finiteArray = np.ones(shape=(numWindows))
		for iW in xrange(0, numWindows):
			if not np.isfinite(data[iW]).all():
				finiteArray[iW] = 0
		numGoodWindows = sum(finiteArray)
		if numGoodWindows == numWindows:
			return numWindows, data
		self.printWarning("Bad data found...number of windows reduced from {} to {}".format(numWindows, numGoodWindows))
		goodWindowIndices = np.where(finiteArray == 1)
		return numGoodWindows, data[goodWindowIndices]

	def prepareLinearEqn(self, data):
		# prepare observations and regressors for linear processing
		numWindows = data.shape[0]
		numWindows, data = self.checkForBadValues(numWindows, data)
		totalSize = self.getOutSize() + self.getInSize()
		# for each output variable, have ninput regressor variables
		# let's construct our arrays
		obs = np.empty(shape=(self.getOutSize(), totalSize*numWindows), dtype="complex")
		reg = np.empty(shape=(self.getOutSize(), totalSize*numWindows, self.getInSize()), dtype="complex")
		for iW in xrange(0, numWindows):
			iOffset = iW*totalSize
			for i in xrange(0, self.getOutSize()):
				for j in xrange(0, totalSize):
					# this is the observation row where,i is the observed output
					obs[i, iOffset + j] = data[iW, self.getInSize() + i, j]
					for k in xrange(0, self.getInSize()):
						reg[i, iOffset + j, k] = data[iW, k, j]
		return numWindows, obs, reg

	def prepareLinearEqn2(self, data):
		# prepare observations and regressors for linear processing
		numWindows = data.shape[0]
		numWindows, data = self.checkForBadValues(numWindows, data)
		totalSize = self.getOutSize() + self.getInSize()
		crossSize = len(self.getCrossChannels())
		# for each output variable, have ninput regressor variables
		# let's construct our arrays
		obs = np.empty(shape=(self.getOutSize(), crossSize*numWindows), dtype="complex")
		reg = np.empty(shape=(self.getOutSize(), crossSize*numWindows, self.getInSize()), dtype="complex")
		for iW in xrange(0, numWindows):
			iOffset = iW*crossSize
			for i in xrange(0, self.getOutSize()):
				for j, crossChan in enumerate(self.crossChannels):
					# this is the observation row where,i is the observed output
					crossIndex = self.allChannels.index(crossChan)
					obs[i, iOffset + j] = data[iW, self.getInSize() + i, crossIndex]
					for k in xrange(0, self.getInSize()):
						reg[i, iOffset + j, k] = data[iW, k, crossIndex]
		return numWindows, obs, reg	

	# def calcConfidenceIntervals(self, obs, reg, out):
			

	def robustProcess(self, numWindows, obs, reg):
		# do mmestimate robust processing for a single evaluation frequency
		crossSize = len(self.getCrossChannels())
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		# solve
		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, scale, weights = mmestimateModel(predictors, observation, intercept=self.getIntercept())

			# now take the weights, apply to the observations and predictors, stack the appropriate rows and test
			observation2 = np.zeros(shape=(crossSize), dtype="complex")
			predictors2 = np.zeros(shape=(crossSize, self.getInSize()), dtype="complex")
			for iChan in xrange(0, crossSize):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*crossSize, crossSize)
				weightsLim = weights[indexArray]
				# weightsLim = weightsLim/np.sum(weightsLim) # normalise weights to 1
				observation2[iChan] = np.sum(obs[i, indexArray]*weightsLim)/numWindows
				# now for the regressors
				for j in xrange(0, self.getInSize()):
					predictors2[iChan, j] = np.sum(reg[i, indexArray, j]*weightsLim)/numWindows
			out, resids, scale, weights = mmestimateModel(predictors2, observation2, intercept=self.getIntercept())

			if self.getIntercept():
				output[i] = out[1:]
			else:
				output[i] = out
		return output

	def robustProcessStack(self, numWindows, obs, reg):
		# loop over the outputs
		output = np.zeros(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		totalSize = self.getOutSize() + self.getInSize()		
		for i in xrange(0, self.getOutSize()):
			# lets out some easier lettering
			y = obs[i]
			A = reg[i]
			# get some sizes
			n = A.shape[0]
			p = A.shape[1]		
			# first calculate the leverage weights
			# this is based on the hat matrix
			q, r = linalg.qr(A)			
			Pdiag = np.empty(shape=(n), dtype="float")
			for iRow in xrange(0, n):
				Pdiag[iRow] = np.absolute(np.sum(q[iRow,:]*np.conjugate(q[iRow,:]))).real
			del q, r
			Pdiag = Pdiag/np.max(Pdiag)
			leverageScale = sampleMAD0(Pdiag)
			leverageWeights = getRobustLocationWeights(Pdiag/leverageScale, "huber")	
		
			# Begin with stacking the data and solving
			observation = np.zeros(shape=(totalSize), dtype="complex")
			predictors = np.zeros(shape=(totalSize, self.getInSize()), dtype="complex")
			for iChan in xrange(0, totalSize):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*totalSize, totalSize)
				observation[iChan] = np.sum(y[indexArray])/numWindows
				# now for the regressors
				for j in xrange(0, self.getInSize()):
					predictors[iChan, j] = np.sum(A[indexArray, j])/numWindows
			initParams, residsStacked, scaleStacked, weightsStacked = mmestimateModel(predictors, observation, intercept=False)
			# calculate out scale and weights
			resids = y - np.dot(A, initParams)
			scale = sampleMAD0(resids)
			weights = getRobustLocationWeights(resids/scale, "huber")*leverageWeights

			# now get m-estimates and do the process again
			maxiter = 50
			iteration = 0
			while iteration < maxiter:

				# now stack with the weights and solve again
				observation = np.zeros(shape=(totalSize), dtype="complex")
				predictors = np.zeros(shape=(totalSize, self.getInSize()), dtype="complex")
				for iChan in xrange(0, totalSize):
					# now need to have my indexing array
					indexArray = np.arange(iChan, numWindows*totalSize, totalSize)
					weightsLim = weights[indexArray]
					# weightsLim = weightsLim/np.sum(weightsLim) # normalise weights to 1
					observation[iChan] = np.sum(y[indexArray]*weightsLim)/numWindows
					# now for the regressors
					for j in xrange(0, self.getInSize()):
						predictors[iChan, j] = np.sum(A[indexArray, j]*weightsLim)/numWindows
				paramsNew, residsStacked, scaleStacked, weightsStacked = mmestimateModel(predictors, observation)

				# now calculate residsNew etc
				residsNew = y - np.dot(A, paramsNew)

				if np.sum(np.absolute(residsNew)) < eps():
					# then return everything here
					params = paramsNew
					resids = residsNew					
					break 

				scale = sampleMAD0(residsNew)
				# standardise and calculate weights
				weightsNew = getRobustLocationWeights(residsNew/scale, "huber")*leverageWeights
				# increment iteration and save weightsNew
				iteration = iteration + 1
				weights = weightsNew
				params = paramsNew
				# check to see whether the change is smaller than the tolerance
				# use the R method of checking change in residuals (can check change in params)
				changeResids = linalg.norm(residsNew-resids)/linalg.norm(residsNew)
				if changeResids < eps():
					# update residuals
					resids = residsNew
					break
				# update residuals
				resids = residsNew

			# another go with tukey weights
			# return to original solution
			resids = y - np.dot(A, initParams)	
			weights = getRobustLocationWeights(resids/scale, "bisquare")*leverageWeights

			while iteration < maxiter:
				# now stack with the weights and solve again
				observation = np.zeros(shape=(totalSize), dtype="complex")
				predictors = np.zeros(shape=(totalSize, self.getInSize()), dtype="complex")
				for iChan in xrange(0, totalSize):
					# now need to have my indexing array
					indexArray = np.arange(iChan, numWindows*totalSize, totalSize)
					weightsLim = weights[indexArray]
					# weightsLim = weightsLim/np.sum(weightsLim) # normalise weights to 1
					observation[iChan] = np.sum(y[indexArray]*weightsLim)/numWindows
					# now for the regressors
					for j in xrange(0, self.getInSize()):
						predictors[iChan, j] = np.sum(A[indexArray, j]*weightsLim)/numWindows
				paramsNew, residsStacked, scaleStacked, weightsStacked = mmestimateModel(predictors, observation)

				# now calculate residsNew etc
				residsNew = y - np.dot(A, paramsNew)

				if np.sum(np.absolute(residsNew)) < eps():
					# then return everything here
					params = paramsNew
					resids = residsNew					
					break 

				scale = sampleMAD0(residsNew)
				# standardise and calculate weights
				weightsNew = getRobustLocationWeights(residsNew/scale, "bisquare")*leverageWeights
				# increment iteration and save weightsNew
				iteration = iteration + 1
				weights = weightsNew
				params = paramsNew
				# check to see whether the change is smaller than the tolerance
				# use the R method of checking change in residuals (can check change in params)
				changeResids = linalg.norm(residsNew-resids)/linalg.norm(residsNew)
				if changeResids < eps():
					# update residuals
					resids = residsNew
					break
				# update residuals
				resids = residsNew	

			output[i] = params

		return output		

	def robustProcessOLS(self, numWindows, obs, reg):
		# ordinary least squares solution
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		# solve
		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, squareResid, rank, s = olsModel(predictors, observation, intercept=self.getIntercept())
			if self.getIntercept():
				output[i] = out[1:]
			else:
				output[i] = out
		return output

	def robustProcessCM(self, numWindows, obs, reg):
		# do the robust processing for a single evaluation frequency
		crossSize = len(self.getCrossChannels())
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		varOutput = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="float")				
		# solve
		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, weights = chatterjeeMachler(predictors, observation, intercept=self.getIntercept())

			# now take the weights, apply to the observations and predictors, stack the appropriate rows and test
			observation2 = np.zeros(shape=(crossSize), dtype="complex")
			predictors2 = np.zeros(shape=(crossSize, self.getInSize()), dtype="complex")
			for iChan in xrange(0, crossSize):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*crossSize, crossSize)
				weightsLim = weights[indexArray]
				# weightsLim = weightsLim/np.sum(weightsLim) # normalise weights to 1
				observation2[iChan] = np.sum(obs[i, indexArray]*weightsLim)/numWindows
				# now for the regressors
				for j in xrange(0, self.getInSize()):
					predictors2[iChan, j] = np.sum(reg[i, indexArray, j]*weightsLim)/numWindows
			out, resids, weights = chatterjeeMachler(predictors2, observation2, intercept=self.getIntercept())

			# now calculate out the varainces - have the solution out, have the weights
			# recalculate out the residuals with the final solution
			# calculate standard deviation of residuals
			# and then use chatterjee machler formula to estimate variances
			# this needs work - better to use an empirical bootstrap method, but this will do for now
			resids = np.absolute(observation - np.dot(predictors, out))
			scale = sampleMAD0(resids) # some measure of standard deviation, rather than using the standard deviation
			residsVar = scale*scale
			varPred = np.dot(hermitianTranspose(predictors), weights*predictors)
			varPred = np.linalg.inv(varPred) # this is a pxp matrix
			varOut = 1.91472*residsVar*varPred
			print varOut
			varOut = np.diag(varOut).real # this should be a real number
			print varOut			

			if self.getIntercept():
				output[i] = out[1:]
				varOutput[i] = varOut[1:]
			else:
				output[i] = out
				varOutput[i] = varOut

		return output, varOutput

	def robustProcessCMMod(self, numWindows, obs, reg):
		# do the robust processing for a single evaluation frequency
		crossSize = len(self.getCrossChannels())
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		# solve
		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, weights = chatterjeeMachlerMod(predictors, observation, intercept=self.getIntercept())

			# now take the weights, apply to the observations and predictors, stack the appropriate rows and test
			observation2 = np.zeros(shape=(crossSize), dtype="complex")
			predictors2 = np.zeros(shape=(crossSize, self.getInSize()), dtype="complex")
			for iChan in xrange(0, crossSize):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*crossSize, crossSize)
				weightsLim = weights[indexArray]
				# weightsLim = weightsLim/np.sum(weightsLim) # normalise weights to 1
				observation2[iChan] = np.sum(obs[i, indexArray]*weightsLim)/numWindows
				# now for the regressors
				for j in xrange(0, self.getInSize()):
					predictors2[iChan, j] = np.sum(reg[i, indexArray, j]*weightsLim)/numWindows
			out, resids, weights = chatterjeeMachlerMod(predictors2, observation2, intercept=self.getIntercept())

			if self.getIntercept():
				output[i] = out[1:]
			else:
				output[i] = out
		return output

	def robustProcessCMHadi(self, numWindows, obs, reg):
		# chatterjee machler robust processing wth Hadi weights
		crossSize = len(self.getCrossChannels())
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		# solve
		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, weights = chatterjeeMachlerHadi(predictors, observation, intercept=self.getIntercept())

			# now take the weights, apply to the observations and predictors, stack the appropriate rows and test
			observation2 = np.zeros(shape=(crossSize), dtype="complex")
			predictors2 = np.zeros(shape=(crossSize, self.getInSize()), dtype="complex")
			for iChan in xrange(0, crossSize):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*crossSize, crossSize)
				weightsLim = weights[indexArray]
				# weightsLim = weightsLim/np.sum(weightsLim) # normalise weights to 1
				observation2[iChan] = np.sum(obs[i, indexArray]*weightsLim)/numWindows
				# now for the regressors
				for j in xrange(0, self.getInSize()):
					predictors2[iChan, j] = np.sum(reg[i, indexArray, j]*weightsLim)/numWindows
			out, resids, weights = chatterjeeMachler(predictors2, observation2, intercept=self.getIntercept())

			if self.getIntercept():
				output[i] = out[1:]
			else:
				output[i] = out
		return output

	def robustProcessReduced(self, data):
		# do the robust processing for a single evaluation frequency
		# this is an array with size numWindows, numIn
		# here only the input channels are used as the multiplying channels
		numWindows = data.shape[0]
		totalChans = self.getOutSize() + self.getInSize()
		# check for bad values
		numWindows, data = self.checkForBadValues(numWindows, data)
		# for each output variable, have ninput regressor variables
		# let's construct our arrays
		obs = np.empty(shape=(self.getOutSize(), self.getInSize()*numWindows), dtype="complex")
		reg = np.empty(shape=(self.getOutSize(), self.getInSize()*numWindows, self.getInSize()), dtype="complex")
		for iW in xrange(0, numWindows):
			iOffset = iW*self.getOutSize()
			for i in xrange(0, self.getOutSize()):
				for j in xrange(0, self.getInSize()):
					# this is the observation row where,i is the observed output
					obs[i, iOffset + j] = data[iW, self.getInSize() + i, j]
					for k in xrange(0, self.getInSize()):
						reg[i, iOffset + j, k] = data[iW, k, j]

		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")

		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, weights = chatterjeeMachlerMod(predictors, observation, intercept=self.getIntercept())
			if self.getIntercept():
				output[i] = out[1:]
			else:
				output[i] = out
		return output

	def stackedProcess(self, data):
		# then do various sums
		numWindows = data.shape[0]
		crossSize = len(self.getCrossChannels())
		# unweighted sum (i.e. normal solution)
		unWeightedSum = np.sum(data, axis=0)
		unWeightedSum = unWeightedSum/numWindows

		# for each output variable, have ninput regressor variables
		# let's construct our arrays
		obs = np.empty(shape=(self.getOutSize(), crossSize), dtype="complex")
		reg = np.empty(shape=(self.getOutSize(), crossSize, self.getInSize()), dtype="complex")
		for i in xrange(0, self.getOutSize()):
			for j, crossChan in enumerate(self.crossChannels):
				crossIndex = self.allChannels.index(crossChan)
				obs[i, j] = unWeightedSum[self.getInSize() + i, crossIndex]
				for k in xrange(0, self.getInSize()):
					reg[i, j, k] = unWeightedSum[k, crossIndex]

		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")

		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, scale, weights = mmestimateModel(predictors, observation, intercept=self.getIntercept())
			if self.getIntercept():
				output[i] = out[1:]
			else:
				output[i] = out
		return output

	###################
	### UTIL FUNCTIONS
	###################
	def writeTF(self, prepend, freq, data, **kwargs):
		inChans = self.getInChannels()
		outChans = self.getOutChannels()
		rootpath = self.proj.getTransDataPathSite(self.getOutSite())
		filename = "{}_fs{:d}_{}".format(self.getOutSite(), int(self.decParams.getSampleFreq()), prepend)
		datapath = os.path.join(rootpath, "{}".format(int(self.getWinSelector().getSampleFreq())))
		checkAndMakeDir(datapath)
		outfile = os.path.join(datapath, filename)
		outF = open(outfile, "w")
		# first write evalfreq
		outF.write("Evaluation frequencies\n")
		outF.write("{}\n".format(arrayToString(freq)))
		# now need to write out the other data
		for i in xrange(0, self.getOutSize()):
			for j in xrange(0, self.getInSize()):
				outF.write("IJ = {},{} = {}{}\n".format(i, j, outChans[i], inChans[j]))
				dataArray = np.empty(shape=(len(freq)), dtype="complex")
				for ifreq in xrange(0, len(freq)):
					dataArray[ifreq] = data[ifreq][i,j]
				outF.write("{}\n".format(arrayToString(dataArray)))

		# check for variances
		if "variances" in kwargs:
			# write out variances
			variances = kwargs["variances"]
			for i in xrange(0, self.getOutSize()):
				for j in xrange(0, self.getInSize()):
					outF.write("Var IJ = {},{} = {}{}\n".format(i, j, outChans[i], inChans[j]))
					varArray = np.empty(shape=(len(freq)), dtype="float")
					for ifreq in xrange(0, len(freq)):
						varArray[ifreq] = variances[ifreq][i,j]
					outF.write("{}\n".format(arrayToString(varArray)))				
		outF.close()

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("SINGLE SITE PROCESSOR INFO BEGIN")
		self.printText("####################")
		self.printText("In Site = {}".format(self.inSite))
		self.printText("In Channels = {}".format(self.getInChannels()))
		self.printText("Out Site = {}".format(self.outSite))
		self.printText("Out Channels = {}".format(self.getOutChannels()))
		self.printText("####################")
		self.printText("SINGLE SITE PROCESSOR INFO END")
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Single Site Processor Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Single Site Processor Warning", warnStr)
