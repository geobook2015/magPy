"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
The remote reference processor calculates different types of spectra
Inherits from single site processor
Just does remote reference computations
"""
import os
import numpy as np
import scipy.signal as signal
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# utils
from utilsFreq import *
from utilsIO import *
from utilsRobust import *
from utilsProcess import *
# import ProcessorSingleSite
from processorSingleSite import ProcessorSingleSite

class ProcessorRemoteReference(ProcessorSingleSite):

	###################
	### SET DEFAULTS
	##################
	def setDefaults(self):
		# inputs
		self.inSite = ""
		self.inChannels = []
		self.outSite = ""
		self.outChannels = []
		self.remoteSite = ""
		self.remoteChannels = []
		# evaluation frequency data
		self.evalFreq = []
		self.evalFreqEqns = []
		# smoothing options
		self.win = "hanning"
		self.winSmooth = -1
		# intercept options
		self.intercept = False
		# output filename
		self.prepend = ""

	###################
	### GET GENERAL INFO
	##################
	def getRemoteSite(self):
		return self.remoteSite

	def getRemoteChannels(self):
		return self.remoteChannels

	def getRemoteSize(self):
		return self.remoteSize

	###################
	### SET REMOTE REFERENCE
	##################
	def setRemote(self, remoteSite, remoteChannels):
		self.remoteSite = remoteSite
		self.remoteChannels = remoteChannels
		self.remoteSize = len(remoteChannels)
		self.printText("Remote reference set with site {} and channels {}".format(self.remoteSite, self.remoteChannels))

	###################
	### PROCESS - ONLY THIS FUNCTION IS DIFFERENT
	##################
	def process(self):
		# different types of solution
		evalFreqEqnsTest = []
		evalFreqEqnsTest2 = []
		evalFreqEqnsTest3 = []
		evalFreqEqnsTest4 = []
		evalFreqVarsTest4 = []
		evalFreqEqnsTest5 = []
		# for each decimation level
		# read in the shared windows from all sites
		# for each evaluation frequency, store the data from each window
		# and then at the end, perform robust processing
		numLevels = self.getDecParams().getNumLevels()
		inChans = self.getInChannels()
		outChans = self.getOutChannels()
		dataChans = inChans + outChans
		remoteChans = self.getRemoteChannels()
		for iDec in xrange(0, numLevels):
			# print out some info
			self.printText("Processing decimation level {}".format(iDec))

			fs = self.getWinSelector().getDecParams().getSampleFreqLevel(iDec)
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
			totalChans = self.getInSize() + self.getOutSize()
			numEvalFreq = len(evalFreq)
			dataSize = self.getWinSelector().getDataSize(iDec)
			freq = np.linspace(0, fs/2, dataSize)
			# get the window smoothing params
			smoothLen = self.getWindowSmooth(datasize=dataSize)

			# create the data array
			# for each evaluation frequency
			# keep the spectral power information for all windows
			evalFreqData = np.empty(shape=(numEvalFreq, numWindows, totalChans, self.getRemoteSize()), dtype="complex")

			# an array for the in and out channels fourier data
			winDataArray = np.empty(shape=(totalChans, dataSize), dtype="complex")
			# an array for the remote reference fourier data
			winRemoteArray = np.empty(shape=(self.getRemoteSize(), dataSize), dtype="complex")
			# an array for the power spectra data
			winSpectraMatrix = np.empty(shape=(totalChans, self.getRemoteSize(), dataSize), dtype="complex")

			# loop over shared windows
			localWin = 0
			global2local = {}
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

				# now get the remote reference data - assume this does not equal input or output
				remoteSF, remoteReader = self.getWinSelector().getSpecReaderForWindow(self.getRemoteSite(), iDec, iWin)
				remoteData = remoteReader.readBinaryWindowGlobal(iWin)

				# get data into the right part of the arrays
				for i in xrange(0, self.getInSize()):
					winDataArray[i] = inData[inChans[i]]
				for i in xrange(0, self.getOutSize()):
					winDataArray[self.getInSize() + i] = outData[outChans[i]]
				for i in xrange(0, self.getRemoteSize()):
					winRemoteArray[i] = remoteData[remoteChans[i]]

				# and now can fill the parts of the matrix
				# recall, smooth the power spectra
				for iD, dataChan in enumerate(dataChans):
					for iR, remoteChan in enumerate(remoteChans):
						# calculate each one, cannot use complex symmetry here
						# cannot use conjugate symmetry like with the single site processor
						winSpectraMatrix[iD,iR] = smooth1d(winDataArray[iD] * np.conjugate(winRemoteArray[iR]), smoothLen, self.getWindow())
				# after running through all windows, calculate evaluation frequencies
				# calculate frequency array
				evalFreqData[:, localWin] = self.calcEvalFrequencyData(freq, evalFreq, winSpectraMatrix)

				# increment local window
				localWin = localWin + 1

			# now all the data has been collected
			# for each evaluation frequency, do the robust processing
			# and get the evaluation frequency data
			evalFreqEqns = []
			for eIdx in xrange(0, numEvalFreq):
				self.printText("Processing evaluation frequency = {:.6f} [Hz], period = {:.6f} [s]".format(evalFreq[eIdx], 1.0/evalFreq[eIdx]))
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
				# use process reduced - only the input channels from the remote reference
				# print "Prepare linear equation"
				numSolveWindows, obs, reg = self.prepareLinearEqn(evalFreqData[eIdx, localWinIndices])
				# print "Robust process"
				self.evalFreqEqns.append(self.robustProcess(numSolveWindows, obs, reg))
				# print "Robust process stacking solve"				
				# evalFreqEqnsTest.append(self.robustProcessStack(numSolveWindows, obs, reg))
				# print "Robust OLS"
				# evalFreqEqnsTest2.append(self.robustProcessOLS(numSolveWindows, obs, reg))
				# print "Robust stacked"
				# evalFreqEqnsTest3.append(self.stackedProcess(evalFreqData[eIdx, localWinIndices]))
				# print "Robust CM"
				out, var = self.robustProcessCM(numSolveWindows, obs, reg)
				evalFreqEqnsTest4.append(out)
				evalFreqVarsTest4.append(var)
				# evalFreqEqnsTest4.append(self.robustProcessCM(numSolveWindows, obs, reg))
				# evalFreqEqnsTest5.append(self.robustProcessCMMod(numSolveWindows, obs, reg))

		# write out all the data
		self.writeTF(self.getPrepend() + "_mmest", self.evalFreq, self.evalFreqEqns)
		# self.writeTF(self.getPrepend() + "_mestStack", self.evalFreq, evalFreqEqnsTest)
		# self.writeTF(self.getPrepend() + "_ols", self.evalFreq, evalFreqEqnsTest2)
		# self.writeTF(self.getPrepend() + "_stacked", self.evalFreq, evalFreqEqnsTest3)
		self.writeTF(self.getPrepend() + "_cm", self.evalFreq, evalFreqEqnsTest4, variances=evalFreqVarsTest4)
		# self.writeTF(self.getPrepend() + "_cmMod", self.evalFreq, evalFreqEqnsTest5)

	###################
	### SOLVER ROUTINES
	###################
	def prepareLinearEqn(self, data):
		# prepare observations and regressors for linear processing
		numWindows = data.shape[0]
		numWindows, data = self.checkForBadValues(numWindows, data)
		# for each output variable, have ninput regressor variables
		# let's construct our arrays
		obs = np.empty(shape=(self.getOutSize(), self.getRemoteSize()*numWindows), dtype="complex")
		reg = np.empty(shape=(self.getOutSize(), self.getRemoteSize()*numWindows, self.getInSize()), dtype="complex")
		for iW in xrange(0, numWindows):
			iOffset = iW*self.getRemoteSize()
			for i in xrange(0, self.getOutSize()):
				for j in xrange(0, self.getRemoteSize()):
					# this is the observation row where,i is the observed output
					obs[i, iOffset + j] = data[iW, self.getInSize() + i, j]
					for k in xrange(0, self.getInSize()):
						reg[i, iOffset + j, k] = data[iW, k, j]
		return numWindows, obs, reg

	def robustProcessStack(self, numWindows, obs, reg):
		# loop over the outputs
		output = np.zeros(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
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
			observation = np.zeros(shape=(self.getRemoteSize()), dtype="complex")
			predictors = np.zeros(shape=(self.getRemoteSize(), self.getInSize()), dtype="complex")
			for iChan in xrange(0, self.getRemoteSize()):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*self.getRemoteSize(), self.getRemoteSize())
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
				observation = np.zeros(shape=(self.getRemoteSize()), dtype="complex")
				predictors = np.zeros(shape=(self.getRemoteSize(), self.getInSize()), dtype="complex")
				for iChan in xrange(0, self.getRemoteSize()):
					# now need to have my indexing array
					indexArray = np.arange(iChan, numWindows*self.getRemoteSize(), self.getRemoteSize())
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
				observation = np.zeros(shape=(self.getRemoteSize()), dtype="complex")
				predictors = np.zeros(shape=(self.getRemoteSize(), self.getInSize()), dtype="complex")
				for iChan in xrange(0, self.getRemoteSize()):
					# now need to have my indexing array
					indexArray = np.arange(iChan, numWindows*self.getRemoteSize(), self.getRemoteSize())
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

	def robustProcess(self, numWindows, obs, reg):
		# do the mmestimate robust processing for a single evaluation frequency
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")

		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, scale, weights = mmestimateModel(predictors, observation, intercept=self.getIntercept())

			# now take the weights, apply to the observations and predictors, stack the appropriate rows and test
			observation2 = np.zeros(shape=(self.getRemoteSize()), dtype="complex")
			predictors2 = np.zeros(shape=(self.getRemoteSize(), self.getInSize()), dtype="complex")
			for iChan in xrange(0, self.getRemoteSize()):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*self.getRemoteSize(), self.getRemoteSize())
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

	def robustProcessCM(self, numWindows, obs, reg):
		# do the chatterjeeMachlerMod robust processing for a single evaluation frequency
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
			observation2 = np.zeros(shape=(self.getRemoteSize()), dtype="complex")
			predictors2 = np.zeros(shape=(self.getRemoteSize(), self.getInSize()), dtype="complex")
			for iChan in xrange(0, self.getRemoteSize()):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*self.getRemoteSize(), self.getRemoteSize())
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
		# do the chatterjeeMachlerMod robust processing for a single evaluation frequency
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		# solve
		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, weights = chatterjeeMachlerMod(predictors, observation, intercept=self.getIntercept())

			# now take the weights, apply to the observations and predictors, stack the appropriate rows and test
			observation2 = np.zeros(shape=(self.getRemoteSize()), dtype="complex")
			predictors2 = np.zeros(shape=(self.getRemoteSize(), self.getInSize()), dtype="complex")
			for iChan in xrange(0, self.getRemoteSize()):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*self.getRemoteSize(), self.getRemoteSize())
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
		# do the robust processing for a single evaluation frequency
		# create array for output
		output = np.empty(shape=(self.getOutSize(), self.getInSize()), dtype="complex")
		# solve
		for i in xrange(0, self.getOutSize()):
			observation = obs[i,:]
			predictors = reg[i,:,:]
			# save the output
			out, resids, weights = chatterjeeMachlerHadi(predictors, observation, intercept=self.getIntercept())

			# now take the weights, apply to the observations and predictors, stack the appropriate rows and test
			observation2 = np.zeros(shape=(self.getRemoteSize()), dtype="complex")
			predictors2 = np.zeros(shape=(self.getRemoteSize(), self.getInSize()), dtype="complex")
			for iChan in xrange(0, self.getRemoteSize()):
				# now need to have my indexing array
				indexArray = np.arange(iChan, numWindows*self.getRemoteSize(), self.getRemoteSize())
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

	def robustProcessOLS(self, numWindows, obs, reg):
		# do the robust processing for a single evaluation frequency
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

	def stackedProcess(self, data):
		# then do various sums
		numWindows =  data.shape[0]
		numWindows, data = self.checkForBadValues(numWindows, data)
		# unweighted sum (i.e. normal solution)
		unWeightedSum = np.sum(data, axis=0)
		unWeightedSum = unWeightedSum/numWindows
		# for each output variable, have ninput regressor variables
		# let's construct our arrays
		obs = np.empty(shape=(self.getOutSize(), self.getRemoteSize()), dtype="complex")
		reg = np.empty(shape=(self.getOutSize(), self.getRemoteSize(), self.getInSize()), dtype="complex")
		for i in xrange(0, self.getOutSize()):
			for j in xrange(0, self.getRemoteSize()):
				obs[i, j] = unWeightedSum[self.getInSize() + i, j]
				for k in xrange(0, self.getInSize()):
					reg[i, j, k] = unWeightedSum[k, j]
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

	def checkRemote(self):
		check = True
		check = check and self.getRemoteSize() == self.getInSize()
		check = check and self.getRemoteChannels() == self.getInChannels()
		return check

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("REMOTE REFERENCE PROCESSOR INFO BEGIN")
		self.printText("####################")
		self.printText("In Site = {}".format(self.getInSite()))
		self.printText("In Channels = {}".format(self.getInChannels()))
		self.printText("Out Site = {}".format(self.getOutSite()))
		self.printText("Out Channels = {}".format(self.getOutChannels()))
		self.printText("Remote Site = {}".format(self.getRemoteSite()))
		self.printText("Remote Channels = {}".format(self.getRemoteChannels()))
		self.printText("####################")
		self.printText("REMOTE REFERENCE PROCESSOR INFO END")
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Remote Reference Processor Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Remote Reference Processor Warning", warnStr)
