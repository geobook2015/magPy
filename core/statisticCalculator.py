"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
Univariate Statistics
Controls saving and access to various statistics
The number of statistics per window has to be provided by the user
This is better for long-term usage
"""
import numpy as np
import scipy.stats as stats
import scipy.interpolate as interp
from copy import deepcopy
# import utils
from utilsIO import *
from utilsRobust import *
from utilsProcess import smooth1d

# frequency statistics
# have a single statistics for each evaluation frequency
# and for each window
# meaning, in total, there are nwindow*nfreq statistics
class StatisticCalculator(object):

	###################
	### CONSTRUCTOR
	##################
	def __init__(self):
		# default evaluation frequencies
		self.evalFreq = []
		# power smoothing vals
		self.winLen = 13
		self.winType = "hanning"		
		# set some defaults
		self.inChans = ["Hx", "Hy"]	
		self.inSize = len(self.inChans)
		self.outChans = ["Ex", "Ey"]
		self.outSize = len(self.outChans)
		self.specChans = self.inChans + self.outChans
		self.remoteChans = self.inChans
		self.psdChans = ["Ex", "Ey", "Hx", "Hy"]
		self.cohPairs = [["Ex", "Hx"], ["Ex", "Hy"], ["Ey", "Hx"], ["Ey", "Hy"]]
		self.polDirs = [["Ex", "Ey"],["Hx","Hy"]]	
		# set data presets
		self.spec = {}	
		# self.specChans = []
		# output data and marker for transfer function calculated
		self.tfCalculated = False
		self.remoteCalculated = False	
		self.intercept = False	
		self.outData = {}

	###################
	### GET FUNCTIONS
	###################
	def getEvalFreq(self):
		return deepcopy(self.evalFreq)

	def getInChans(self):
		return deepcopy(self.inChans)

	def getOutChans(self):
		return deepcopy(self.outChans)
		
	def getSpecChans(self):
		return deepcopy(self.specChans)

	def getRemoteChans(self):
		return deepcopy(self.remoteChans)

	def getPSDChans(self):
		return deepcopy(self.psdChans)

	def getCohPairs(self):
		return deepcopy(self.cohPairs)

	def getPolDirs(self):
		return deepcopy(self.polDirs)

	def getWinLen(self):
		return self.winLen

	def getWinType(self):
		return self.winType

	def getSpectra(self):
		return self.spec

	def getIntercept(self):
		return self.intercept

	# note: autopowers are real
	def getAutoPower(self, chan):
		idx = self.specChans.index(chan)
		# then return the autopower
		return self.spectralMatrix[idx, idx].real

	def getAutoPowerEval(self, chan, eIdx):
		idx = self.specChans.index(chan)
		# then return the autopower
		return self.evalMatrix[idx, idx, eIdx].real				

	def getCrossPower(self, chan1, chan2):
		idx1 = self.specChans.index(chan1)
		idx2 = self.specChans.index(chan2)
		# then return the autopower
		return self.spectralMatrix[idx1, idx2]

	def getCrossPowerEval(self, chan1, chan2, eIdx):
		idx1 = self.specChans.index(chan1)
		idx2 = self.specChans.index(chan2)
		# then return the autopower
		return self.evalMatrix[idx1, idx2, eIdx]	

	def getOutData(self):
		return deepcopy(self.outData)		

	###################
	### SET FUNCTIONS
	###################
	def setInChans(self, inChans):
		self.inChans = inChans
		self.inSize = len(self.inChans)

	def setOutChans(self, outChans):
		self.outChans = outChans
		self.outSize = len(self.outChans)

	def setRemoteChans(self, remoteChans):
		self.remoteChans = remoteChans

	def setPSDChans(self, psdChans):
		self.psdChans = psdChans

	def setCohPairs(self, cohPairs):
		self.cohPairs = cohPairs

	def setPolDirs(self, polDirs):
		self.polDirs = polDirs

	def setSpectra(self, freq, spectra, evalFreq):
		self.freq = freq
		self.spec = spectra
		self.evalFreq = evalFreq		
		# self.specChans = sorted(self.spec.keys())
		self.numChans = len(self.specChans)
		self.dataSize = self.spec[self.specChans[0]].size
		# calculate the power matrix
		self.calculateSpectralMatrix()
		self.calculateEvalMatrix()		
		# clear the out dictionary and set that transfer function not calculated
		self.prepareOutDict()

	def setIntercept(self, intercept):
		self.intercept = intercept

	###################
	### INITIAL HELPER FUNCTIONS
	### SPEED UP OTHER CALCULATIONS
	###################
	def calculateSpectralMatrix(self):
		# create the 3d array
		self.spectralMatrix = np.empty(shape=(self.numChans, self.numChans, self.dataSize), dtype="complex")
		# now need to go through the chans
		for i in xrange(0, self.numChans):
			for j in xrange(i, self.numChans):
				chan1 = self.specChans[i]
				chan2 = self.specChans[j]
				self.spectralMatrix[i, j] = smooth1d(self.spec[chan1]*np.conjugate(self.spec[chan2]), self.winLen, self.winType)
				if i == j:
					self.spectralMatrix[j, i] = np.conjugate(self.spectralMatrix[i, j]) # conjugate symmtry

	def calculateEvalMatrix(self):
		# create the array
		self.evalMatrix = np.empty(shape=(self.numChans, self.numChans, len(self.evalFreq)), dtype="complex")
		for i in xrange(0, self.numChans):
			for j in xrange(i, self.numChans):
				self.evalMatrix[i,j] = self.interpolateToEvalFreq(self.spectralMatrix[i,j])
				if i != j:
					self.evalMatrix[j, i] = np.conjugate(self.evalMatrix[i, j]) # conjugate symmtry 

	###################
	### ADD REMOTE SPEC
	### AND REMOTE GET FUNCTIONS
	###################
	def addRemoteSpec(self, remoteSpec, **kwargs):
		self.remoteSpec = remoteSpec
		if "remotechans" in kwargs:
			self.remoteChans = kwargs["remotechans"]
		# now calculate some remote reference related values
		self.calculateRemoteSpectralMatrix()
		self.calculateRemoteEvalMatrix()
		self.calculateReferenceSpectralMatrix()
		self.calculateReferenceEvalMatrix()

	def calculateRemoteSpectralMatrix(self):
		# create the 3d array
		numRemoteChans = len(self.remoteChans)
		self.remoteSpectralMatrix = np.empty(shape=(numRemoteChans, numRemoteChans, self.dataSize), dtype="complex")
		# now need to go through the chans
		for i in xrange(0, numRemoteChans):
			for j in xrange(i, numRemoteChans):
				chan1 = self.remoteChans[i]
				chan2 = self.remoteChans[j]
				self.remoteSpectralMatrix[i, j] = smooth1d(self.remoteSpec[chan1]*np.conjugate(self.remoteSpec[chan2]), self.winLen, self.winType)
				if i == j:
					self.remoteSpectralMatrix[j, i] = np.conjugate(self.remoteSpectralMatrix[i, j]) # conjugate symmtry

	def calculateRemoteEvalMatrix(self):
		# create the array
		numRemoteChans = len(self.remoteChans)
		self.remoteEvalMatrix = np.empty(shape=(numRemoteChans, numRemoteChans, len(self.evalFreq)), dtype="complex")
		for i in xrange(0, numRemoteChans):
			for j in xrange(i, numRemoteChans):
				self.remoteEvalMatrix[i,j] = self.interpolateToEvalFreq(self.remoteSpectralMatrix[i,j])
				if i != j:
					self.remoteEvalMatrix[j, i] = np.conjugate(self.remoteEvalMatrix[i, j]) # conjugate symmtry 

	def calculateReferenceSpectralMatrix(self):
		# cannot use conjugate symmetry in this case
		self.referenceSpectralMatrix = np.empty(shape=(self.numChans, len(self.remoteChans), self.dataSize), dtype="complex")
		for i, chan1 in enumerate(self.specChans):
			for j, chan2 in enumerate(self.remoteChans):
				self.referenceSpectralMatrix[i,j] = smooth1d(self.spec[chan1]*np.conjugate(self.remoteSpec[chan2]), self.winLen, self.winType)

	def calculateReferenceEvalMatrix(self):
		self.referenceEvalMatrix = np.empty(shape=(self.numChans, len(self.remoteChans), len(self.evalFreq)), dtype="complex")
		for i, chan1 in enumerate(self.specChans):
			for j, chan2 in enumerate(self.remoteChans):
				self.referenceEvalMatrix[i,j] = self.interpolateToEvalFreq(self.referenceSpectralMatrix[i,j])

	def getRemoteAutoPower(self, chan):
		idx = self.remoteChans.index(chan)
		return self.remoteSpectralMatrix[idx, idx].real

	def getRemoteAutoPowerEval(self, chan, eIdx):
		idx = self.remoteChans.index(chan)
		return self.remoteEvalMatrix[idx, idx, eIdx].real

	def getRemoteCrossPower(self, chan1, chan2):
		idx1 = self.remoteChans.index(chan1)
		idx2 = self.remoteChans.index(chan2)
		return self.remoteSpectralMatrix[idx1, idx2]

	def getRemoteCrossPowerEval(self, chan1, chan2, eIdx):
		idx1 = self.remoteChans.index(chan1)
		idx2 = self.remoteChans.index(chan2)
		return self.remoteSpectralMatrix[idx1, idx2, eIdx]		

	def getReferenceCrossPower(self, dataChan, remoteChan):
		idx1 = self.specChans.index(dataChan)
		idx2 = self.remoteChans.index(remoteChan)
		return self.referenceSpectralMatrix[idx1,idx2]

	def getReferenceCrossPowerEval(self, dataChan, remoteChan, eIdx):
		idx1 = self.specChans.index(dataChan)
		idx2 = self.remoteChans.index(remoteChan)
		return self.referenceEvalMatrix[idx1, idx2, eIdx]	

	###################
	### HELPER FUNCTION - dictionaries and interpolate to eval freq
	###################
	def interpolateToEvalFreq(self, data):
		interpFunc = interp.interp1d(self.freq, data)
		interpData = interpFunc(self.evalFreq)
		return interpData

	def prepareOutDict(self):
		self.outData = {}
		for e in self.evalFreq:
			self.outData[e] = {}
		# set various calculated flags to false
		self.tfCalculated = False
		self.remoteCalculated = False					

	###################
	### HELPER FUNCTION - return based on name of stat
	###################
	def getDataForStatName(self, statName):
		if statName == "absvalEqn":
			return self.winAbsVal()
		if statName == "psd":
			return self.winPSD()
		elif statName == "coherence":
			return self.winCoherence()
		elif statName == "poldir":
			return self.winPolarisations()
		elif statName == "partialcoh":
			return self.winPartials()
		elif statName == "transFunc" or statName == "resPhase":
			if self.tfCalculated:
				return self.getOutData()
			return self.winTransferFunction()
		elif statName == "coherenceRR":
			return self.winRemoteCoherence()
		elif statName == "coherenceRREqn":
			return self.winRemoteEqnCoherence()
		elif statName == "absvalRREqn":
			return self.winRemoteAbsVal()
		elif statName == "transFuncRR" or statName == "resPhaseRR":
			if self.remoteCalculated:
				return self.getOutData()
			return self.winRemoteTransferFunction()
		else:
			self.printWarning("Statistic in getDataForStatName not recognised")
			return self.winCoherence()

	###################
	### CALCULATE STATISTICS
	### POWER / COHERENCIES / POLARISATION DIRECTIONS
	###################
	def winPSD(self):
		# calculate PSD - want to divide by length of time too
		freqLen = self.freq.size
		timeLen = (freqLen-1)*2 # minus 1 because time sections are usually even
		fs = self.freq[-1]*2 # sampling frequency
		# and then calculate amount of time
		duration = timeLen/fs	
		# interpolate onto evaluation frequency and output to outData
		for eIdx, eF in enumerate(self.evalFreq):
			for chan in self.getPSDChans():
				key = "psd{}".format(chan)		
				self.outData[eF][key] = self.getAutoPowerEval(chan, eIdx)/duration	
		return self.getOutData()

	def winCoherence(self):
		# now calculate out the relevant coherencies	
		for idx, p in enumerate(self.getCohPairs()):
			c1 = p[0] # chan1
			c2 = p[1] # chan2
			for eIdx, eF in enumerate(self.evalFreq):				
				# now calculate the nominator and denominator	
				cohNom = np.power(np.absolute(self.getCrossPowerEval(c1, c2, eIdx)), 2).real
				cohDenom = self.getAutoPowerEval(c1, eIdx) * self.getAutoPowerEval(c2, eIdx)
				# save in outData
				key = "coh{}".format(c1+c2)		
				self.outData[eF][key] = cohNom/cohDenom
		return self.getOutData()

	def winPolarisations(self):
		# calculate polarisation directions
		for idx, p in enumerate(self.getPolDirs()):
			c1 = p[0] # chan1
			c2 = p[1] # chan2	
			for eIdx, eF in enumerate(self.evalFreq):				
				# now calculate the nominator and denominator	
				cohNom = 2*self.getCrossPowerEval(c1, c2, eIdx).real # take the real part of this
				cohDenom = self.getAutoPowerEval(c1, eIdx) - self.getAutoPowerEval(c2, eIdx)
				# save to out dictionary
				key = "pol{}".format(c1+c2)		
				self.outData[eF][key] = np.arctan(cohNom/cohDenom)*(180.0/np.pi)	
		return self.getOutData()	

	# this is based on paper Weckmann, Magunia Ritter 2005
	def winPartials(self):
		# calculate partial coherencies
		# e.g. Ex, Hx w.r.t Hy
		# this currently only works for impedance tensor calculations
		# do not want to get into higher power partial coherencies
		# get the coherences - these will be required later
		winCoherence = self.winCoherence()	

		for i, outChan in enumerate(self.outChans):
			for eIdx, eFreq in enumerate(self.evalFreq):
				inChan1 = self.inChans[0]
				inChan2 = self.inChans[1]
				xOutIn1 = self.getCrossPowerEval(outChan, inChan1, eIdx)
				xOutIn2 = self.getCrossPowerEval(outChan, inChan2, eIdx)
				xIn1In2 = self.getCrossPowerEval(inChan1, inChan2, eIdx)
				xIn2In1 = self.getCrossPowerEval(inChan2, inChan1, eIdx)
				# calculate out transFunc components
				denom = self.getAutoPowerEval(inChan1, eIdx)*self.getAutoPowerEval(inChan2, eIdx) - xIn1In2*xIn2In1
				# Z1
				Z1nom = xOutIn1*self.getAutoPowerEval(inChan2, eIdx) - xIn2In1*xOutIn2
				Z1 = Z1nom/denom
				# Z2
				Z2nom = self.getAutoPowerEval(inChan1, eIdx)*xOutIn2 - xIn1In2*xOutIn1
				Z2 = Z2nom/denom
				# calculate bivariate coherency
				rb = Z1*self.getCrossPowerEval(inChan1, outChan, eIdx) + Z2*self.getCrossPowerEval(inChan2, outChan, eIdx)
				rb = rb / self.getAutoPowerEval(outChan, eIdx)
				# now calculate out partials
				# calculate partial inChan, outChan1 with respect to outChan2
				cohkey = "coh{}".format(outChan+inChan2) 
				rp1 = (rb - winCoherence[eFreq][cohkey]) / (1.0 - winCoherence[eFreq][cohkey])
				# calculate partial inChan, outChan2 with respect to outChan1
				cohkey = "coh{}".format(outChan+inChan1)
				rp2 = (rb - winCoherence[eFreq][cohkey]) / (1.0 - winCoherence[eFreq][cohkey])
				# now save in outDict
				self.outData[eFreq]["bivar{}".format(outChan)] = rb
				self.outData[eFreq]["par{}".format(outChan+inChan1)] = rp1
				self.outData[eFreq]["par{}".format(outChan+inChan2)] = rp2
		return self.getOutData() 

	# simply save the absolute values of the cross power matrix
	# this is useful for cross plotting
	def winAbsVal(self):
		for eIdx, eFreq in enumerate(self.evalFreq):
			for iChan, chan in enumerate(self.specChans):
				# first do the outchans multiplied by every other channel
				for iOut, outChan in enumerate(self.outChans):			
					absval = np.absolute(self.getCrossPowerEval(outChan, chan, eIdx))
					key = "abs{}{}".format(outChan, chan)
					self.outData[eFreq][key] = absval
					
				# then  do the inchans multiplied by every other channel
				for iIn, inChan in enumerate(self.inChans):			
					absval = np.absolute(self.getCrossPowerEval(inChan, chan, eIdx))
					key = "abs{}{}".format(inChan, chan)
					self.outData[eFreq][key] = absval
		# return the dictionary
		return self.getOutData()
	
	###################
	### CALCULATE STATISTICS
	### TRANSFER FUNCTIONS
	###################
	# calculate components of impedance tensor 
	# separately for each window
	def winTransferFunction(self):
		# now do the work
		totalSize = self.inSize + self.outSize

		# now want to calculate the transfer function for each evaluation frequency
		output = np.empty(shape=(self.evalFreq.size, self.outSize, self.inSize), dtype="complex")		
		for eIdx, eFreq in enumerate(self.evalFreq):
			# solve transfer function
			obs = np.empty(shape=(self.outSize, totalSize), dtype="complex")
			reg = np.empty(shape=(self.outSize, totalSize, self.inSize), dtype="complex")
			for i in xrange(0, self.outSize):
				for j in xrange(0, totalSize):
					# this is the observation row where,i is the observed output
					# idx in the evaluation frequency
					obs[i, j] = self.getCrossPowerEval(self.outChans[i], self.specChans[j], eIdx)
					for k in xrange(0, self.inSize):
						reg[i, j, k] = self.getCrossPowerEval(self.inChans[k], self.specChans[j], eIdx)
			
			for i in xrange(0, self.outSize):
				observation = obs[i,:]
				predictors = reg[i,:,:]	
				# now do the solution
				out, resids, squareResid, rank, s = olsModel(predictors, observation, intercept=self.getIntercept())
				# out, resids, scale, weights	= mmestimateModel(predictors, observation, intercept=False)
				# not interested in the intercept (const) term
				if self.getIntercept():
					output[eIdx, i] = out[1:]
				else:
					output[eIdx, i] = out

			# calculate components of transfer function and res and phase    	
			for i in xrange(0, self.outSize):	
				for j in xrange(0, self.inSize):
					period = 1.0/eFreq
					res = 0.2 * period * np.power(np.absolute(output[eIdx, i, j]), 2)
					phase = np.angle(output[eIdx, i, j], deg=True)
					keyRes = self.outChans[i] + self.inChans[j] + "Res"
					keyPhase = self.outChans[i] + self.inChans[j] + "Phase"
					self.outData[eFreq][keyRes] = res
					self.outData[eFreq][keyPhase] = phase 
					# add the components
					keyReal = self.outChans[i] + self.inChans[j] + "Real"
					keyImag = self.outChans[i] + self.inChans[j] + "Imag"
					self.outData[eFreq][keyReal] = output[eIdx, i, j].real
					self.outData[eFreq][keyImag] = output[eIdx, i, j].imag	
		# set transfer function calculated as true
		# saves having to do it again
		self.tfCalculated = True
		return self.getOutData()

	###################
	### CALCULATE STATISTICS
	### REMOTE REFERENCE
	###################
	def winRemoteCoherence(self):
		# this is the coherence of ExHxR, ExHyR, EyHxR, EyHyR, HxHxR, HxHyR, HyHxR, HyHyR
		# now let's calculate coherency
		# abs(crosspower(A,B))^2/autopower(A)*autpower(B)
		for dataChan in self.specChans:
			for remoteChan in self.remoteChans:
				key = "{}{}RR".format(dataChan, remoteChan)				
				for eIdx, eFreq in enumerate(self.evalFreq):
					cohNom = np.power(np.absolute(self.getReferenceCrossPowerEval(dataChan, remoteChan, eIdx)),2)
					cohDenom = self.getAutoPowerEval(dataChan, eIdx)*self.getRemoteAutoPowerEval(remoteChan, eIdx)
					coh = cohNom/cohDenom
					self.outData[eFreq][key] = coh			
		return self.getOutData()

	def winRemoteEqnCoherence(self):
		# now calculate out the relevant coherencies
		# here we calculate the coherency between <Ex,HyR> and <Hy,HyR> for example
		for iOut, outChan in enumerate(self.outChans):
			for iIn, inChan in enumerate(self.inChans):
				for iRemote, remoteChan in enumerate(self.remoteChans):
					# calculate powers
					c1c1 = smooth1d(self.getReferenceCrossPower(outChan, remoteChan)*np.conjugate(self.getReferenceCrossPower(outChan, remoteChan)), self.winLen, self.winType) 
					c2c2 = smooth1d(self.getReferenceCrossPower(inChan, remoteChan)*np.conjugate(self.getReferenceCrossPower(inChan, remoteChan)), self.winLen, self.winType)
					c1c2 = smooth1d(self.getReferenceCrossPower(outChan, remoteChan)*np.conjugate(self.getReferenceCrossPower(inChan, remoteChan)), self.winLen, self.winType)
					# now interpolate
					c1c1 = self.interpolateToEvalFreq(c1c1)
					c2c2 = self.interpolateToEvalFreq(c2c2)
					c1c2 = self.interpolateToEvalFreq(c1c2)
					# now calculate the nominator and denominator	
					cohNom = np.power(np.absolute(c1c2), 2)
					cohDenom = c1c1*c2c2
					coh = cohNom/cohDenom # cast as float - discard complex part (complex part should be zero anyway)
					# now need the coherencies for the evaluation frequencies
					# this involves interpolation
					key = "{}{}R-{}{}R".format(outChan, remoteChan, inChan, remoteChan)				
					for iFreq, eFreq in enumerate(self.evalFreq):
						self.outData[eFreq][key] = coh[iFreq]
		return self.getOutData()	

	def winRemoteAbsVal(self):
		for eIdx, eFreq in enumerate(self.evalFreq):
			for iOut, outChan in enumerate(self.outChans):
				for iRemote, remoteChan in enumerate(self.remoteChans):
					absOut = np.absolute(self.getReferenceCrossPowerEval(outChan, remoteChan, eIdx))
					keyOut = "abs{}{}R".format(outChan, remoteChan)
					self.outData[eFreq][keyOut] = absOut
					for iIn, inChan in enumerate(self.inChans):
						absIn = np.absolute(self.getReferenceCrossPowerEval(inChan, remoteChan, eIdx))
						keyIn = "abs{}{}R".format(inChan, remoteChan)
						self.outData[eFreq][keyIn] = absIn
		return self.getOutData()

	def winRemoteTransferFunction(self):
		output = np.empty(shape=(self.evalFreq.size, self.outSize, self.inSize), dtype="complex")		
		for eIdx, eFreq in enumerate(self.evalFreq):
			# solve transfer function
			obs = np.empty(shape=(self.outSize, self.inSize), dtype="complex")
			reg = np.empty(shape=(self.outSize, self.inSize, self.inSize), dtype="complex")
			for i, outChan in enumerate(self.outChans):
				for j, remoteChan in enumerate(self.remoteChans):
					# this is the observation row where,i is the observed output
					# eIdx in the evaluation frequency
					obs[i, j] = self.getReferenceCrossPowerEval(outChan, remoteChan, eIdx)
					for k, inChan in enumerate(self.inChans):
						reg[i, j, k] = self.getReferenceCrossPowerEval(inChan, remoteChan, eIdx)
			
			for i in xrange(0, self.outSize):
				observation = obs[i,:]
				predictors = reg[i,:,:]	
				# now do the solution
				out, resids, squareResid, rank, s = olsModel(predictors, observation, intercept=self.getIntercept())
				# out, resids, scale, weights	= mmestimateModel(predictors, observation, intercept=False)				
				# not interested in the intercept (const) term
				if self.getIntercept():
					output[eIdx, i] = out[1:]
				else:
					output[eIdx, i] = out		

			# calculate components of transfer function and res and phase    	
			for i in xrange(0, self.outSize):	
				for j in xrange(0, self.inSize):
					period = 1.0/eFreq
					res = 0.2 * period * np.power(np.absolute(output[eIdx, i, j]), 2)
					phase = np.angle(output[eIdx, i, j], deg=True)
					keyRes = self.outChans[i] + self.inChans[j] + "ResRR"
					keyPhase = self.outChans[i] + self.inChans[j] + "PhaseRR"
					self.outData[eFreq][keyRes] = res
					self.outData[eFreq][keyPhase] = phase 
					# add the components
					keyReal = self.outChans[i] + self.inChans[j] + "RealRR"
					keyImag = self.outChans[i] + self.inChans[j] + "ImagRR"
					self.outData[eFreq][keyReal] = output[eIdx, i, j].real
					self.outData[eFreq][keyImag] = output[eIdx, i, j].imag	
		# set transfer function calculated as true
		# saves having to do it again
		self.remoteCalculated = True
		return self.getOutData()

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("STATISTIC CALCULATOR INFO BEGIN")		
		self.printText("####################")	
		self.printText("Default options")
		self.printText("\tInput Chans = {}".format(listToString(self.getInChans())))
		self.printText("\tOutput Chans = {}".format(listToString(self.getOutChans())))
		self.printText("\tRemote Chans = {}".format(listToString(self.getRemoteChans())))		
		self.printText("\tPowers = {}".format(listToString(self.getPSDChans())))
		self.printText("\tCoherence pairs = {}".format(listToString(self.getCohPairs())))
		self.printText("\tPartial coherence = {}".format(listToString(self.getPolDirs())))
		if len(self.getEvalFreq()) == 0:
			self.printText("Evaluation frequencies = {}")
		else:
			self.printText("Evaluation frequencies = {}".format(arrayToString(self.getEvalFreq())))
		self.printText("####################")
		self.printText("STATISTIC CALCULATOR INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Statistic Calculator Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Statistic Calculator Warning", warnStr)	

