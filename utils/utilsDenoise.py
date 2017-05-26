from sklearn.decomposition import FastICA, PCA
import pywt
# import mlpy.wavelet as wavelet

###########################
### WAVELETS
###########################
def waveletDenoise(data, stdThresh, fs):
	# this is for weak signal processing
	# anything away from a certain std deviation globally gets removed
	wav = 'bior6.8'
	mode = 'per'
	strong = {}
	for c in data:
		coeffs = pywt.wavedec(data[c], wav, mode)	
		# want to calculate global statistics
		# coeffs here is a list of coefficients
		allData = coeffs[0]
		for idx, coeff in enumerate(coeffs):
			if idx < 1:
				continue
			allData = np.concatenate((allData, coeff))
	
		# statistics
		absAll = np.absolute(allData)
		avg = np.average(absAll)
		mo = np.median(absAll)
		std = np.std(absAll)
		# thresh = avg + stdThresh*std
		thresh = mo + stdThresh*std
		allDataSize = allData.size

		# now threshold
		for idx, coeff in enumerate(coeffs):
			absCoeff = np.absolute(coeffs[idx])
			size = absCoeff.size
			for j in xrange(0, size):
				if absCoeff[j] < thresh:
					coeffs[idx][j] = 0
		
		# inverse trans
		strong[c] = pywt.waverec(coeffs, wav, mode)	
		data[c] = data[c] - strong[c]
		
		allDataProcessed = coeffs[0]
		for idx, coeff in enumerate(coeffs):
			if idx < 1:
				continue
			allDataProcessed = np.concatenate((allDataProcessed, coeff))		
		
		plt.figure()		
		#plt.scatter(np.arange(0, allDataSize), absAll)
		plt.scatter(np.arange(0, allDataSize), np.absolute(allDataProcessed))
		plt.plot(np.ones(shape=(allDataSize))*avg)
		plt.plot(np.ones(shape=(allDataSize))*thresh)		
	plt.show()		
	return data, strong	

	
def waveletDenoiseStrong(data, fs):
	# this is for weak signal processing
	# anything away from a certain std deviation globally gets removed
	wav = 'bior6.8'
	mode = 'per'
	out = {}
	for c in data:
		coeffs = pywt.wavedec(data[c], wav, mode)	
		# want to calculate global statistics
		# coeffs here is a list of coefficients
		allData = coeffs[0]
		for idx, coeff in enumerate(coeffs):
			if idx < 1:
				continue
			allData = np.concatenate((allData, coeff))
	
		# statistics
		absAll = np.absolute(allData)
		# now drop the lowest 20% of coefficients
		absAll = np.sort(absAll)
		allDataSize = allData.size
		threshIndex = int(allDataSize*0.5)
		thresh = absAll[threshIndex]
		# try using robust techniques
		# robustLocation, robustScale = mestimate(allData)
		# # calculate the MAD
		# madScale = sampleMedian(np.absolute(allData - robustLocation))
		# madScale = madScale/0.6745
		# thresh = robustLocation +  3*madScale
		
		# now threshold
		for idx, coeff in enumerate(coeffs):
			absCoeff = np.absolute(coeffs[idx])
			size = absCoeff.size
			for j in xrange(0, size):
				if absCoeff[j] < thresh:
					coeffs[idx][j] = 0
		
		# inverse trans
		out[c] = pywt.waverec(coeffs, wav, mode)	
		
		allDataProcessed = coeffs[0]
		for idx, coeff in enumerate(coeffs):
			if idx < 1:
				continue
			allDataProcessed = np.concatenate((allDataProcessed, coeff))		
		
		# plt.figure()		
		# #plt.scatter(np.arange(0, allDataSize), absAll)
		# plt.scatter(np.arange(0, allDataSize), np.absolute(allData))
		# plt.plot(np.ones(shape=(allDataSize))*thresh)		
	#plt.show()		
	return out	

def waveletDenoiseWeak(data, fs):
	# this is for weak signal processing
	# anything away from a certain std deviation globally gets removed
	wav = 'bior6.8'
	mode = 'per'
	out = {}
	for c in data:
		coeffs = pywt.wavedec(data[c], wav, mode)	
		# want to calculate global statistics
		# coeffs here is a list of coefficients
		allData = coeffs[0]
		for idx, coeff in enumerate(coeffs):
			if idx < 1:
				continue
			allData = np.concatenate((allData, coeff))
	
		# statistics
		absAll = np.absolute(allData)
		# now drop the lowest 20% of coefficients
		# absAll = np.sort(absAll)
		allDataSize = allData.size
		# threshIndex = int(allDataSize*0.6)
		# thresh = absAll[threshIndex]
		robustLocation, robustScale = mestimate(allData)
		# calculate the MAD
		madScale = sampleMedian(np.absolute(allData - robustLocation))
		madScale = madScale/0.6745
		thresh = robustLocation + 3*madScale

		# now threshold
		for idx, coeff in enumerate(coeffs):
			absCoeff = np.absolute(coeffs[idx])
			size = absCoeff.size
			for j in xrange(0, size):
				if absCoeff[j] < thresh:
					coeffs[idx][j] = 0
		
		# inverse trans
		out[c] = data[c] - pywt.waverec(coeffs, wav, mode)	
		
		allDataProcessed = coeffs[0]
		for idx, coeff in enumerate(coeffs):
			if idx < 1:
				continue
			allDataProcessed = np.concatenate((allDataProcessed, coeff))		
		
		# plt.figure()		
		# #plt.scatter(np.arange(0, allDataSize), absAll)
		# plt.scatter(np.arange(0, allDataSize), np.absolute(allData))
		# plt.plot(np.ones(shape=(allDataSize))*thresh)	
		# plt.plot(np.ones(shape=(allDataSize))*robustLocation)	
		# plt.plot(np.ones(shape=(allDataSize))*robustScale)				
	#plt.show()		
	return out		
	
def uwtDenoise(data):
	# this is for weak signal processing
	# anything away from a certain std deviation globally gets removed
	wav = 'd'
	k = 20
	strong = {}
	stdThresh = 0.05
	for c in data:
		length = data[c].size
		coeffs = wavelet.uwt(signal.detrend(data[c]), wav, k)	
		x = np.arange(length)
		plt.figure(figsize=(20,20))
		plt.subplot(2,2,1)
		plt.plot(x, data[c])
		plt.subplot(2,2,3)		
		for i in xrange(0, coeffs.shape[0]):
			plt.plot(x, coeffs[i,:], label=i)
		plt.legend()
		plt.title(c)
		# try some statistical thresholding
		#absCoeff = np.absolute(coeffs)
		avg = np.average(coeffs)
		std = np.std(coeffs)
		thresh = avg + std*stdThresh
		plt.plot(x, np.ones(shape=(length))*avg)
		plt.plot(x, np.ones(shape=(length))*thresh)

		for i in xrange(0, coeffs.shape[0]):
			for j in xrange(0, coeffs.shape[1]):
				if np.abs(coeffs[i,j]) < thresh:
					coeffs[i,j] = 0
		print coeffs.shape
		# for i in xrange(0, coeffs.shape[0]):
		# 	if i in [9]:
		# 		continue
		# 	coeffs[i,:] = coeffs[0,:]*0

		highamp = wavelet.iuwt(coeffs, wav, k)
		out = data[c] - highamp
		plt.subplot(2,2,1)
		plt.plot(x, out)
		plt.subplot(2,2,2)
		plt.plot(x, out)		

		# plot ffts
		plt.subplot(2,2,4)
		plt.plot(np.log10(np.absolute(forwardFFT(data[c]))))
		plt.plot(np.log10(np.absolute(forwardFFT(highamp))))
		plt.plot(np.log10(np.absolute(forwardFFT(out))))


		# want to calculate global statistics
		# allData = coeffs[0]
		# for idx, coeff in enumerate(coeffs):
		# 	if idx < 1:
		# 		continue
		# 	allData = np.concatenate((allData, coeff))
	
		# # statistics
		# absAll = np.absolute(allData)
		# avg = np.average(absAll)
		# mo = np.median(absAll)
		# std = np.std(absAll)
		# # thresh = avg + stdThresh*std
		# thresh = mo + stdThresh*std
		# allDataSize = allData.size

		# # now threshold
		# for idx, coeff in enumerate(coeffs):
		# 	absCoeff = np.absolute(coeffs[idx])
		# 	size = absCoeff.size
		# 	for j in xrange(0, size):
		# 		if absCoeff[j] < thresh:
		# 			coeffs[idx][j] = 0
		
		# # inverse trans
		# strong[c] = pywt.waverec(coeffs, wav, mode)	
		# data[c] = data[c] - strong[c]
		
		# allDataProcessed = coeffs[0]
		# for idx, coeff in enumerate(coeffs):
		# 	if idx < 1:
		# 		continue
		# 	allDataProcessed = np.concatenate((allDataProcessed, coeff))		
		
		# plt.figure()		
		# #plt.scatter(np.arange(0, allDataSize), absAll)
		# plt.scatter(np.arange(0, allDataSize), np.absolute(allDataProcessed))
		# plt.plot(np.ones(shape=(allDataSize))*avg)
		# plt.plot(np.ones(shape=(allDataSize))*thresh)		
	plt.show()		
	return data, strong		

# def mlpyDenoise(data, fs):
	# # this is for weak signal processing
	# # anything away from a certain std deviation globally gets removed
	# wav = 'd'
	# k = 20
	# strong = {}
	# stdThresh = 3
	# probFreq = [8.0, 16.0, 32.0, 64.0]
	# for c in data:
		# length = data[c].size
		# coeffs = wavelet.dwt(signal.detrend(data[c]), wav, k)	
		# x = np.arange(length)
		# plt.figure(figsize=(20,20))
		# plt.subplot(2,2,1)
		# plt.plot(x, data[c])
		# plt.subplot(2,2,3)		
		# plt.plot(coeffs)
		# plt.title(c)
		# # try some statistical thresholding
		# # The output data has the following form
		# # J = \log_2(n). 
		# # k = 0 ... (2^j)-1. The total number of levels is 
		# levels = int(np.log2(length))
		# levelSizes = [1]
		# coeffList = [coeffs[0]]
		# ranges = [[0,1]]
		# levelFreq = [0]
		# statsSubset = np.array(coeffs[0])
		# # add the first smoothing term
		# ic = 1
		# # now lets separate the levels
		# for il in xrange(0, levels):
			# # calculate size 2^j
			# size = np.power(2, il)
			# levelSizes.append(size)
			# # get subsection of array
			# coeffList.append(coeffs[ic:ic+size])
			# ranges.append([ic, ic+size])
			# levelFreq.append((size*fs)/(2*length))
			# if levelFreq[-1] not in probFreq:
				# print levelFreq[-1]
				# # add to stats subset
				# statsSubset = np.hstack((statsSubset, coeffList[-1]))
			# ic = ic+size

		# print statsSubset.size
		# avg = np.average(statsSubset)
		# std = np.std(statsSubset)
		# avgAbs = np.average(np.absolute(statsSubset))
		# thresh = avg + std*stdThresh
		# plt.plot(x, np.ones(shape=(length))*avg)
		# plt.plot(x, np.ones(shape=(length))*thresh)

		# print levelSizes
		# print levelFreq
		# print ranges
		# for r, freq in zip(ranges, levelFreq):
			# # threshold only in prob frequencies
			# if freq in probFreq:
				# for i in xrange(r[0], r[1]):
					# if np.abs(coeffs[i]) > thresh:
						# coeffs[i] = np.sign(coeffs[i])*avgAbs

		# # highamp = wavelet.idwt(coeffs, wav, k)
		# # out = data[c] - highamp
		# out = wavelet.idwt(coeffs, wav, k)
		# plt.subplot(2,2,2)
		# plt.plot(x, out)		

		# # plot ffts
		# plt.subplot(2,2,4)
		# plt.plot(np.log10(np.absolute(forwardFFT(data[c]))), label='input')
		# #plt.plot(np.log10(np.absolute(forwardFFT(highamp))), label='highamp')
		# plt.plot(np.log10(np.absolute(forwardFFT(out))), label='output')	
		# plt.legend()
	# plt.show()		
	# return data, strong		




def cwtDenoise(data):
	for c in data:
		widths = np.arange(1, 31, 1)
		cwtmatr = signal.cwt(data[c], signal.ricker, widths)
		plt.figure()
		plt.pcolor(np.log10(np.absolute(cwtmatr)))
		plt.colorbar()
	plt.show()
		
	
def harmonicICA(fs, length, data, freq):
	# frequency is a list of frequencies of harmnoic noise
	# create the arrays
	numFreq = len(freq)*2 # multiply by 2 for sine and cosine
	icaIn = np.zeros(shape=(numFreq + 1, length), dtype='float')
	# generate time array
	dt = 1/fs
	timeArray = np.arange(0, (length)*dt, dt)
	# generate the frequency sin and cosine curves
	# plt.figure()
	# splot = 1
	# for c in data:	
		# plt.subplot(4, 1, splot)
		# plt.plot(timeArray, data[c])
		# splot = splot + 1
		
	# plt.figure()
	for idx, f in enumerate(freq):
		sinCurve = np.sin(2*math.pi*f*timeArray)
		cosCurve = np.cos(2*math.pi*f*timeArray)
		icaIn[idx*2, :] = sinCurve
		icaIn[idx*2 + 1, :] = cosCurve
		# plot
		# plt.subplot(numFreq, 1, idx*2 + 1)
		# plt.plot(timeArray, icaIn[idx*2, :])
		# plt.subplot(numFreq, 1, idx*2 + 2)		
		# plt.plot(timeArray, icaIn[idx*2 + 1, :])
	# iterate over the channels
	ica = FastICA(n_components=3)
	#ica = FastICA()
	for c in data:
		plt.figure()
		# copy data into final array
		icaIn[-1,:] = data[c]
		test = np.transpose(icaIn)		
		# perform ICA
		output = ica.fit_transform(test)
		output = np.transpose(output)
		print output.shape
		plt.subplot(6,1,1)
		plt.plot(timeArray, data[c])
		plt.subplot(6,1,2)
		plt.plot(timeArray, output[0,:])		
		plt.subplot(6,1,3)
		plt.plot(timeArray, output[1,:])		
		plt.subplot(6,1,4)
		plt.plot(timeArray, output[2,:])	
		# plt.subplot(6,1,5)
		# plt.plot(timeArray, output[3,:])		
		# plt.subplot(6,1,6)
		# plt.plot(timeArray, output[4,:])			
	plt.show()


def waveletSpectraDenoise(specData):
	# this is for weak signal processing
	# anything away from a certain std deviation globally gets removed
	wav = "db4"
	k = 20
	strong = {}
	stdThresh = 0.05
	for c in specData:
		length = specData[c].size
		coeffs = pywt.dwt(specData[c], wav)
		print coeffs	
		x = np.arange(length)
		# plt.figure(figsize=(20,20))
		# plt.subplot(2,2,1)
		# plt.plot(x, data[c])
		# plt.subplot(2,2,3)		
		# for i in xrange(0, coeffs.shape[0]):
		# 	plt.plot(x, coeffs[i,:], label=i)
		# plt.legend()
		# plt.title(c)
		# # try some statistical thresholding
		# #absCoeff = np.absolute(coeffs)
		# avg = np.average(coeffs)
		# std = np.std(coeffs)
		# thresh = avg + std*stdThresh
		# plt.plot(x, np.ones(shape=(length))*avg)
		# plt.plot(x, np.ones(shape=(length))*thresh)

		# for i in xrange(0, coeffs.shape[0]):
		# 	for j in xrange(0, coeffs.shape[1]):
		# 		if np.abs(coeffs[i,j]) < thresh:
		# 			coeffs[i,j] = 0
		# print coeffs.shape
		# # for i in xrange(0, coeffs.shape[0]):
		# # 	if i in [9]:
		# # 		continue
		# # 	coeffs[i,:] = coeffs[0,:]*0

		# highamp = wavelet.iuwt(coeffs, wav, k)
		# out = data[c] - highamp
		# plt.subplot(2,2,1)
		# plt.plot(x, out)
		# plt.subplot(2,2,2)
		# plt.plot(x, out)		

		# # plot ffts
		# plt.subplot(2,2,4)
		# plt.plot(np.log10(np.absolute(forwardFFT(data[c]))))
		# plt.plot(np.log10(np.absolute(forwardFFT(highamp))))
		# plt.plot(np.log10(np.absolute(forwardFFT(out))))


		# want to calculate global statistics
		# allData = coeffs[0]
		# for idx, coeff in enumerate(coeffs):
		# 	if idx < 1:
		# 		continue
		# 	allData = np.concatenate((allData, coeff))
	
		# # statistics
		# absAll = np.absolute(allData)
		# avg = np.average(absAll)
		# mo = np.median(absAll)
		# std = np.std(absAll)
		# # thresh = avg + stdThresh*std
		# thresh = mo + stdThresh*std
		# allDataSize = allData.size

		# # now threshold
		# for idx, coeff in enumerate(coeffs):
		# 	absCoeff = np.absolute(coeffs[idx])
		# 	size = absCoeff.size
		# 	for j in xrange(0, size):
		# 		if absCoeff[j] < thresh:
		# 			coeffs[idx][j] = 0
		
		# # inverse trans
		# strong[c] = pywt.waverec(coeffs, wav, mode)	
		# data[c] = data[c] - strong[c]
		
		# allDataProcessed = coeffs[0]
		# for idx, coeff in enumerate(coeffs):
		# 	if idx < 1:
		# 		continue
		# 	allDataProcessed = np.concatenate((allDataProcessed, coeff))		
		
		# plt.figure()		
		# #plt.scatter(np.arange(0, allDataSize), absAll)
		# plt.scatter(np.arange(0, allDataSize), np.absolute(allDataProcessed))
		# plt.plot(np.ones(shape=(allDataSize))*avg)
		# plt.plot(np.ones(shape=(allDataSize))*thresh)		
		# plt.show()		
	return data, strong	