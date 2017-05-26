"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
Univariate Statistics
Controls saving and access to various statistics
The number of statistics per window has to be provided by the user
This is better for long-term usage
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.dates import DateFormatter, DayLocator, AutoDateLocator, AutoDateFormatter
from datetime import datetime, timedelta
# utils
from utilsWindow import *
from utilsIO import *
from utilsPlotter import *

class Plotter(object):

	def __init__(self):
		self.initialiseDefaults()
		self.usePlotFonts() 

	def initialiseDefaults(self):
		self.default = defaultPlotOptions()

	def usePlotFonts(self):
		self.plotFonts = getPlotFonts()

	# bigger font size for presentations
	def usePresentationFonts(self):
		self.plotFonts = getPresentationFonts()				

	##############
	### GET FUNCTIONS
	#############
	def getDefaultOptions(self):
		return self.default

	def getPlotFonts(self):
		return self.plotFonts		

	##############
	### SET DEFAULTS
	#############
	def setDefaultOptions(self, default):
		self.default = default

	def setDefaultSingle(self, key, val):
		self.default[key] = val

	def setPlotFonts(self, fonts):
		self.plotFonts = fonts

	def setPlotFontsSingle(self, key, val):
		self.plotFonts[key] = value

	##############
	### DEAL WITH KEYWORDS
	#############	
	def parseKeywords(self, keywords):	
		options = {}
		for w in self.default:
			options[w] = self.default[w]
			if w in keywords:
				options[w] = keywords[w]
		return options

	##############
	### PLOT FORMATTING
	#############	
	def applyPlotOptions(self, ax, options, **kwargs):
		plotFonts = self.plotFonts
		if "plotFonts" in kwargs:
			plotsFonts = kwargs["plotFonts"]
		ax.set_xlabel(options["xlabel"], fontsize=plotFonts["axisLabel"])
		ax.set_ylabel(options["ylabel"], fontsize=plotFonts["axisLabel"])
		if options["logx"]:
			ax.set_xscale("log")
		if options["logy"]:
			ax.set_yscale("log")
		if options["xlim"]:
			plt.xlim(options["xlim"])
		if options["ylim"]:
			plt.ylim(options["ylim"])
		if options["grid"]:	
			ax.grid(True)

	def dateTicks(self, gIndices, dates, timeNum):
		numVals = len(gIndices)
		if timeNum >= numVals:
			timeNum = numVals
		plotIndices = []
		for i in xrange(0, timeNum):
			plotIndices.append(int(i*numVals*1.0/(timeNum-1)))
		plotIndices[-1] = numVals-1
		ticks = []
		tickLabels = []	
		for i in plotIndices:
			ticks.append(gIndices[i])
			tickLabels.append(dates[i].strftime('%m-%d %H:%M:%S'))	
		return ticks, tickLabels

	##############
	### PLOT FUNCTIONS
	#############
	def close(self):
		plt.close("all")

	##############
	### COLORMAPS
	##############	
	def cmap2dTime(self):
		# return plt.cm.rainbow
		return plt.cm.viridis

	def cmap2dOther(self):
		return plt.cm.rainbow

	def cmapDiscretise(cmap, N):
	    """Return a discrete colormap from the continuous colormap cmap.
	    
	        cmap: colormap instance, eg. cm.jet. 
	        N: number of colors.
	    
	    Example
	        x = resize(arange(100), (5,100))
	        djet = cmap_discretize(cm.jet, 5)
	        imshow(x, cmap=djet)
	    """
	    
	    if type(cmap) == str:
	        cmap = get_cmap(cmap)
	    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
	    colors_rgba = cmap(colors_i)
	    indices = np.linspace(0, 1., N+1)
	    cdict = {}
	    for ki,key in enumerate(('red','green','blue')):
	        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
	    # Return colormap object.
	    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

	def cmapMap(self, function, cmap):
	    """ Applies function (which should operate on vectors of shape 3:
	    [r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.
	    """
	    cdict = cmap._segmentdata
	    step_dict = {}
	    # Firt get the list of points where the segments start or end
	    for key in ('red','green','blue'):         step_dict[key] = map(lambda x: x[0], cdict[key])
	    step_list = sum(step_dict.values(), [])
	    step_list = array(list(set(step_list)))
	    # Then compute the LUT, and apply the function to the LUT
	    reduced_cmap = lambda step : array(cmap(step)[0:3])
	    old_LUT = array(map( reduced_cmap, step_list))
	    new_LUT = array(map( function, old_LUT))
	    # Now try to make a minimal segment definition of the new LUT
	    cdict = {}
	    for i,key in enumerate(('red','green','blue')):
	        this_cdict = {}
	        for j,step in enumerate(step_list):
	            if step in step_dict[key]:
	                this_cdict[step] = new_LUT[j,i]
	            elif new_LUT[j,i]!=old_LUT[j,i]:
	                this_cdict[step] = new_LUT[j,i]
	        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
	        colorvector.sort()
	        cdict[key] = colorvector

	    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

	def cmapTruncate(self, cmap, minval=0.0, maxval=1.0, n=100):
	    new_cmap = colors.LinearSegmentedColormap.from_list(
	        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
	        cmap(np.linspace(minval, maxval, n)))
	    return new_cmap

	def cmapFlip(self, items, ncol):
	    return itertools.chain(*[items[i::ncol] for i in range(ncol)])	

	###################
	### DEBUG
	##################
	def printInfo(self):
		self.printText("####################")
		self.printText("PLOTTER INFO BEGIN")		
		self.printText("####################")	
		self.printText("Default plot options = {}".format(self.default))	
		self.printText("Default plot fonts = {}".format(self.plotFonts))		
		self.printText("####################")
		self.printText("PLOTTER INFO END")		
		self.printText("####################")

	def printText(self, infoStr):
		generalPrint("Plotter Info", infoStr)

	def printWarning(self, warnStr):
		warningPrint("Plotter Warning", warnStr)


##############
### FREQUENCY PLOTTER
#############	
# class PlotterSpectra(Plotter):



##############
### STATISTIC PLOTTER
#############	
class PlotterStatistic(Plotter):

	def __init__(self):
		self.initialiseDefaults()
		self.usePlotFonts() 
		self.winSelector = False
		self.numBins = 150

	##################
	### DATA PARAMS
	### FOR CALCULATING DATE ARRAYS
	##################	
	def setDataInfo(self, reftime, iDec, fsDec, winParams, statWindows):
		self.reftime = reftime
		self.iDec = iDec
		self.fsDec = fsDec		
		self.winParams = winParams
		self.statWindows = statWindows
		self.windowSize = self.winParams.getWindowSize(iDec)
		self.windowOverlap = self.winParams.getOverlap(iDec)		
		# stat window data is an array of global indices for the data
		self.dataDatesList = [] # dates for each dataset
		self.dataDates = []
		self.dataIndices = []
		for gIndices in self.statWindows:
			startdates, enddates = gArray2datetime(gIndices, self.reftime, self.fsDec, self.windowSize, self.windowOverlap)
			self.dataDates = self.dataDates + list(startdates)
			self.dataIndices = self.dataIndices + list(gIndices)	
			# save the dates in a list for each dataset
			self.dataDatesList.append(startdates)
		self.minDate = min(self.dataDates)
		self.maxDate = max(self.dataDates)	
		self.minIndex = min(self.dataIndices)
		self.maxIndex = max(self.dataIndices)
		# now create the plot indices
		# have these separately to support non consecutive global indices
		self.plotIndices = np.arange(self.minIndex, self.maxIndex + 1)
		self.plotDatesStart, self.plotDatesEnd = gArray2datetime(self.plotIndices, self.reftime, self.fsDec, self.windowSize, self.windowOverlap)	

	def addWinSelector(self, winSelector):
		self.winSelector = winSelector

	##################
	### PLOTS
	##################	
	def selectWindows(self, data, eIdx):
		# limit data by selecting windows
		if not self.winSelector:
			# then return the data and the dataIndices and dataDates belonging to the class
			return data, self.dataIndices, self.dataDates

		# otherwise get the windows for this evaluation frequency
		windowsToUse = list(self.winSelector.getWindowsForFreq(self.iDec, eIdx))
		dataSelected = []
		indicesSelected = []
		datesSelected = []
		for idx, d in enumerate(data):
			# now need to find the data indices of the limited data
			# inSet is a boolean array with Trues where the indices are in windowsToUse
			inSet = np.in1d(self.statWindows[idx], windowsToUse)
			if not np.any(inSet):
				continue # no data

			# now select the data
			dataSelected = dataSelected + list(d[inSet])
			indicesSelected = indicesSelected + list(self.statWindows[idx][inSet])
			datesSelected = datesSelected + list(self.dataDatesList[idx][inSet])
		return dataSelected, indicesSelected, datesSelected

	def scatter1d(self, fig, ax, data, options):
		data, dataIndices, dataDates = self.selectWindows(data, options["eIndex"])

		if options["notime"]:
			xarr = np.arange(0, len(dataIndices))
			scat = plt.scatter(xarr, data, c=dataIndices, edgecolors='none', marker='o', s=12, cmap=self.cmap2dTime())
		else:
			scat = plt.scatter(dataDates, data, c=dataIndices, edgecolors='none', marker='o', s=12, cmap=self.cmap2dTime())
		scat.set_clim([self.minIndex, self.maxIndex])
		# format the axes
		self.applyPlotOptions(ax, options)
		if not options["xlim"] and not options["notime"]:
			plt.xlim([self.minDate, self.maxDate])
			ax.format_xdata = DateFormatter("%H-%M-%S")
			# fig.autofmt_xdate()	# date formats	
		if options["notime"]:
			plt.xlim([xarr[0], xarr[-1]])
		return scat

	def histogram(self, fig, ax, data, options):
		histData, dataIndices, dataDates = self.selectWindows(data, options["eIndex"])
		if len(histData) == 0:
			return # need some way of dealing with this when there is no histData
		# get the bounds of the histogram data
		xlim = [np.min(histData), np.max(histData)]
		if options["ylim"]:
			xlim = options["ylim"]
		if options["logy"]:
			histData = np.log10(histData)
			xlim = np.log10(xlim)
			# remove nans and infinities
			histData = histData[np.isfinite(histData)]
		n, bins, patches = plt.hist(histData, self.numBins, range=xlim, facecolor="red", alpha=0.75)
		plt.xlim(xlim)
		plt.xlabel(options["ylabel"])
		plt.ylabel("Count")
		plt.grid(True)	

	def addColorbar(self, plot, cax, options, plotFonts):
		ticks, tickLabels = self.dateTicks(self.plotIndices, self.plotDatesStart, options["timeNum"])
		cb = plt.colorbar(plot, cax=cax)
		cb.set_ticks(ticks)
		cb.set_ticklabels(tickLabels)	
		cax.set_title(options["cTitle"], y=1.02, fontsize=plotFonts["title"])	

	def plotData1d(self, winStats, statData, **kwargs):
		# plot options	
		options = self.parseKeywords(kwargs)						
		plotFonts = self.getPlotFonts()

		# calculate the number of rows and columns
		numRows = int(np.ceil(1.0*len(winStats)/options["colsPerRow"])) # round up
		numCols = options["colsPerRow"]
		if options["hist"]:
			numRows = numRows*2
		fig = plt.figure(figsize=(options["colSize"]*numCols, options["rowSize"]*numRows))
		st = plt.suptitle(options["title"], fontsize=plotFonts["suptitle"])
		st.set_y(0.98)	
		for idx, ws in enumerate(winStats):	
			# collect all data		
			allY = []
			for dataIdx, data in enumerate(statData):	
				if options["eIndex"] >= 0:
					data = np.squeeze(data[:,options["eIndex"],:]) # otherwise, window statistic						
				allY.append(data[:,idx])		
			# plot
			ax = plt.subplot(numRows, numCols, idx+1)
			ax.set_title(ws, fontsize=plotFonts["title"])			
			scat = self.scatter1d(fig, ax, allY, options)			

			if options["hist"]:
				ax = plt.subplot(numRows, numCols, len(winStats) + idx + 1)
				ax.set_title("{} histogram".format(ws), fontsize=plotFonts["title"])
				self.histogram(fig, ax, allY, options)

		# layout
		fig.tight_layout()
		# shift subplots down:					
		fig.subplots_adjust(top=0.90, right=0.85)
		# do the overall colorbar
		cax = fig.add_axes([0.88, 0.07, 0.03, 0.83])		
		self.addColorbar(scat, cax, options, plotFonts)
		return fig

	def plotData1d_2comp(self, winStats, statData, **kwargs):
		# plot options	
		options = self.parseKeywords(kwargs)										
		altoptions = self.parseKeywords(kwargs)	
		altoptions["logy"] = altoptions["logyalt"]
		altoptions["ylim"] = altoptions["ylimalt"]
		plotFonts = self.getPlotFonts()

		# calculate the number of rows and columns - hard set this
		numRows = 2
		numCols = len(winStats)/2
		if options["hist"]:
			numRows = numRows*2		
		fig = plt.figure(figsize=(options["colSize"]*numCols, options["rowSize"]*numRows))	
		st = plt.suptitle(options["title"], fontsize=plotFonts["suptitle"])
		st.set_y(0.98)
		# loop over window stats and plot
		for idx in np.arange(0, len(winStats), 2):	
			# collect all data				
			allY1 = []
			allY2 = []	
			for dataIdx, data in enumerate(statData):
				if options["eIndex"] >= 0:
					data = np.squeeze(data[:,options["eIndex"],:]) # otherwise, window statistic								
				allY1.append(data[:,idx])
				allY2.append(data[:,idx+1])		

			# now plot
			ax = plt.subplot(numRows, numCols, (idx/2)+1)
			ax.set_title(winStats[idx], fontsize=plotFonts["title"])
			scat = self.scatter1d(fig, ax, allY1, options)				
			# plot the other part		
			ax = plt.subplot(numRows, numCols, (idx/2)+1+numCols)
			ax.set_title(winStats[idx+1], fontsize=plotFonts["title"])	
			scat = self.scatter1d(fig, ax, allY2, altoptions)	

			# histograms
			if options["hist"]:
				ax = plt.subplot(numRows, numCols, len(winStats)+(idx/2)+1)
				ax.set_title("{} histogram".format(winStats[idx]), fontsize=plotFonts["title"])
				self.histogram(fig, ax, allY1, options)		
				# and the other part
				ax = plt.subplot(numRows, numCols, len(winStats)+(idx/2)+1+numCols)				
				ax.set_title("{} histogram".format(winStats[idx+1]), fontsize=plotFonts["title"])
				self.histogram(fig, ax, allY2, altoptions)					
		
		# layout					
		fig.tight_layout()
		# shift subplots down:		
		fig.subplots_adjust(top=0.92, right=0.85)	
		# do the overall colorbar
		cax = fig.add_axes([0.88, 0.07, 0.03, 0.83])		
		self.addColorbar(scat, cax, options, plotFonts)
		return fig		

	def plotData2d(self, winStats, statData, **kwargs):
		# plot options	
		options = self.parseKeywords(kwargs)						
		plotFonts = self.getPlotFonts()

		# calculate the number of rows and columns
		numRows = int(np.ceil(0.5*len(winStats)/options["colsPerRow"])) # round up
		numCols = options["colsPerRow"]
		fig = plt.figure(figsize=(options["colSize"]*numCols, options["rowSize"]*numRows))	
		st = plt.suptitle(options["title"], fontsize=plotFonts["suptitle"])
		st.set_y(0.98)
		# loop over window stats and plot
		for idx in np.arange(0, len(winStats), 2):	
			allX = []
			allY = []		
			for data in statData:
				if options["eIndex"] >= 0:
					data = np.squeeze(data[:,options["eIndex"],:]) # otherwise, window statistic			
				allX.append(data[:,idx])
				allY.append(data[:,idx+1])

			# select data if given (indicesY and indiceX should be the same. the dates too)
			allX, indicesX, datesX = self.selectWindows(allX, options["eIndex"])
			allY, indicesY, datesY = self.selectWindows(allY, options["eIndex"])

			# now plot
			ax = plt.subplot(numRows, numCols, (idx/2)+1)
			ws = "{} and {}".format(winStats[idx], winStats[idx+1])	
			ax.set_title(ws, fontsize=plotFonts["title"])			
			scat = plt.scatter(allX, allY, c=indicesY, edgecolors='none', marker='o', s=12, cmap=self.cmap2dTime())
			scat.set_clim([self.minIndex, self.maxIndex])
			if len(allX) > 0 and options["gradients"]:
				# draw a set of gradients on the scatter plot
				minX = np.min(allX)
				maxX = np.max(allX)
				maxY = np.max(allY)
				for idx, g in enumerate(defaultGrads()):
					gminY = minX*g
					if gminY > maxY:
						continue
					gmaxY = maxX*g
					plotEndX = maxX
					if gmaxY > maxY:
						plotEndX = 1.0*maxY/g
					# don't plot if too close to the axes
					if plotEndX < 0.07*maxX:
						continue
					if plotEndX*g < 0.07*maxY:
						continue
					# otherwise plot
					plt.plot([minX, plotEndX], [minX*g, plotEndX*g], color="lightgray")
					if idx%3 == 0:
						plt.text(plotEndX-0.1*plotEndX, g*(plotEndX-0.075*plotEndX), "m = {}".format(g))
					elif idx%3 == 1:
						plt.text(plotEndX-0.1*plotEndX, g*(plotEndX-0.125*plotEndX), "m = {}".format(g))
					else:
						plt.text(plotEndX-0.1*plotEndX, g*(plotEndX-0.175*plotEndX), "m = {}".format(g))
			self.applyPlotOptions(ax, options)		

		# layout					
		fig.tight_layout()
		# shift subplots down:		
		fig.subplots_adjust(top=0.92, right=0.85)	
		# do the overall colorbar
		cax = fig.add_axes([0.88, 0.07, 0.03, 0.83])		
		self.addColorbar(scat, cax, options, plotFonts)
		return fig

	def plotData2d_1d(self, winStats, statData1, statData2, **kwargs):
		# plot options	
		options = self.parseKeywords(kwargs)						
		plotFonts = self.getPlotFonts()

		# calculate the number of rows and columns
		numRows = int(np.ceil(0.5*len(winStats)/options["colsPerRow"])) # round up
		numCols = options["colsPerRow"]
		fig = plt.figure(figsize=(options["colSize"]*numCols, options["rowSize"]*numRows))	
		st = plt.suptitle(options["title"], fontsize=plotFonts["suptitle"])
		st.set_y(0.98)
		# loop over window stats and plot
		for idx in np.arange(0, len(winStats), 2):	
			allX = []
			allY = []
			all2 = []		
			for data1, data2 in zip(statData1, statData2):
				if options["eIndex"] >= 0:
					data1 = np.squeeze(data1[:,options["eIndex"],:]) # otherwise, window statistic	
					data2 = np.squeeze(data2[:,options["eIndex"],:])
				allX.append(data1[:,idx])
				allY.append(data1[:,idx+1])
				all2.append(data2[:,idx/2])		

			# select data if given
			allX, indicesX, datesX = self.selectWindows(allX, options["eIndex"])
			allY, indicesY, datesY = self.selectWindows(allY, options["eIndex"])
			all2, indices2, dates2 = self.selectWindows(all2, options["eIndex"])

			# now plot
			ax = plt.subplot(numRows, numCols, (idx/2)+1)
			ws = "{} and {}".format(winStats[idx], winStats[idx+1])	
			ax.set_title(ws, fontsize=plotFonts["title"])			
			scat = plt.scatter(allX, allY, c=all2, edgecolors='none', marker='o', s=12, cmap=self.cmap2dOther())
			scat.set_clim(options["clim"])
			self.applyPlotOptions(ax, options)			

		# layout
		fig.tight_layout()
		# shift subplots down:					
		fig.subplots_adjust(top=0.92, right=0.85)
		# do the overall colorbar
		cax = fig.add_axes([0.88, 0.04, 0.03, 0.88])
		cb = plt.colorbar(scat , cax=cax)	
		cax.set_title(options["cTitle"], y=1.02, fontsize=plotFonts["title"])	
		return fig	

