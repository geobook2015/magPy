# utility functions for plotting
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.dates import DateFormatter,DayLocator, AutoDateLocator, AutoDateFormatter
from datetime import datetime, timedelta
# utils
from utilsWindow import *
from utilsIO import *
		
##############
### OPTIONS
#############	
def defaultPlotOptions():
	default = {}
	default["colsPerRow"] = 2
	default["suptitle"] = "Suptitle"
	default["title"] = "Title"
	default["logx"] = False
	default["logxalt"] = False
	default["logy"] = False
	default["logyalt"] = False
	default["grid"] = True
	default["xlim"] = False
	default["xlimalt"] = False
	default["ylim"] = False
	default["ylimalt"] = False
	default["xlabel"] = ""
	default["ylabel"] = ""
	default["timeTicks"] = AutoDateLocator()
	default["timeFormat"] = default["timeTicks"]
	default["timeNum"] = 8
	default["cTitle"] = "Time"
	default["clim"] = [0.0, 1.0]	
	default["rowSize"] = 7
	default["colSize"] = 8
	default["notime"] = False
	default["hist"] = False
	default["gradients"] = False
	default["eIndex"] = -1
	return default

def defaultGrads():
	return [0.0001, 0.0003, 0.0005, 0.0007, 0.001, 0.003, 0.005, 0.007, 0.01, 0.03, 
		0.05, 0.07, 0.1, 0.3, 0.5, 0.7, 1,3,5,7,10,30,50,70,100,300,500,700,1000,3000,5000,7000,10000]	

##############
### FONT SIZES
#############	
def getPlotFonts():
	plotFonts = {}
	plotFonts["suptitle"] = 18
	plotFonts["title"] = 16
	plotFonts["axisLabel"] = 16
	plotFonts["axisTicks"] = 14
	plotFonts["legend"] = 12
	return plotFonts

def getPresentationFonts():
	plotFonts = {}
	plotFonts["suptitle"] = 22
	plotFonts["title"] = 20
	plotFonts["axisLabel"] = 20
	plotFonts["axisTicks"] = 20
	plotFonts["legend"] = 20
	return plotFonts

def getPaperFonts():
	plotFonts = {}
	plotFonts["suptitle"] = 18
	plotFonts["title"] = 17
	plotFonts["axisLabel"] = 16
	plotFonts["axisTicks"] = 15
	plotFonts["legend"] = 15
	return plotFonts	

##############
### COLORMAPS
#############	
def colorbar2dTime():
	return plt.cm.viridis

def colorbar2dOther():
	return plt.cm.rainbow
