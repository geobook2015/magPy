# utility functions for windowing
from datetime import datetime, timedelta
import calendar
import numpy as np
import math

def getDefaultMinWindowSize():
	return 512

def getDefaultMinOverlapSize():
	return 128

def gIndex2datetime(gIndex, refTime, fs, windowSize, windowOverlap):
	# global index 0 starts at refTime
	timeOffset = 1.0*(windowSize-windowOverlap)/fs
	totalOffset = gIndex*timeOffset
	startTime = refTime + timedelta(seconds=totalOffset)
	# windowSize - 1 because inclusive of start sample
	endTime = startTime + timedelta(seconds=(windowSize-1)/fs)
	return startTime, endTime

def gArray2datetime(gArray, refTime, fs, windowSize, windowOverlap):
	arrSize = gArray.size
	startTime = np.zeros(shape=(arrSize), dtype=datetime)
	endTime = np.zeros(shape=(arrSize), dtype=datetime)
	for i in xrange(0, arrSize):
		startTime[i], endTime[i] = gIndex2datetime(gArray[i], refTime, fs, windowSize, windowOverlap)
	return startTime, endTime	

def gIndex2timestamp(gIndex, refTime, fs, windowSize, windowOverlap):
	# global index 0 starts at refTime
	startTime, endTime = gIndex2datetime(gIndex, refTime, fs, windowSize, windowOverlap)
	return calendar.timegm(startTime.timetuple()), calendar.timegm(endTime.timetuple())

def gArray2timestamp(gArray, refTime, fs, windowSize, windowOverlap):
	arrSize = gArray.size
	startTime = np.zeros(shape=(arrSize), dtype=datetime)
	endTime = np.zeros(shape=(arrSize), dtype=datetime)
	for i in xrange(0, arrSize):
		startTime[i], endTime[i] = gIndex2timestamp(gArray[i], refTime, fs, windowSize, windowOverlap)
	return startTime, endTime	
	
def datetime2gIndex(refTime, inTime, fs, windowSize, windowOverlap):
	# need to return the next one close
	# calculate 
	deltaRefStart = inTime - refTime
	winStartIncrement = (windowSize-windowOverlap)/fs
	# calculate number of windows started before reference time
	# and then by taking the ceiling, find the global index of the first window in the data
	gIndex = int(math.ceil(deltaRefStart.total_seconds()/winStartIncrement))
	# calculate start time of first global window
	offsetSeconds = gIndex*winStartIncrement
	# calculate the first window time
	firstWindowTime = refTime + timedelta(seconds=offsetSeconds)
	return gIndex, firstWindowTime	

