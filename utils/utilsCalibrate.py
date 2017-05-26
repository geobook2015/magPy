# utility functions for calibration
import numpy as np
import math
import xml.etree.ElementTree as ET

def getKnownCalibrationFormats():
	calExt = ['TXT', 'RSP', 'RSPX']
	calFormats = ['metronix', 'rsp', 'rspx']
	return calExt, calFormats

def getCalName(format, ext, sensor, serial, chopper):
	if format == 'metronix':
		return metronixName(ext, sensor, serial, chopper)
	elif format == 'rsp':
		return rspName(ext, sensor, serial, chopper)
	elif format == 'rspx':
		return rspxName(ext, sensor, serial, chopper)
	else:
		return metronixName(ext, sensor, serial, chopper)


def metronixName(ext, sensor, serial, chopper):
	return '{}{}.{}'.format(sensor, serial, ext)

def rspName(ext, sensor, serial, chopper):
	board = 'HF'
	if chopper:
		board = 'LF'	
	sensorNum = int(sensor[3:5])
	return 'Metronix_Coil-----TYPE-{:03d}_{}-ID-{:06d}.{}'.format(sensorNum, board, serial, ext)

def rspxName(ext, sensor, serial, chopper):
	board = 'HF'
	if chopper:
		board = 'LF'
	sensorNum = int(sensor[3:5])
	return 'Metronix_Coil-----TYPE-{:03d}_{}-ID-{:06d}.{}'.format(sensorNum, board, serial, ext)


# these are functions for reading calibration files
# all data is returned with units
# F [Hz]
# Magnitude [mV/nT]
# Phase [radians]

def readCalFile(format, filepath, sensor, serial, chopper):
	if format == 'metronix':
		return metronixData(filepath, chopper)
	elif format == 'rsp':
		return rspData(filepath)
	elif format == 'rspx':
		return rspxData(filepath)
	else:
		# return a unit response
		return defaultCalibration()

# Metronix file format
# both chopper on
# units are F [Hz], Magnitude [V/nT*Hz], Phase [deg]
def metronixData(filepath, chopper):
	# no static gain - already included
	staticGain = 1
	# open file
	f = open(filepath,'r')
	lines = f.readlines()
	numLines = len(lines)
	f.close()
	# variables to save line numbers
	chopperOn = 0
	chopperOff = 0
	# find locations for chopperOn and chopperOff
	for il in xrange(0, numLines):
		# remove whitespace and new line characters
		lines[il] = lines[il].strip()
		if 'Chopper On' in lines[il]:
			chopperOn = il
		if 'Chopper Off' in lines[il]:
			chopperOff = il

	# get the part of the file required depending on chopper on or off
	dataLines = []
	dataOn = chopperOff
	if chopper:
		dataOn = chopperOn
	# get the data - starting from the next line
	il = dataOn + 1
	while (il < numLines and lines[il] != ''):
		# save line then increment
		dataLines.append(lines[il])
		il = il + 1
	
	# get the data as an array
	data = linesToArray(dataLines)
	# sort and extend
	data = sortCalData(data)
	data = extendCalData(data)
	# unit manipulation
	# change V/(nT*Hz) to mV/nT
	data[:,1] = data[:,1]*data[:,0]*1000
	#data[:,1] = data[:,1]*1000
	# change phase to radians
	data[:,2] = data[:,2]*(math.pi/180)		
	return data, staticGain	


# units are F [Hz], Magnitude [mv/nT], Phase [deg]
def rspData(filepath):
	# open file
	f = open(filepath,'r')
	lines = f.readlines()
	numLines = len(lines)
	f.close()

	staticGain = 1
	dataOn = 0
	for il in xrange(0, numLines):
		# remove whitespace and new line characters
		lines[il] = lines[il].strip()
		# find static gain value
		if "StaticGain" in lines[il]:
			staticGain = float(lines[il].split()[1])
		if "FREQUENCY" in lines[il]:
			dataOn = il
	dataLines = []
	il = dataOn + 2
	# get part of file desired
	while (il < numLines and lines[il] != ''):
		# save line then increment
		dataLines.append(lines[il])
		il = il + 1

	# get the data as an array
	data = linesToArray(dataLines)	
	# change phase to radians and apply static gain
	data[:,1] = data[:,1]*staticGain
	data[:,2] = data[:,2]*(math.pi/180)
	# sort and extend	
	data = sortCalData(data)
	data = extendCalData(data)
	return data, staticGain		

# units are F [Hz], Magnitude [mv/nT], Phase [deg]
def rspxData(filepath):
	# this is xml format - use EL tree
	tree = ET.parse(filepath)
	root = tree.getroot()
	# static gain
	staticGain = 1
	if root.find("StaticGain") is not None:
		staticGain = float(root.find("StaticGain").text)		
	# get the calibration data
	dataList = []
	for resp in root.findall("ResponseData"):
		dataList.append([float(resp.get("Frequency")), float(resp.get("Magnitude")), float(resp.get("Phase"))])
	# now create array
	data = np.array(dataList)
	# change phase to radians and apply static gain
	data[:,1] = data[:,1]*staticGain	
	data[:,2] = data[:,2]*(math.pi/180)
	# sort and extend	
	data = sortCalData(data)
	data = extendCalData(data)		
	return data, staticGain

# sort the calData
def sortCalData(data):
	# sort such that it goes from low frequency to high
	return data[data[:,0].argsort()]

# extend the calData - the data should already be sorted
def extendCalData(data):
	# add a line at the top (low frequency) extending the calibration information
	data = np.vstack((np.array([0.0000001, data[0,1], data[0,2]]), data))
	# add a line at the top (high frequency) extending the calibration information	
	data = np.vstack((data, np.array([100000000, data[-1,1], data[-1,2]])))
	return data	


def linesToArray(dataLines):
	# data to columns
	numData = len(dataLines) 
	for il in xrange(0, numData):
		dataLines[il] = dataLines[il].split()
	return np.array(dataLines, dtype=float)

# default calibration if none found
def defaultCalibration():
	return [1]*10, 1

def unitCalibration():
	unitCal = [
		[-100000000, 1, 0],
		[0, 1, 0], 
		[100000000, 1, 0]
	]
	return np.array(unitCal)



