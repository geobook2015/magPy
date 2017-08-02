# utility functions for readers
import os
import sys
import itertools
import glob
from datetime import datetime

###################
### DIRECTORY CHECKING
###################
# move this out in the future
def getDataDirectoryFormats():
	return ["meas", "run", "phnx"]

# get both files and directories in directory
def getDirectoryContents(path):
	if not checkDirExistence(path):
		return [], [] # return empty lists
	dirList = os.listdir(path) # current directory
	dirs = []
	files = []
	for d in dirList:
		if os.path.isdir(os.path.join(path,d)):
			dirs.append(d)
		else:
			files.append(d)
	return dirs, removeHiddenFiles(files)

# get files in directory
def getFilesInDirectory(path):
	dirs, files = getDirectoryContents(path)
	return files

# get directories in directory
def getDirsInDirectory(path):
	dirs, files = getDirectoryContents(path)
	return dirs

def removeHiddenFiles(files):
	filesNew = []
	for f in files:
		if f[0] != '.':
			filesNew.append(f)
	return(filesNew)

def getDataDirsInDirectory(path):
	dirs = getDirsInDirectory(path)
	dirsData = []
	formats = getDataDirectoryFormats()
	for d in dirs:
		for format in formats:
			if format in d:
				dirsData.append(d)
	return dirsData

# checks before writing
def checkDirExistence(path):
	if not os.path.exists(path):
		return False
	return True

def makeDir(path):
	os.makedirs(path)

def checkAndMakeDir(path):
	if not checkDirExistence(path):
		makeDir(path)

# put this into a function because it allows dealing with error messages in one place
def checkFilepath(path):
	if not os.path.exists(path):
		generalPrint("Utils IO", "File path {} could not be found.".format(path))
		return False
	return True

###################
### PRINTING FUNCTIONS
###################
# printing functions
# this should help tranisition to python 3
def generalPrint(pre, info):
	print "{} {}: {}".format(datetime.now().strftime("%H:%M:%S"), pre, info)

def warningPrint(pre, info):
	print "{} {}: {}".format(datetime.now().strftime("%H:%M:%S"), pre, info)

# useful output functions
def arrayToString(data):
	outputStr = ''
	for d in data:
		outputStr = outputStr + "{:.8f}, ".format(d)
	outputStr = outputStr.strip()
	outputStr = outputStr.rstrip(",")
	return outputStr

def arrayToStringTab(data):
	outputStr = ''
	for d in data:
		outputStr = outputStr + "{:.8f}\t".format(d)
	outputStr = outputStr.strip()
	return outputStr

def arrayToStringSci(data):
	outputStr = ''
	for d in data:
		outputStr = outputStr + "{:.6e}, ".format(d)
	outputStr = outputStr.strip()
	return outputStr.rstrip(",")

def arrayToStringInt(data):
	outputStr = ''
	for d in data:
		if isinstance(d, float):
			d = int(d)
		outputStr = outputStr + "{:d}, ".format(d)
	outputStr = outputStr.strip()
	return outputStr.rstrip(",")

def listToString(lst):
	outputStr = ""
	for val in lst:
		outputStr = outputStr + "{}, ".format(val)
	outputStr = outputStr.strip()
	return outputStr.rstrip(",")

def list2rangesFormatter(start, end, step):
	return "{}-{}:{}".format(start, end, step)

def list2ranges(data):
	lst = data
	if type(data) is set:
		lst = list(data)
	lst = sorted(lst)
	n = len(lst)
	result = []
	resultVals = []
	scan = 0
	while n - scan > 2:
		step = lst[scan + 1] - lst[scan]
		if lst[scan + 2] - lst[scan + 1] != step:
			result.append(str(lst[scan]))
			resultVals.append([lst[scan], lst[scan], 0])
			scan += 1
			continue

		for j in xrange(scan+2, n-1):
			if lst[j+1] - lst[j] != step:
				result.append(list2rangesFormatter(lst[scan], lst[j], step))
				resultVals.append([lst[scan], lst[j], step])
				scan = j+1
				break
		else:
			result.append(list2rangesFormatter(lst[scan], lst[-1], step))
			resultVals.append([lst[scan], lst[-1], step])
			return ",".join(result)

	if n - scan == 1:
		result.append(str(lst[scan]))
		resultVals.append([lst[scan], lst[scan], 0])
	elif n - scan == 2:
		result.append(",".join(itertools.imap(str, lst[scan:])))
		resultVals.append([lst[scan], lst[scan], 0])
		resultVals.append([lst[scan+1], lst[scan+1], 0])

	return ",".join(result)
