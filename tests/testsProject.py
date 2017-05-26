#!/usr/bin/python
import os
import sys
sys.path.append(os.path.join("..", "core"))
sys.path.append(os.path.join("..", "utils"))
# import project
from project import Project
# import utils
from utilsIO import *

def setupTestProject():
	projectPath = os.path.join("..","..","testProject")
	projectName = "test"
	# the reftime for the test project
	refTime = "2016-01-18 00:00:00"

	# need to check and make project path
	checkAndMakeDir(projectPath)

	proj = Project()
	proj.setTimeDataPath(os.path.join(projectPath, "timeData"))
	proj.setSpecDataPath(os.path.join(projectPath, "specData"))
	proj.setStatDataPath(os.path.join(projectPath, "statData"))
	proj.setTransDataPath(os.path.join(projectPath, "transFuncData"))
	proj.setCalDataPath(os.path.join(projectPath, "calData"))
	proj.setRefTime(refTime)
	proj.initialiseProject()
	proj.saveProjectFile(os.path.join(projectPath, "{}.prj".format(projectName)))
	proj.printInfo()
	sites = proj.getAllSites()

	# lets print out information about the sites
	for s in sites:
		proj.printSiteInfo(s)
		timeFiles = proj.getSiteTimeFiles(s)
		for tFile in timeFiles:
			proj.printMeasInfo(s, tFile)

def setupEthiopiaProject():
	projectPath = os.path.join("..","..","ethiopiaProject")
	projectName = "ethiopia"
	# the reftime for ethiopia project
	refTime = "2012-02-10 00:00:00"

	# need to check and make project path
	checkAndMakeDir(projectPath)

	proj = Project()
	proj.setTimeDataPath(os.path.join(projectPath, "timeData"))
	proj.setSpecDataPath(os.path.join(projectPath, "specData"))
	proj.setStatDataPath(os.path.join(projectPath, "statData"))
	proj.setTransDataPath(os.path.join(projectPath, "transFuncData"))
	proj.setCalDataPath(os.path.join(projectPath, "calData"))
	proj.setRefTime(refTime)
	proj.initialiseProject()
	proj.saveProjectFile(os.path.join(projectPath, "{}.prj".format(projectName)))
	proj.printInfo()
	sites = proj.getAllSites()

	# lets print out information about the sites
	for s in sites:
		proj.printSiteInfo(s)
		timeFiles = proj.getSiteTimeFiles(s)
		for tFile in timeFiles:
			proj.printMeasInfo(s, tFile)			


##########
### RUN
##########
# setupTestProject()
setupEthiopiaProject()