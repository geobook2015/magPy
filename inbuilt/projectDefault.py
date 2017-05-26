#!/usr/bin/python

# imports
import sys
import os
sys.path.append(os.path.join('..', 'core'))
sys.path.append(os.path.join('..', 'utils'))
import numpy as np
from datetime import datetime
# import my classes
from project import Project
# import utils
from utilsIO import *

# make the default project setup
# user has to supply project path
# name is optional
def projectMakeDefault(projectPath, refTime, **kwargs):
	generalPrint("ProjectMakeDefault", "Creating a new project in path: {}".format(projectPath))
	name = "mtProj"
	
	if "name" in kwargs:
		name = kwargs["name"]
	
	# need to check and make project path
	checkAndMakeDir(projectPath)
	
	proj = Project()
	proj.setTimeDataPath(os.path.join(projectPath, "timeData"))
	proj.setSpecDataPath(os.path.join(projectPath, "specData"))
	proj.setStatDataPath(os.path.join(projectPath, "statData"))
	proj.setTransDataPath(os.path.join(projectPath, "transFuncData"))
	proj.setCalDataPath(os.path.join(projectPath, "calData"))
	proj.setImageDataPath(os.path.join(projectPath, "images"))	
	proj.setRefTime(refTime)
	proj.initialiseProject()
	proj.saveProjectFile(os.path.join(projectPath, "{}.prj".format(name)))
	proj.printInfo()
	sites = proj.getAllSites()

	# lets print out information about the sites
	for s in sites:
		proj.printSiteInfo(s)
	
	return proj

def projectLoad(projectPath, projectFile):
	generalPrint("ProjectLoad", "Loading project file: {}".format(os.path.join(projectPath, projectFile)))	
	proj = Project()
	proj.loadProjectFile(os.path.join(projectPath, projectFile))
	proj.printInfo()
	sites = proj.getAllSites()

	# lets print out information about the sites
	for s in sites:
		proj.printSiteInfo(s)

	return proj
	
