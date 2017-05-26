# utility functions for readers
import os
import sys
import itertools
import glob
# import the different readers
from dataReaderATS import DataReaderATS
from dataReaderSpam import DataReaderSPAM
from dataReaderInternal import DataReaderInternal

# return the reader for a particular time data directory
def getDataReader(datapath):
	# check the file endings in the datapath
	# internal format
	headerF = glob.glob(os.path.join(datapath,"*.hdr"))
	headerF = headerF + glob.glob(os.path.join(datapath,"*.HDR"))
	if len(headerF) > 0: # then internal format
		return DataReaderInternal(datapath)
	headerF = glob.glob(os.path.join(datapath,"*.xml"))
	headerF = headerF + glob.glob(os.path.join(datapath,"*.XML"))
	if len(headerF) > 0: # then ats format
		return DataReaderATS(datapath)
	headerF = glob.glob(os.path.join(datapath,"*.xtrx"))
	headerF = headerF + glob.glob(os.path.join(datapath,"*.XTRX"))
	if len(headerF) > 0: # then xtrx format
		return DataReaderSPAM(datapath)
	headerF = glob.glob(os.path.join(datapath,"*.xtr"))
	headerF = headerF + glob.glob(os.path.join(datapath,"*.XTR"))
	if len(headerF) > 0: # then xtr format
		return DataReaderSPAM(datapath)
	# if nothing found, return false
	return False