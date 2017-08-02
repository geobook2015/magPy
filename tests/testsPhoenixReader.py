import os
import sys
sys.path.append(os.path.join("..", "core"))
sys.path.append(os.path.join("..", "inbuilt"))
sys.path.append(os.path.join("..", "utils"))
import numpy as np
import math
from datetime import datetime, timedelta
import struct
# import readers
from dataReaderPhoenix import DataReaderPhoenix
from dataReaderInternal import DataReaderInternal
# import writers
from dataWriterInternal import DataWriterInternal
# import inbuilt
from projectDefault import *
from projectViewTime import *
# import utils
from utilsProcess import *
from utilsIO import *
# graphing
import matplotlib.pyplot as plt

# def readCoil(coilFile):

# coilPath = os.path.join("..", "..", "Data", "riftVolc", "202", "COIL1547.CLC")
# coilFile = open(coilPath, "rb")
# print struct.unpack("20b", coilFile.read(20))
# print struct.unpack("12s", coilFile.read(12))
# print struct.unpack("12s", coilFile.read(12))
# print struct.unpack("12s", coilFile.read(12))
# print struct.unpack("8s", coilFile.read(8))
# print struct.unpack("12s", coilFile.read(12))
# print struct.unpack("d", coilFile.read(8))
# print struct.unpack("d", coilFile.read(8))
# print struct.unpack("7s", coilFile.read(7))
# print struct.unpack("18s", coilFile.read(18))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("f", coilFile.read(4))
# print struct.unpack("d", coilFile.read(8))
# print struct.unpack("d", coilFile.read(8))
# print struct.unpack("d", coilFile.read(8))
# print struct.unpack("d", coilFile.read(8))
# print struct.unpack("d", coilFile.read(8))
# print struct.unpack("500s", coilFile.read(500))
# coilFile.close()

### test the data reader
dataPath = os.path.join("..", "..", "Data", "riftVolc", "202")
dataReader = DataReaderPhoenix(dataPath)
dataReader.printInfo()
dataReader.printDataFileList()
print dataReader.getSamplesRatesTS()
print dataReader.getNumberSamplesTS()
dataReader.printTableFile()
# startTime = "2017-04-07 23:00:00"
# endTime = "2017-04-08 01:00:00"
# data = dataReader.getUnscaledData(startTime, endTime)
# plt.figure()
# for idx, chan in enumerate(data.keys()):
#     plt.subplot(dataReader.getNumChannels(), 1, idx+1)
#     plt.title(chan)
#     plt.plot(data[chan]-np.average(data[chan]))
# plt.show()

### now try and reformat
# outpath = os.path.join("..", "..", "Data", "riftVolc", "202_reformat")
# dataReader.reformat(outpath)

### create a project
# projectPath = (os.path.join("..", "..", "Data", "riftVolcProject"))
# projectMakeDefault(projectPath, "2017-04-07 06:00:00")
# proj = projectLoad(projectPath, "mtProj.prj")

### let's look at some time
# projectViewTime(proj, "2017-04-08 02:00:00", "2017-04-08 04:30:00", freqs=[15], save=True, chans=["Ex", "Ey", "Hx", "Hy", "Hz"])
# projectViewTime(proj, "2017-04-07 09:16:00", "2017-04-07 09:16:16", freqs=[150], save=True, chans=["Ex", "Ey", "Hx", "Hy", "Hz"])
# projectViewTime(proj, "2017-04-07 09:33:00", "2017-04-07 09:33:01", freqs=[2400], save=True, chans=["Ex", "Ey", "Hx", "Hy", "Hz"])
