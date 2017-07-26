import os
import sys
sys.path.append(os.path.join("..", "core"))
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
import matplotlib.pyplot as plt
# import utils
from utilsProcess import *
from utilsIO import *
# graphing
import matplotlib.pyplot as plt

# def removeControl(inString):
#     import re
#     # return re.sub(r'[\x00-\x1f\x7f-\x9f]', '', inString)
#     return re.sub(r'[\x00-\x05\xff]', '', inString)
#
# def twosComplement(dataBytes):
#     # two's complement 24-bit integer, little endian, unsigned and signed
#     unsigned = struct.unpack('<I', dataBytes + '\x00')[0]
#     signed = unsigned if not (unsigned & 0x800000) else unsigned - 0x1000000
#     return signed
#
# def readTag(dataFile):
#     second = struct.unpack('b', dataFile.read(1))[0]
#     minute = struct.unpack('b', dataFile.read(1))[0]
#     hour = struct.unpack('b', dataFile.read(1))[0]
#     day = struct.unpack('b', dataFile.read(1))[0]
#     month = struct.unpack('b', dataFile.read(1))[0]
#     year = struct.unpack('b', dataFile.read(1))[0]
#     dayOfWeek = struct.unpack('b', dataFile.read(1))[0]
#     century = struct.unpack('b', dataFile.read(1))[0]
#     print "{:02d}:{:02d}:{:02d} {:02d}/{:02d}/{:02d}".format(hour, minute, second, day, month, year)
#     # serial number
#     serialNum = struct.unpack('h', dataFile.read(2))
#     # num scans
#     numScans = struct.unpack('h', dataFile.read(2))[0]
#     print "numScans = {}".format(numScans)
#     # channels per scan
#     numChans = struct.unpack('b', dataFile.read(1))[0]
#     print "numChans = {}".format(numChans)
#     # tag length
#     tagLength = struct.unpack('b', dataFile.read(1))
#     # status code
#     statusCode = struct.unpack('b', dataFile.read(1))
#     # bit-wise saturation flags
#     saturationFlag = struct.unpack('b', dataFile.read(1))
#     # reserved
#     reserved = struct.unpack('b', dataFile.read(1))
#     # sample length
#     sampleLength = struct.unpack('b', dataFile.read(1))
#     # sample rate
#     sampleRate = struct.unpack('h', dataFile.read(2))
#     # units of sample rate: 0 = Hz, 1 = minute, 2 = hour, 3 = day
#     sampleUnits = struct.unpack('b', dataFile.read(1))
#     # clock status
#     clockStatus = struct.unpack('b', dataFile.read(1))
#     # clock error in micro seconds
#     clockError = struct.unpack('i', dataFile.read(4))
#     # reserved
#     res1 = struct.unpack('b', dataFile.read(1))
#     res2 = struct.unpack('b', dataFile.read(1))
#     res3 = struct.unpack('b', dataFile.read(1))
#     res4 = struct.unpack('b', dataFile.read(1))
#     res5 = struct.unpack('b', dataFile.read(1))
#     res6 = struct.unpack('b', dataFile.read(1))
#
#     return numScans, numChans
#
# def readRecord(dataFile, numChans, numScans):
#     print "samples to read = {}".format(numChans*numScans)
#     data = np.zeros(shape=(numChans, numScans), dtype='int')
#     for scan in xrange(0, numScans):
#         for chan in xrange(0, numChans):
#             dataBytes = dataFile.read(3)
#             data[chan, scan] = twosComplement(dataBytes)
#     return data
#
# def readTable(tablePath):
#     tableFile = open(tablePath, "rb")
#     tableData = {}
#     for i in xrange(0, 138):
#         # header = tableFile.read(4)
#         # value = tableFile.read(21)
#         header = struct.unpack('12s', tableFile.read(12))
#         header = removeControl(header[0])
#         # value = struct.unpack('21b', tableFile.read(21))
#         ints1 = ["EGN", "EGNC", "HGN", "MTSR", "L2NS", "L3NS", "L4NS", "TXPR", "TBVO", "TBVI", "INIT", "RQST", "MODE", "AQST", "HSMP", "CCLS",
#             "TEMP", "TMAX", "CHEX", "CHEY", "CHHX", "CHHY", "CHHZ", "NREF", "CCLT", "PZLT", "NSAT", "OCTR", "CLST", "TALS", "TCMB", "TOTL"
#         ]
#         ints2 = ["SNUM", "MXSC", "SATR", "BADR", "BAT1","BAT2","BAT3", "EXR", "EYR", "ELEV", "SRL2", "SRL3", "SRL4", "SRL5", "LPFR", "LFRQ"]
#         ints4 = ["STDE", "STDH"]
#         ints8 = []
#         ints1_4 = []
#         ints1_8 = ["TDSP", "LFIX", "TSYN", "STIM", "ETIM", "HTIM", "ETMH", "NUTC", "FTIM", "LTIM"]
#         floats = []
#         doubles = ["EXAC", "EXDC", "EYAC", "EYDC", "HXAC", "HXDC", "HYAC", "HYDC", "HZAC", "HZDC", "EXNR", "EXPR", "EYNR", "EYPR", "GNDR",
#             "TSTV", "FSCV", "CCMN", "CCMX", "HATT", "HAMP", "LFIX", "EXLN", "EYLN", "TSTR", "INPR", "CFMN", "CFMX", "HNOM"
#         ]
#         if header in ints1:
#             value = struct.unpack('b', tableFile.read(1))
#             tableFile.read(12)
#             test = value[0]
#         elif header in ints2:
#             value = struct.unpack('h', tableFile.read(2))
#             tableFile.read(11)
#             test = value[0]
#         elif header in ints4:
#             value = struct.unpack('i', tableFile.read(4))
#             tableFile.read(9)
#             test = value[0]
#         elif header in ints8:
#             value = struct.unpack('q', tableFile.read(8))
#             tableFile.read(5)
#             test = value[0]
#         elif header in ints1_4:
#             value = struct.unpack('4b', tableFile.read(4))
#             tableFile.read(9)
#             test = value
#         elif header in ints1_8:
#             value = struct.unpack('8b', tableFile.read(8))
#             tableFile.read(5)
#             test = value
#         elif header in floats:
#             value = struct.unpack('f', tableFile.read(4))
#             tableFile.read(9)
#             test = value[0]
#         elif header in doubles:
#             value = struct.unpack('d', tableFile.read(8))
#             tableFile.read(5)
#             test = value[0]
#         else:
#             value = struct.unpack('13s', tableFile.read(13))
#             test = removeControl(value[0])
#             continue
#         # test
#         # print "{} = {} - {}".format(header, value, test)
#         print "{} = {}".format(header, test)
#         tableData[header] = value
#     tableFile.close()
#     return tableData
#
# def readCoil(coilPath):
#     coilFile = open(coilPath, "rb")
#     coilData = {}
#     print struct.unpack('20b', coilFile.read(20))
#     print struct.unpack('12s', coilFile.read(12))
#     print struct.unpack('12s', coilFile.read(12))
#     print struct.unpack('20b', coilFile.read(20))
#     print struct.unpack('12s', coilFile.read(12))
#     print struct.unpack('d', coilFile.read(8))
#     print struct.unpack('d', coilFile.read(8))
#     print struct.unpack('44s', coilFile.read(44))
#     print struct.unpack('f', coilFile.read(4))
#     print struct.unpack('f', coilFile.read(4))
#     # value = struct.unpack('72s', coilFile.read(72))
#     # print value
#     # for i in xrange(0, 50):
#     #     # header = tableFile.read(4)
#     #     # value = tableFile.read(21)
#     #     # header = struct.unpack('8s', coilFile.read(8))
#     #     # value = struct.unpack('22s', coilFile.read(22))
#     #     # print "{} = {} - {}".format(header, value, value[0])
#     #     for offset in xrange(0,500):
#     #         coilFile.seek(0)
#     #         coilFile.read(offset)
#     #         value = struct.unpack('f', coilFile.read(4))[0]
#     #         if (value > 0.0000001 and value < 100000000):
#     #             print "offset {} - float {}".format(offset, value)
#     #     #coilData[header] = value
#     coilFile.close()
#     return coilData
#
# # read table file
# print "READING THE TABLE FILE"
# tablePath = os.path.join("..", "..", "Data", "riftVolc", "202", "2300407B.TBL")
# tabelVals = readTable(tablePath)
# exit()
#
# # coil file
# # print "READING THE COIL FILE"
# # coilPath = os.path.join("..", "..", "Data", "riftVolc", "202", "COIL2616.CLC")
# # coilVals = readCoil(coilPath)
# # exit()
#
# # now let's populate the
#
# # read data file
# dataPath = os.path.join("..", "..", "Data", "riftVolc", "202", "2300407B.TS5")
# dataFile = open(dataPath, "rb")
# scan = 0
# while scan < 2000:
#     print "------------------\nRECORD = {}".format(scan)
#     numScans, numChans = readTag(dataFile)
#     if scan == 0:
#         dataRecord = readRecord(dataFile, numChans, numScans)
#     else:
#         dataRecord = np.hstack((dataRecord, readRecord(dataFile, numChans, numScans)))
#     scan += 1
#     print "SHAPE = {}\n------------------".format(dataRecord.shape)
# dataFile.close()
#
# # now what I want to do is read the record and then the next
# # if the records are consecutive, do not write out a file
# # if they are non-consecutive, then write out a datafile and then move on
# # this will be phoenix to internal convertor
#
# # now try and plot
# plt.figure()
# for chan in xrange(0, numChans):
#     plt.subplot(numChans, 1, chan + 1)
#     # plt.plot(dataRecord[chan, 5000:20000] - np.average(dataRecord[chan, 5000:20000]))
#     plt.plot(dataRecord[chan] - np.average(dataRecord[chan]))
# plt.show()

# test the data reader
dataPath = os.path.join("..", "..", "Data", "riftVolc", "202")
dataReader = DataReaderPhoenix(dataPath)
dataReader.printInfo()
dataReader.printDataFileList()
print dataReader.getSamplesRatesTS()
print dataReader.getNumberSamplesTS()
dataReader.printTableFile()
startTime = "2017-04-07 23:00:00"
endTime = "2017-04-08 01:00:00"
data = dataReader.getUnscaledData(startTime, endTime)
# data = dataReader.getUnscaledSamples(startSample=1000, endSample=10000)
plt.figure()
for idx, chan in enumerate(data.keys()):
    plt.subplot(dataReader.getNumChannels(), 1, idx+1)
    plt.title(chan)
    plt.plot(data[chan]-np.average(data[chan]))
plt.show()
