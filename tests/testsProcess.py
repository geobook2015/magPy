import sys
import os
sys.path.append(os.path.join('..','utils'))
import numpy as np
import matplotlib.pyplot as plt
from utilsProcess import *
import copy


def testRemoveZeros():
	arr = np.arange(0,2000)*2.0
	arr[50:80] = arr[50:80]*0
	arr[1600:1650] = arr[1600:1650]*0
	test = copy.deepcopy(arr)
	test = removeZerosSingle(test)
	# now plot
	x = np.arange(0,2000)
	fig = plt.figure()
	plt.plot(x, arr, label="Original")
	plt.plot(x, test, label="Remove zeros")
	plt.legend()

def testRemoveNans():
	size = 100
	arr = np.arange(0,100)*2.0
	arr[50:80] = np.nan
	test = copy.deepcopy(arr)
	test = removeNansSingle(test)
	# now plot
	x = np.arange(0,100)
	fig = plt.figure()
	plt.plot(x, arr, "x-", markersize=12, label="Original")
	plt.plot(x, test, "o--", label="Remove nans")
	plt.legend()

# run tests
testRemoveZeros() 
testRemoveNans()
plt.show()