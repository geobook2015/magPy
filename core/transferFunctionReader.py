# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 08:18:04 2016

@author: npop
"""
from os.path import basename, splitext
import numpy as np
import math

class TransferFunctionReader(object):

    ###################
    ### CONSTRUCTOR
    ##################
    def __init__(self, path):
        self.filepath = path
        self.data = {}
        self.readFile() 
        # magnetic permeability in nT . m / A
        # recall, E data is in mV/m and H data is in nT
        # so Z = E/H is in mV / m . nT
        # so 
        # units of resistance = Ohm = V / A
        self.magPerm = 400*math.pi      
        
    ###################
    ### GET FUNCTIONS
    ##################
    def getFilePath(self):
        return self.filepath

    def getFileName(self):
        base = basename(self.filepath)
        return splitext(base)[0]

    def getFreqs(self):
        return self.freq

    def getPeriods(self):
        return self.period

    def getComponent(self, component):
        return self.data[component]

    def getVariance(self, component):
        return self.data["Var" + component]

    def getResAndPhase(self, component):
        # the data was in
        # E = mV
        # H = nT
        data = self.getComponent(component)
        res = 0.2 * self.period * np.power(np.absolute(data), 2)
        phase = np.angle(data)
        # check, can i unwrap into specific quadrant 
        phase = np.unwrap(phase)
        # convert to degrees
        phase = phase*180/math.pi
        if component == "ExHx" or component == "ExHy":
            # do modulo 360, see if this helps
            phase = np.mod(phase, 360)        
            phase = phase - 180
        return res, phase

    def getResAndPhaseErrors(self, component):
        # the data was in
        # E = mV
        # H = nT
        data = self.getComponent(component)
        res, phase = self.getResAndPhase(component)
        var = self.getVariance(component)
        resError = 1.96*np.sqrt(2 * self.period * res * var / 5.0 )      
        phaseError = 1.96*(180/math.pi) * (np.sqrt( var / 2.0 )/np.absolute(data))   
        # # calculate uncertainties on apparent resistivity and phase
        # drho = absData*var
        # dphase = np.arcsin(var/absData)
        # now calculate confidence intervals
        return resError, phaseError   

    # def getResAndPhaseErrors(self, component):
    #     # the data was in
    #     # E = mV
    #     # H = nT
    #     resVar, phaseVar = self.getResAndPhaseVariances(component) 
    #     # this assumes gaussian
    #     return resVar*1.96, phaseVar*1.96         
    
    ###################
    ### SET FUNCTIONS
    ##################
    def setPath(self, path):
        self.filepath = path
        self.readFile()
    
    ###################
    ### SET FUNCTIONS
    ##################   
    def readFile(self):
        key = ''
        with open(self.filepath, 'r') as inFile:
            for line in inFile:
                line = line.strip()
                if line == '':
                    continue

                if "Evaluation frequencies" in line:
                    key = 'evalFreq'
                    continue
                if "Var" in line:  
                    split = line.split("=")
                    key = "Var" + split[2].strip()
                    continue                  
                if "IJ" in line:
                    split = line.split("=")
                    key = split[2].strip()
                    continue             
                # otherwise, append the data
                if key == 'evalFreq':
                    self.freq = np.loadtxt(line.strip().split(","), dtype='float') 
                    self.period = np.reciprocal(self.freq)
                else: 
                    self.data[key] = np.loadtxt(line.strip().split(","), dtype='complex') 





	