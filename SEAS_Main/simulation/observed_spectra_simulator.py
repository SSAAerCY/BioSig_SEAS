#!/usr/bin/env python
#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

Functions related to simulating observed spectra based on calculated theoratical spectra

for 0.8, need a full scale conversion of all list into dicts
instead of normalized_xxx, let's have a dict with pressure_layers as keys and relevent data as data

Takes in a simulated theoretical spectra and add observational effects


"""
import os
import sys
import numpy as np
import time
from scipy import interpolate
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)



import SEAS_Aux.cross_section.hapi as hp

import SEAS_Main.observation_effects.observation_noise as noise


class OS_Simulator():
    
    
    def __init__(self, user_input):
        
        self.user_input = user_input
     
    def add_noise(self, nu, trans):
        pass
    
    def calculate_convolve(self, nu, trans):
        #if self.user_input["Observation"]["Convolve"] == "true":
        amount = float(self.user_input["Observation"]["Convolve_Amount"])
        nu,Transit_Signal,i1,i2,slit = hp.convolveSpectrum(nu,trans,SlitFunction=hp.SLIT_RECTANGULAR,Resolution=amount,AF_wing=20.0)
        
        return nu,Transit_Signal
    
    def spectra_window(self, nu, coef, type="A",threshold=200.,span=100.,min_signal=0):
        
        if type == "A":
        
            window = []
            
            start = True
            win = [0,0]
            
            self.stuff = []
            for i,n in enumerate(nu):
                
                if coef[i] < threshold and start == True:
                    win[0] = n
                    self.stuff.append(n)
                    start = False
                if coef[i] > threshold and start == False:
                    
                    self.stuff.append(n)
                    if n>=win[0]+span:
                        win[1] = n
                        start = True
                        window.append(np.array(win))
                    else:
                        win[0] = n
                        start = True
        
        elif type == "T":
            
            window = []
            start = True
            win = [0,0]
                        
            Min = min_signal
            Max = max(coef)#self.max_signal
            if threshold > 1:
                threshold = 1000/threshold
            
            
            threshold = Min+(Max-Min)*threshold
            self.threshold = threshold
            
            self.stuff = []
            for i,n in enumerate(nu):
                
                if coef[i] < threshold and start == True:
                    win[0] = n
                    self.stuff.append(n)
                    start = False
                if coef[i] > threshold and start == False:
                    
                    self.stuff.append(n)
                    if n>=win[0]+span:
                        win[1] = n
                        start = True
                        window.append(np.array(win))
                    else:
                        win[0] = n
                        start = True

        if self.user_input["Save"]["Window"]["save"] == "true":
            with open(os.path.join(self.user_input["Save"]["Window"]["path"],
                                   self.user_input["Save"]["Window"]["name"]),"w") as f:
                for i in window:
                    f.write("%s-%s\n"%(i[0],i[1]))
        
        return window

    def analyze_spectra_detection(self,nu,trans,bio_trans,method="max"):
        """
        How to implement area under curve?
        """
        
        noise_level = 10
        comp = 2
        detection = False
        Detected = []
        
        for i in self.nu_window:
            detected = False
            reference =  trans[list(nu).index(i[0]):list(nu).index(i[1])]
            signal = bio_trans[list(nu).index(i[0]):list(nu).index(i[1])]
            
            
            # above certain ppm
            difference = max(signal-reference)*10**6
            # above certain comparision
            comparison = max((signal-self.min_signal)/(reference-self.min_signal))
        
            if difference > 3*noise_level:
                detection = True
                detected = True
            if comparison > comp:
                detection = True
                detected = True
                
            
            Detected.append(detected)
                
        
        return detection, Detected
        
        





