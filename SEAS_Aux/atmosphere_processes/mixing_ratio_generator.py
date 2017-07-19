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
This is a logistic module that generates the mixing ratio files
This module does not handle the physics and chemistry that goes into 
determining the actual mixing ratio of the atmosphere

how do you handle variable mixing ratios?


"""

import os
import sys
import numpy as np

from SEAS_Utils import to_float

class mixing_ratio_generator():
    
    
    def __init__(self,
                 ratio_input,
                 filler = "N2",
                 pressures = [100000,10000,1000,100,10,1,0.1,0.01,0.001,0.0001,0.00001],
                 path = "../../input/atmosphere_data/Mixing_Ratio",
                 name = "Temp.txt"
                 ):
        
        
        self.ratio_input = ratio_input
        self.filler = filler
        self.pressures = pressures
        self.path = path
        self.name = name
        
    
    def generate(self):
        
        # create pressures
        self.data = [["Pressure"]]
        for P in self.pressures:
            self.data.append([str(P)])
        
        Surface_Pressure   = self.pressures[0]
        
        # add each molecules
        for Molecule in self.ratio_input:
            
            self.data[0].append(Molecule)
            Surface_Ratio  = to_float(self.ratio_input[Molecule]["Surface_Ratio"])
            Type           = self.ratio_input[Molecule]["Type"]
            Transition     = self.ratio_input[Molecule]["Transition"]
            Start_Pressure = to_float(self.ratio_input[Molecule]["Start_Pressure"]) 
            End_Pressure   = to_float(self.ratio_input[Molecule]["End_Pressure"]) 
            End_Ratio      = to_float(self.ratio_input[Molecule]["End_Ratio"]) 
            
            for j,pres in enumerate(self.pressures):
                    
                if Type == "constant":
                    self.data[j+1].append(str(Surface_Ratio))
                elif Type in ["decrease","increase"]:
                    if pres >= Start_Pressure:
                        self.data[j+1].append(str(Surface_Ratio))
                    elif pres <= End_Pressure:
                        self.data[j+1].append(str(End_Ratio))                   
                    else:
                        current = Surface_Ratio+(End_Ratio-Surface_Ratio)*(1-np.log10(pres/End_Pressure)/np.log10(Start_Pressure/End_Pressure))
                        self.data[j+1].append(str(current))
        
        # assuming single filler for now
        self.data[0].append(self.filler)
        for k,ratio in enumerate(self.data[1:]):
            total_ratio = sum([float(x) for x in ratio[1:]])
            if total_ratio > 100:
                print "Total Mixing Ratio exceed maximum, check mixing ratio generation"
                print self.data[0]
                print ratio
                sys.exit()
            self.data[k+1].append(str(100-total_ratio)) 
            
        return self.data
    
    def save(self):
        
        save_path = os.path.join(self.path,self.name)
        # need a file/path checker
        
        with open(save_path,"w") as file:
            for i in self.data:
                print i
                file.write(" ".join(i))
                file.write("\n")
            
            file.save()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    