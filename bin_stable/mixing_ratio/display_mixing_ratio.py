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
Display the current mixing ratio used
 

"""

import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data
from SEAS_Utils.common_utils.data_loader import multi_column_file_loader
#import SEAS_Utils.common_utils.data_plotter as plt

import matplotlib.pyplot as plt

if __name__ == "__main__":
    

    user_input = config.Configuration("../a.Main/user_input_dev.cfg")
    
    if user_input["Simulation_Control"]["Mixing_Ratio"] == "File":
        filename = os.path.join(Mixing_Ratio_Data,user_input["Simulation_Control"]["Mixing_Ratio_Name"])
        print filename
        
        data = multi_column_file_loader(filename,type="mixed")
        
        pressure = data[0]
        mixing_ratio = data[1:]
        print pressure
        print mixing_ratio
        
        P_0 = float(pressure[1])
        for i,MR in enumerate(mixing_ratio):
            
            plt.plot(MR[1:],-np.log(np.array(pressure[1:], dtype="float")/P_0))
            break

        plt.xscale("log")
        #plt.yscale("log")
        #plt.gca().invert_yaxis()
        plt.show()











    
    