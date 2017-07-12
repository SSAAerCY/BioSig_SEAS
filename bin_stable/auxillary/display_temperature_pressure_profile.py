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
Display the current TP Profile used
 

"""

import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import TP_Profile_Data
from SEAS_Utils.common_utils.data_loader import two_column_file_loader
import SEAS_Utils.common_utils.data_plotter as plt

if __name__ == "__main__":
    

    user_input = config.Configuration("user_input_dev.cfg")
    
    if user_input["Simulation_Control"]["TP_Profile"] == "File":
        filename = os.path.join(TP_Profile_Data,user_input["Simulation_Control"]["TP_Profile_Name"])
        
        x,y = two_column_file_loader(filename)
        plt.simple_plot(y,np.log10(x))














    
    