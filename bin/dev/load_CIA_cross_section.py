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

Read HITRAN CIA data and plot it

This reader currently works for a subset of the simulation

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.DIRs import HITRAN_CIA
from SEAS_Aux.atmosphere_processes.CIA_selector import cia_selector
from SEAS_Utils.common_utils.data_loader import two_column_file_loader,HITRAN_CIA_loader
#import matplotlib.pyplot as plt

class HITRAN_CIA_data():
    
    
    def __init__(self,filename):
        
        self.filename = filename
    
    def load(self):
        
        pass




if __name__ == "__main__":
    
    with open(os.path.join(HITRAN_CIA,"H2-H2_2011.cia")) as f:

        data = f.read().split("\n")
        header = data[0]
        print header
        
        formula      = header[:20].strip()
        numin        = float(header[20:30].strip())
        numax        = float(header[30:40].strip())
        data_point   = int(header[40:47].strip())
        temperature  = float(header[47:54].strip())
        maximum_xsec = float(header[54:64].strip())
        resolution   = float(header[64:70].strip())
        comments     = header[76:97].strip()
        reference    = header[97:].strip()

        
        if data[-1] == "":
            data = data[:-1]
        
        
        data_grid = np.reshape(np.array(data),(-1,data_point+1))
        
        for i in data_grid:
            print i[0]





