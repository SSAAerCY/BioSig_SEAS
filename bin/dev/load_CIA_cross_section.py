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
from SEAS_Main.atmosphere_effects.collision_induced_absorption import HITRAN_CIA_data_processor, select_molecule_cia
#import matplotlib.pyplot as plt


if __name__ == "__main__":
    

    molecule_list = ["H2"]
    data_file = select_molecule_cia(molecule_list)
    
    data = {}
    for filename in data_file:
        processor = HITRAN_CIA_data_processor(HITRAN_CIA, filename)
        temp,nu,xsec = processor.load()
        data[filename] = {}
        
        data[filename]["Temperature"] = temp
        data[filename]["nu"] = nu
        data[filename]["xsec"] = xsec
            
    print filename
    print data[filename]["Temperature"]


