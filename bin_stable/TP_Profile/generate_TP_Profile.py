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
This code test the generation of TP Profiles
"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.atmosphere_processes.TP_profile_generator as TPgen
import SEAS_Utils.common_utils.configurable as config

if __name__ == "__main__":
    
    TP_input = config.Configuration("TP_selection.cfg")["Test Isothermal Atmosphere"]


    simulator = TPgen.temperature_pressure_profile_generator(TP_input, name="isothermal_300K.txt")
    simulator.generate()
    
    simulator.save()
    
    
    
    
    
    
    
    
    
    
    
    