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
This code test the generation of mixing_ratios
"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.atmosphere_processes.mixing_ratio_generator as mix
import SEAS_Utils.common_utils.configurable as config

if __name__ == "__main__":
    
    ratio_input = config.Configuration("mixing_ratio_selection.cfg")


    simulator = mix.mixing_ratio_generator(ratio_input)
    simulator.generate()
    simulator.save()
    
    
    
    
    
    
    
    
    
    
    
    