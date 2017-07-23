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
This simulation aimed towards testing H2 atmospheres with CIA added
"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Main.simulation.observer as obs
import SEAS_Utils.common_utils.configurable as config
import SEAS_Utils.common_utils.data_plotter as plt


if __name__ == "__main__":
    
    user_input = config.Configuration("user_input_dev.cfg")
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation"
    user_input["Simulation_Control"]["DB_Name"]             = "cross_sec_simulation.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "isothermal_300K.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "H2&He.txt"


    user_input["Atmosphere_Effects"]["CIA"] = "true"

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "Temp_H2&He_Cross_Section.npy"

    simulation = theory.TS_Simulator(user_input)
    
    Raw_TS = simulation.simulate_CIA()
    
    
    
    
    
    