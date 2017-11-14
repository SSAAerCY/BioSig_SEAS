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

need a better system for reading in what molecules we should generate.
"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.atmosphere_processes.mixing_ratio_generator as mix
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import HITRAN_Lines


def generate_h2():
    
    ratio_input = config.Configuration("selection/h2_only_selection.cfg")

    simulator = mix.mixing_ratio_generator(ratio_input,
                                           filler=True,
                                           filler_molecule="He",
                                           pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                           name="H2&He.txt",
                                           overwrite=True)
    simulator.generate()
    simulator.save()

def generate_water():

    ratio_input = config.Configuration("water_only_selection.cfg")
    simulator = mix.mixing_ratio_generator(ratio_input,filler=False,name="selection/water_only2.txt")
    simulator.generate()
    simulator.save()
    
    
def generate_iso_earth():

    ratio_input = config.Configuration("selection/well_mixed_earth.cfg")

    simulator = mix.mixing_ratio_generator(ratio_input,
                                           filler=True,
                                           filler_molecule="N2",
                                           pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                           name="iso_earth.txt",
                                           overwrite=True)
    simulator.generate()
    simulator.save()


def generate_h2_and_more():
    
    ratio_input = config.Configuration("selection/h2_and_more.cfg")

    simulator = mix.mixing_ratio_generator(ratio_input,
                                           filler=True,
                                           filler_molecule="He",
                                           pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                           name="H2&More.txt",
                                           overwrite=True)
    simulator.generate()
    simulator.save() 



def generate_all_hitran():
    
    ratio_input = config.Configuration("selection/all_HITRAN.cfg")

    simulator = mix.mixing_ratio_generator(ratio_input,
                                           filler=True,
                                           filler_molecule="He",
                                           pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                           name="All_HITRAN.txt",
                                           overwrite=True)
    simulator.generate()
    simulator.save() 
    
def generate_all_hitran_selection():
    
    formulas = os.listdir(HITRAN_Lines)
    
    for formula in formulas:
        try:
            if ".data" in " ".join(os.listdir(os.path.join(HITRAN_Lines,formula))):
                    print """[%s]
    Surface_Ratio  = 1
    End_Ratio      = None
    Type           = constant
    Transition     = None
    Start_Pressure = None
    End_Pressure   = None\n"""%formula  
        except:
            pass
    

if __name__ == "__main__":
    

    #generate_h2_and_more()
    #generate_all_hitran_selection()
    generate_all_hitran()
    
    
    
    
    
    
    
    
    
    