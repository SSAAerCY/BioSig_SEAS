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

This is an demonstration code for calculating cross sections


"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.cross_section.cross_section_calculator as csc
import SEAS_Utils.common_utils.data_plotter as plt
from SEAS_Utils.common_utils.DIRs import Temp_DIR,HITRAN_Water_lines
from SEAS_Utils.common_utils.timer import simple_timer



def simple_test_inputs():

    d_path      = HITRAN_Water_lines
    r_path      = Temp_DIR
    molecule    = "H2O"
    component   = [1,1,1]
    numin       = 3000
    numax       = 40000
    step        = 0.1
    P           = 100000.
    T           = 300.
    
    return 

def simplist_test():

    d_path      = HITRAN_Water_lines
    r_path      = Temp_DIR
    molecule    = "H2O"
    component   = [1,1,1]
    numin       = 3000
    numax       = 40000
    step        = 0.1
    P           = 100000.
    T           = 300.
    
    
    calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step,P,T)
    
    Timer = simple_timer()
    
    
    data = calc.read_data()
    
    nu, coef = calc.personal_calculator(data)
    
    plt.simple_plot(nu,coef)

def simple_time_test():
    
    d_path      = HITRAN_Water_lines
    r_path      = Temp_DIR
    molecule    = "H2O"
    component   = [1,1,1]
    numin       = 3000
    numax       = 40000
    step        = 0.1
    P           = 100000.
    T           = 300.
    
    
    calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step,P,T)
    
    Timer = simple_timer()
    
    
    data = calc.read_data()
    for i in range(10):
        nu, coef = calc.personal_calculator(data)
        print "hi", Timer.elapse()
    
    plt.simple_plot(nu,coef)
    print Timer.total()    

def test_hapi_calc():
    """
    
    Doesn't really work fast at the moment....
    
    """
    path = "../../input/absorption_data/HITRAN_Line_List/O3"
    molecule = ["O3",3,1]
    numin = 3900
    numax = 4050
    step = 1
    
    cross_section = csc.cross_section_calculator(path,molecule,numin,numax,step)
    nu, coef = cross_section.hapi_calculator()
    
    plt.simple_plot(nu,coef)
    



if __name__ == "__main__":

    d_path      = HITRAN_Water_lines
    r_path      = Temp_DIR
    molecule    = "H2O"
    component   = [1,1,1]
    numin       = 3000
    numax       = 4000
    step        = 0.1
    P           = 100000.
    T           = 300.
    
    calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step,P,T)
    
    Timer = simple_timer()

    data = calc.read_data()
    
    nu, coef = calc.personal_calculator(data)
    
    plt.simple_plot(nu,coef)



