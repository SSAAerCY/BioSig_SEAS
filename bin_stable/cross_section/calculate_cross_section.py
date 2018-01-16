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

This is an demonstration code for calculating cross sections from line lists


"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.cross_section.cross_section_calculator as csc
#import SEAS_Utils.common_utils.data_plotter as plt
from SEAS_Utils.common_utils.DIRs import Temp_DIR, HITRAN_Water_Lines,HITRAN_Lines
from SEAS_Utils.common_utils.timer import simple_timer

import matplotlib.pyplot as plt


def simple_test_inputs():

    d_path      = HITRAN_Water_Lines
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

    d_path      = HITRAN_Water_Lines
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
    
    #plt.simple_plot(nu,coef)

def simple_time_test():
    
    d_path      = HITRAN_Water_Lines
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
    
    #plt.simple_plot(nu,coef)
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
    
    #plt.simple_plot(nu,coef)






if __name__ == "__main__":

    d_path      = os.path.join(HITRAN_Lines,"CH4")
    r_path      = Temp_DIR
    molecule    = "CH4"
    component   = [1,1,1]
    numin       = 2500
    numax       = 5000
    step        = 2
    P           = 10.
    T           = 500.
    
    # 400-2000, 2000-10000, 10000-30000
    
    calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step,P,T)
    
    Timer = simple_timer()

    data = calc.read_data()
    
    nu, coef = calc.personal_calculator(data)
    wav = 10000./nu
    
    #plt.yscale("log")
    plt.plot(nu,coef)
    plt.show()
    
    k = 1.38*10**-23
    n = P/(k*T)
    l = (2*3.14*15320*1000*200*1000)**0.5
    
    print l
    import numpy as np
    
    from scipy import stats
    trans = np.exp(-n*coef*l)
    bin_means, bin_edges, binnumber = stats.binned_statistic(wav, trans, bins=500)
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2

    np.save("temp.txt",[bin_centers,bin_means])
    
    data = np.load("temp.txt")
    bin_centers, bin_means = data


    #plt.plot(bin_centers,bin_means)
    #plt.plot(10000./nu,np.exp(-n*coef*l))
    #plt.show()
    #plt.simple_plot(nu,coef)



