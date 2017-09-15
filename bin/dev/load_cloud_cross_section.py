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
Test Cloud simulation 

"""


import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Main.atmosphere_effects.cloud import Physical_Cloud_Simulator



def test_cloud(lambd,radius,index):
               
    CS = Physical_Cloud_Simulator(lambd,radius)
    mat_abs,mat_sca,mat_qext,x = CS.spect(index)
    mat_sigma=CS.GetSigma(mat_sca)
    
    for i in mat_sigma.T:
        plt.plot(lambd,np.array(i)/100,'.-')

    CS.plot()



if __name__ == "__main__":
    
    
    """
    mean partical size, normal distribution
    refractive indices change with wavelength
    
    index = [real,imag]
    
    
    we need a particle number density
    pressure and temperature dependence?
    
    condensation rate? point? temperature? molecule?
    
    """
    
    index  = 1.33
    lambd  = np.arange(400,30000,10)
    lambd = 10000./lambd
    radius = [1]#np.arange([1,2,3,4]) # this should be a log normal distribution


    test_cloud(lambd,radius,index)
