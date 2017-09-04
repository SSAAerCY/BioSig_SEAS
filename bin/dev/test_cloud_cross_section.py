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
Code to test Gabi's cloud code


mie scattering due to different cloud partical radius?
is it molecule dependent?



"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Main.atmosphere_effects.cloud import Cloud_Simulator


def test_mie(lambd,radius,index):
    #Correct acording to Maetzler 2002
    
    CS = Cloud_Simulator(lambd,radius)
    x =(2.0*np.pi*index*radius)/lambd

    #CS.mie_abcd(1.58+1.0j, 2.1)
    #ext,sca,abso = CS.Mie(5+0.4j,1)
    ext,sca,abso = CS.Mie(x,index)
    print ext, sca, abso

def test_spect(lambd,radius,index):

    CS = Cloud_Simulator(lambd,radius)
    mat_abs,mat_sca,mat_qext,x = CS.spect(index)

    for i in mat_qext.T:
        plt.plot(lambd,i,'.-')

    CS.plot()        
            
def test_cloud(lambd,radius,index):
    
    CS = Cloud_Simulator(lambd,radius)
    mat_abs,mat_sca,mat_qext,x = CS.spect(index)
    mat_sigma=CS.GetSigma(mat_sca)
    
    for i in mat_sigma.T:
        plt.plot(lambd,i,'.-')

    CS.plot()


def test_cloud_cm(lambd,radius,index):
    
    CS = Cloud_Simulator(lambd,radius)
    mat_abs,mat_sca,mat_qext,x = CS.spect(index)
    mat_sigma_cm = np.true_divide(CS.GetSigma(mat_sca),10**(8))
    
    colorrange=np.arange(0.01,1,0.01)
    for count, i in enumerate(mat_sigma_cm.T):
        plt.plot(lambd,i,'.-', label='Radius=%s'%radius[count], color=str(colorrange[count]))
        
    CS.plot()



if __name__ == "__main__":
    
    index  = 1.33
    
    lambd  = 0.56
    radius = 0.2
    test_mie(lambd,radius,index)
    
    lambd  = [1.1,1.2,1.3,1.4]
    radius = np.arange(0.1,2,0.1)
    test_spect(lambd,radius,index)
    
    lambd  = np.arange(0.1,2.5,0.1)
    radius = np.arange(0.1,5,0.1)
    test_cloud(lambd,radius,index) 
    test_cloud_cm(lambd,radius,index)
    
    