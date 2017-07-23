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

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.DIRs import HITRAN_CIA
from SEAS_Main.atmosphere_effects.collision_induced_absorption import select_molecule_cia, HITRAN_CIA_data_processor
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)
    
if __name__ == "__main__":
    
    molecules = ["H2","H2"]
    
    #cia = cia_selector(molecules)
    #filename = cia.select()
    
    P = 100000
    T = 300
    K = 1.38*10**-23
    n = P/(K*T)  # in SI
    nc = n/10**6 # in CGS
    
    with open(os.path.join(HITRAN_CIA,"H2-H2_2011.cia")) as f:
        data = f.read().split("\n")[:9982]
        header = data[0]
        xdata,ydata = [list(x) for x  in zip(*[x.split() for x in data[1:]])]

        ydata = np.array(ydata,dtype=float)
        xdata = np.array(xdata,dtype=float)
        # need a data cropper
        numin = 300
        numax = 10000
        
        head_id,tail_id = 0,0
        head,tail = True,True
        count = 0
        for i in xdata:
            if i>=numin and head:
                head_id = count
                head = False
            if i>=numax and tail:
                tail_id = count
                tail = False
            count +=1
        
        if tail_id < head_id:
            print "error"
            sys.exit()
        
        xdata = xdata[head_id:tail_id]
        ydata = ydata[head_id:tail_id]
        
        # need universal plotting functions

        plt.title("Transmission Spectra for H2 CIA at 200K")
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel("Transmission") 

        ax = plt.gca()
        ax.set_xscale('log')
        plt.tick_params(axis='x', which='minor')
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))     


        
        plt.plot(10000/xdata,np.e**(-ydata*nc**2*100*100000))
        plt.show()
        
    
    
    