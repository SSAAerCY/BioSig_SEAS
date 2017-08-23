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

Example of saving simulated spectra to npy files

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze

import SEAS_Utils as utils
import SEAS_Utils.common_utils.data_plotter as plotter
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.data_loader import load_npy



def simulate_noise(user_input):

    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    sim_plot = plotter.Simulation_Plotter(user_input)

    nu, ref_trans, bio_trans = load_npy(user_input["Save"]["Spectra"]["path"],
                                        user_input["Save"]["Spectra"]["name"])
    
    wav = (10000./nu)
    
    rbin_centers, rbin_mean_I = observer.calculate_bin(wav, ref_trans)
    rbin_mean_error_I, rbin_mean_error_bar = observer.add_noise(rbin_mean_I)

    bbin_centers, bbin_mean_I = observer.calculate_bin(wav, bio_trans)
    bbin_mean_error_I, bbin_mean_error_bar = observer.add_noise(bbin_mean_I)    

    signal_diff = np.array(bbin_mean_I - rbin_mean_I)
    signal_diff2 = np.array(bbin_mean_error_I - rbin_mean_I)
    
    x, SNR = analyzer.spectra_SNR(bbin_centers, signal_diff2)


    maxSNR = max(SNR)

    print "Feature Max SNR: %s. "%maxSNR,
    if maxSNR > 9:
        print "Over 9 sigma"
    elif maxSNR > 6:
        print "Over 6 sigma"        
    elif maxSNR > 3:
        print "Over 3 sigma"
    else:
        print "Not Detected"


    #sim_plot.plot_xy(nu, ref_trans)
    #sim_plot.plot_xy(wav,trans,"ref.",Dtype="um")
    
    #sim_plot.plot_bin(bin_centers, bin_mean_error_I, bin_mean_error_bar)
    #sim_plot.plot_bin(bin_centers, bin_mean_I, bin_mean_error_bar*10**-6)
    #sim_plot.plot_bin(bbin_centers, signal_diff, bbin_mean_error_bar*10**-6)
    #sim_plot.plot_bin(bbin_centers, signal_diff2, bbin_mean_error_bar*10**-6)   
    #sim_plot.plot_SNR(10000./x, SNR[::-1])        
    #sim_plot.show_plot()
    
    
    
    
    
    
    
    

if __name__ == "__main__":
    
    
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename1 = "Test_Earth_with_biosig"

    user_input["Plotting"]["Figure"]["Title"] = "Transit Signal and Atmospheric Window for Simulated Earth Atmosphere"
    user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    

    user_input["Save"]["Spectra"]["path"] = "../../output/Simulated_Spectra"
    user_input["Save"]["Spectra"]["name"] = "%s_Spectra_with_Biosig_Data"%Filename1  


    user_input["Observation_Effects"]["Bin"]["bin_number"] = 1000
    user_input["Observation_Effects"]["Noise"]["error_scale"] = 0.001

    simulate_noise(user_input)
    

    
    
    
    
    
    
        