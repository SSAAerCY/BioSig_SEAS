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

def load_single(user_input):

    nu, ref_trans = load_npy(user_input["Save"]["Spectra"]["path"],
                             user_input["Save"]["Spectra"]["name"])
    
    sim_plot = plotter.Simulation_Plotter(user_input)
    sim_plot.plot_xy(nu,ref_trans,"ref.")
    sim_plot.show_plot()
    

def load_bio(user_input):

    nu, ref_trans, bio_trans = load_npy(user_input["Save"]["Spectra"]["path"],
                                        user_input["Save"]["Spectra"]["name"])
    
    sim_plot = plotter.Simulation_Plotter(user_input)
    sim_plot.plot_xy(nu,bio_trans,"bio.")
    sim_plot.plot_xy(nu,ref_trans,"ref.")
    
    sim_plot.show_plot()

if __name__ == "__main__":
    
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename1 = "Test_Earth_with_biosig"

    user_input["Plotting"]["Figure"]["Title"] = "Transit Signal and Atmospheric Window for Simulated Earth Atmosphere"
    user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    

    user_input["Save"]["Spectra"]["path"] = "../../output/Simulated_Spectra"
    user_input["Save"]["Spectra"]["name"] = "%s_Spectra_with_Biosig_Data"%Filename1
    

    load_bio(user_input)
    
    