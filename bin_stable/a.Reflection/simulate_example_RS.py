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
This is an example simulation for beginner users to try out the code
"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Main.simulation.reflection_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze

import SEAS_Utils as utils
import SEAS_Utils.common_utils.configurable as config
import SEAS_Utils.common_utils.data_plotter as plotter

from SEAS_Utils.common_utils.timer import simple_timer


def simulate(s,o,a):


    s.Timer = simple_timer(4)
    
    #create a base flat spectra based on planetary parameters
    s.Surface_g, s.Base_TS_Value = s.load_astrophysical_properties()
    
    # normalized pressure directly correspond to atmosphere layers.
    s.normalized_pressure = s.load_atmosphere_pressure_layers()
    
    # load mixing ration files and determine what molecules are added to the simulation
    # acquire the mixing ratio at each pressure (need to be interpolated to normalized pressure)
    s.normalized_molecules, s.MR_Pressure, s.molecule_MR = s.load_mixing_ratio()
    s.normalized_abundance = s.interpolate_mixing_ratio()
    
    # load temperature pressure profile
    s.TP_Pressure, s.TP_Temperature = s.load_TP_profile()        
    s.normalized_temperature = s.interpolate_TP_profile()
    
    # calculate the scale height for each layer of the atmosphere
    s.normalized_scale_height = s.calculate_scale_height()
    
    # load molecular cross section for main constituent of the atmosphere
    # will check first to see if the molecule is in the database 
    s.cross_db = s.check_molecules()
    s.nu, s.normalized_cross_section = s.load_molecule_cross_section()   
    
    s.normalized_stellar_spectra = s.load_stellar_spectra()
    
    s.Reference_Reflect_Signal = s.load_atmosphere_geometry_model()
    

    
    s.user_input["Plotting"]["Figure"]["y_multiplier"] = 1
    s.user_input["Plotting"]["Figure"]["Title"] = "Simulated Example Reflection Spectra for Pure Water Atmosphere"
    s.user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    s.user_input["Plotting"]["Figure"]["y_label"] = r"Flux Density (W/m^2)"    
    
    sim_plot = plotter.Simulation_Plotter(s.user_input)

    plt_ref_1 = sim_plot.plot_xy(s.nu,s.normalized_stellar_spectra,"Stellar Spectra")
    plt_ref_2 = sim_plot.plot_xy(s.nu,s.normalized_stellar_spectra*s.Reference_Reflect_Signal,"Reflection Spectra")

    sim_plot.set_legend([plt_ref_1, plt_ref_2])
    if utils.to_bool(user_input["Save"]["Plot"]["save"]):
        sim_plot.save_plot()
    else:
        sim_plot.show_plot()
        
        
        

if __name__ == "__main__":
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Example"
    user_input["Simulation_Control"]["DB_Name"]             = "cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "isothermal_300K.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "H2O_only.txt"

    user_input["Save"]["Plot"]["save"] = False
    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "Temp_H2O_Cross_Section.npy"
    
    

    simulation = theory.RS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    simulate(simulation, observer, analyzer)






