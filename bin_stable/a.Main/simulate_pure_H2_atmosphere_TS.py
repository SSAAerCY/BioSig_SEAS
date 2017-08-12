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

import SEAS_Utils as utils
import SEAS_Utils.common_utils.configurable as config
import SEAS_Utils.common_utils.data_plotter as plotter
from SEAS_Utils.common_utils.timer import simple_timer


def simulate_CIA(s):
    
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

    s.normalized_molecules, s.MR_Pressure, s.molecule_MR = s.load_mixing_ratio()
    s.normalized_abundance = s.interpolate_mixing_ratio()
    
    # load temperature pressure profile
    s.TP_Pressure, s.TP_Temperature = s.load_TP_profile()        
    s.normalized_temperature = s.interpolate_TP_profile()

    # calculate the scale height for each layer of the atmosphere
    s.normalized_scale_height = s.calculate_scale_height()        

    s.cross_db = s.check_molecules()
    s.nu, s.normalized_cross_section = s.load_molecule_cross_section()
    
    s.normalized_rayleigh = s.load_rayleigh_scattering()

    CIA_Enable = utils.to_bool(s.user_input["Atmosphere_Effects"]["CIA"]["enable"])
    if CIA_Enable == True:
        s.CIA_File, s.CIA_Data = s.load_CIA(["H2"])
        s.normalized_CIA = s.interpolate_CIA()
    
    s.Transit_Signal_CIA = s.load_atmosphere_geometry_model(CIA=True)
    s.Transit_Signal = s.load_atmosphere_geometry_model()
    
    
    user_input["Plotting"]["Figure"]["Title"]   = r"Transit Signal and Atmospheric Window for Simulated Atmosphere of %s"%"_".join(s.normalized_molecules)
    user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    
    
    sim_plot = plotter.Simulation_Plotter(s.user_input)
    
    plt_ref_1 = sim_plot.plot_xy(s.nu,s.Transit_Signal,"Molecular")
    plt_ref_2 = sim_plot.plot_xy(s.nu,s.Transit_Signal_CIA,"Molecular+CIA")
    
    sim_plot.set_legend([plt_ref_1, plt_ref_2])
    sim_plot.show_plot()



if __name__ == "__main__":
    
    user_input = config.Configuration("user_input_dev.cfg")
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_simulation.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "isothermal_300K.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "H2&He.txt"

    user_input["Planet"]["R_Planet"] = 10
    user_input["Planet"]["M_Planet"] = 300
    
    
    user_input["Atmosphere_Effects"]["CIA"]["enable"] = "true"

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "Temp_H2&He_Cross_Section.npy"

    simulation = theory.TS_Simulator(user_input)
    
    #Raw_TS = simulation.simulate_CIA()
    
    Raw_TS = simulate_CIA(simulation)
    
    
    
    