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
This is a clean empty example. 
 

"""


import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.configurable as config

import SEAS_Main.simulation.emission_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze


def test_ES_Sim(s,o,a):

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
    

    I = s.load_atmosphere_geometry_model()
    
    nu,ref_emi = o.calculate_convolve(s.nu, I)
    
    #ax = plt.gca()
    
    #ax.set_xscale("log") 
    
    """
    legends = []
    for T in [300,350,400,450]:
        BB = s.blackbody_lam(10000./s.nu, T)
        leg, = plt.plot(10000./s.nu, BB, label="%sK"%T)
        legends.append(leg)
    
    plt.xlabel("Wavelength (um)")
    plt.ylabel("Relative Intensity (erg/s/cm^2/cm/Steradian)")
    plt.title("Emission Spectra from Pure Water Atmosphere")
    
    leg, = plt.plot(10000./nu, ref_emi/np.median(ref_emi),"k", label="ES_ref")
    legends.append(leg)
    plt.legend(handles=legends)
    plt.show()
    """

    legends = []
    BS = s.blackbody_lam(10000./s.nu, 4000)
    nu,BS_ = o.calculate_convolve(s.nu, BS)
    
    for T in [300,350,400,450]:
        BB = s.blackbody_lam(10000./s.nu, T)
        nu, BB_ = o.calculate_convolve(s.nu, BB)
        leg, = plt.plot(10000./nu, BB_/BS_, label="%sK"%T)
        legends.append(leg)
    
    plt.xlabel("Wavelength (um)")
    plt.ylabel("I_P/I_S")
    plt.title("Emission Spectra from Earth Like Atmosphere")
    
    leg, = plt.plot(10000./nu, ref_emi/np.median(ref_emi)/BS_,"k", label="ES_ref")
    legends.append(leg)
    plt.legend(handles=legends)
    plt.show()    
    
    
    
    

if __name__ == "__main__":
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename = "Test_Earth"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "earth.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "earth.txt"

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename

    simulation = theory.ES_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)    
    test_ES_Sim(simulation,observer,analyzer)
    
    
    
    
    
    
    
    
    
    
    