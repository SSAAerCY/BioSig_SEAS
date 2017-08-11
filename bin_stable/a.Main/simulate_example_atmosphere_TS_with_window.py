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

Generate a List of atmospheres and calculate their windows.

Hash is first introduced here for temporary file saving

"""
import os
import sys


DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data

from SEAS_Utils.common_utils.timer import simple_timer
import SEAS_Utils.common_utils.configurable as config
import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe

from SEAS_Utils.common_utils.timer import simple_timer


import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)

def simulate_window(s, o):
    
    
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
    
    s.normalized_rayleigh = s.load_rayleigh_scattering()
    
    print "load time", s.Timer.elapse()

    s.Transit_Signal = s.load_atmosphere_geometry_model()
    nu,trans = o.calculate_convolve(s.nu, s.Transit_Signal)

    s.nu_window = o.spectra_window(nu,trans,"T",0.3, 100.,s.min_signal)
    print "calc time", s.Timer.elapse()

    plt.title("Absorption and Atmospheric Window for Simulated Atmosphere of %s"%"_".join(s.normalized_molecules))
    plt.xlabel(r'Wavelength ($\mu m$)')
    plt.ylabel("absorption")    

    ax = plt.gca()
    ax.set_xscale('log')
    plt.tick_params(axis='x', which='minor')
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))          

    for k in s.nu_window:
        up,down = 10000./k[1],10000./k[0]
        plt.axvspan(up,down,facecolor="k",alpha=0.2)

    plt.plot(10000./nu,trans)
   
    plt.show()
    
    
    return s.Transit_Signal


if __name__ == "__main__":
    
    
    Timer = simple_timer()
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename = "Test_Earth"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "earth.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "earth.txt"

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
    user_input["Save"]["Window"]["path"] = "../../output/Simple_Atmosphere_Window"
    user_input["Save"]["Window"]["name"] = "%s_Window_A1000_S100.txt"%Filename
    
    user_input["Save"]["Plot"] = {}
    user_input["Save"]["Plot"]["path"] = "../../output/Plot_Result"
    user_input["Save"]["Plot"]["name"] = "%s_Plot.png"%Filename
    
    simulation = theory.TS_Simulator(user_input)
    observer = observe.OS_Simulator(user_input)
    
    Raw_TS = simulate_window(simulation, observer)         
    