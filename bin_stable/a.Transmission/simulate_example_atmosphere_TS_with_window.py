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

import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze

import SEAS_Utils as utils
import SEAS_Utils.common_utils.configurable as config
import SEAS_Utils.common_utils.data_plotter as plotter
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data
from SEAS_Utils.common_utils.timer import simple_timer

import SEAS_Aux.atmosphere_processes.TP_profile_generator as TPgen
import SEAS_Aux.atmosphere_processes.mixing_ratio_generator as mix


def simulate_window(s, o, a, Filename=""):
    
    
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



    window_depth = 0.1
    window_width = 100.


    s.nu_window = a.spectra_window(nu,trans,"T",window_depth, window_width,s.min_signal)
    
    print s.nu_window
    
    
    
    
    
    print "calc time", s.Timer.elapse()

    s.user_input["Plotting"]["Figure"]["Title"] = "Absorption and Atmospheric Window for Simulated Atmosphere of %s"%"_".join(s.normalized_molecules)
    s.user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    s.user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    

    s.user_input["Save"]["Window"]["save"] = True
    s.user_input["Save"]["Window"]["path"] = "../../output/Variety_Atmosphere_Window"
    s.user_input["Save"]["Window"]["name"] = "%s_Window_D%s_W%s.txt"%(Filename,window_depth,window_width)

    s.user_input["Save"]["Plot"]["path"] = "../../output/Variety_Atmosphere_Spectra"
    s.user_input["Save"]["Plot"]["name"] = "%s_Window_D%s_W%s.png"%(Filename,window_depth,window_width)

 
    
    sim_plot = plotter.Simulation_Plotter(s.user_input)
    #sim_plot.save_window(s.nu_window)
    sim_plot.plot_xy(nu,trans)
    sim_plot.plot_window(s.nu_window,"k", 0.2)
    sim_plot.save_plot()
    """
    sim_plot.plot_xy(nu,trans)
    sim_plot.plot_window(s.nu_window,"k", 0.2)
    sim_plot.show_plot()
    """
    return s.Transit_Signal

def simulation_prep(atmosphere,name,type="File"):
    
    
    if type == "File":
        Filename = "Test_Earth"
        Molecules = atmosphere["Molecules"]
        MRFile = atmosphere["Mixing_Ratio_File"]
        TPFile = atmosphere["TP_File"]
        
        return Filename, Molecules, TPFile, MRFile


    Type                = atmosphere["Type"]
    Surface_Temperature = atmosphere["Surface_Temperature"]
    Molecules           = atmosphere["Molecules"]
    Mixing_Raio         = atmosphere["Ratio"]
    
    try:
        Filler = atmosphere["Filler"]
    except:
        Filler = "N2"
    Filename = "_".join([name,Type[:3],"_".join(Molecules),str(hash(str(atmosphere.values()))%10**8)])
    
    if Type == "isothermal":
        TP_input = config.Configuration("../../bin_stable/TP_Profile/TP_selection.cfg")["Test Isothermal Atmosphere"]
        TP_Name = "TP_%s.txt"%Filename
        if not os.path.isfile(os.path.join(TP_Profile_Data,TP_Name)):
            TP_simulator = TPgen.temperature_pressure_profile_generator(TP_input, name=TP_Name)
            TP_simulator.generate()
            TP_simulator.save()
    else:
        print "Non isothermal bulk simulation not implemented yet"
        sys.exit()
    
    # generate Mixing Ratio Files
    ratio_input = {}
    for i,m in enumerate(Molecules):
        ratio_input[m] = {}
        ratio_input[m]["Surface_Ratio"]  = Mixing_Raio[i]
        ratio_input[m]["End_Ratio"]      = None
        ratio_input[m]["Type"]           = "constant"
        ratio_input[m]["Transition"]     = None
        ratio_input[m]["Start_Pressure"] = None
        ratio_input[m]["End_Pressure"]   = None 
    
    MR_Name = "MR_%s.txt"%Filename
    if not os.path.isfile(os.path.join(Mixing_Ratio_Data,MR_Name)):
        simulator = mix.mixing_ratio_generator(ratio_input,
                                               filler=True,
                                               filler_molecule=Filler,
                                               pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                               name=MR_Name,
                                               overwrite=True)
        simulator.generate()
        simulator.save()

    return Filename, Molecules,TP_Name, MR_Name


def calculate_simple_window():

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
    analyzer = analyze.Spectra_Analyzer(user_input)
    
    Raw_TS = simulate_window(simulation, observer, analyzer)          


def calculate_window_multiple_atmosphere():

    Timer = simple_timer()
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    atmo_input = config.Configuration("../../bin_stable/a.Main/atmosphere_prototype.cfg")
    atmosphere_types = atmo_input["Simulation"]
    
    for i, atmosphere in enumerate(atmosphere_types):
    
        Filename, Molecules, TP_Name, MR_Name = simulation_prep(atmosphere_types[atmosphere],atmosphere,"Simulation")
        
        print Filename
        
        user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
        user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
        user_input["Simulation_Control"]["TP_Profile_Name"]     = TP_Name
        user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = MR_Name
    
        user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
        
        simulation = theory.TS_Simulator(user_input)
        observer = observe.OS_Simulator(user_input)
        analyzer = analyze.Spectra_Analyzer(user_input)
        
        Raw_TS = simulate_window(simulation, observer, analyzer, Filename)  
    
if __name__ == "__main__":
    
    calculate_window_multiple_atmosphere()
       


       