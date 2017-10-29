"""

Take 2 at trying to simulate a cloudy earth like atmosphere

"""


import os
import sys
import numpy as np
import random
import matplotlib.pyplot as plt

track = "../../.."
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, track))


import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze

import SEAS_Utils as utils
import SEAS_Utils.common_utils.data_plotter as plotter
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
from SEAS_Utils.common_utils.timer import simple_timer


def simulate_NIST(s,o,a):

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
    

    print "load time", s.Timer.elapse()
    
    # load rayleigh scattering
    s.Rayleigh_enable = utils.to_bool(s.user_input["Atmosphere_Effects"]["Rayleigh"]["enable"])
    if s.Rayleigh_enable:
        s.normalized_rayleigh = s.load_rayleigh_scattering()   
        
        
    user_input["Plotting"]["Figure"]["Title"] = "Simulated Earth-like Exoplanet Atmosphere with Physical ZnS Clouds (% condensed)"
    user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    
    user_input["Plotting"]["Figure"]["x_scale"] = "linear" 
    
    sim_plot = plotter.Simulation_Plotter(s.user_input)
    plt_legend = []
         
    # calculate theoretical transmission spectra
    s.Reference_Transit_Signal = s.load_atmosphere_geometry_model()
    nu,ref_trans = o.calculate_convolve(s.nu, s.Reference_Transit_Signal)
    plt_ref = sim_plot.plot_xy(nu,ref_trans,"ref.","r")
    plt_legend.append(plt_ref)     
    
    for condense in [1.,0.5,0.2,0.1,0.05,0.02,0.01]:
    
    
        mu = 0.05
        dev = 0.2
        s.Cloudy_Transit_Signal = s.load_atmosphere_geometry_model_with_cloud4(condense,mu,dev)
        nu,clo_trans = o.calculate_convolve(s.nu, s.Cloudy_Transit_Signal)
        plt_ref = sim_plot.plot_xy(nu,clo_trans,"%s"%(condense*100))
        plt_legend.append(plt_ref)        


    s.nu_window = a.spectra_window(nu,ref_trans,"T",0.3, 100.,s.min_signal)
    sim_plot.plot_window(s.nu_window,"k", 0.2)

    sim_plot.set_legend(plt_legend)

    sim_plot.show_plot()
    
    return s

    



if __name__ == "__main__":
    
    
    Timer = simple_timer()
    
    user_input = config.Configuration(os.path.join(track,"bin_stable/a.Main/user_input_dev.cfg"))
    
    Filename = "Test_Earth"
    Filename1 = "Test_Earth_with_biosig"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "earth.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "earth.txt"

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
    user_input["Save"]["Window"]["path"] = os.path.join(track,"output/Simple_Atmosphere_Window")
    user_input["Save"]["Window"]["name"] = "%s_Window_A1000_S100.txt"%Filename1
    
    user_input["Save"]["Plot"] = {}
    user_input["Save"]["Plot"]["save"] = False    
    user_input["Save"]["Plot"]["path"] = os.path.join(track,"output/Plot_Result")
    user_input["Save"]["Plot"]["name"] = "%s_Plot.png"%Filename1
    
    #user_input["Plotting"]["Figure"]["x_scale"] = "linear"
    
    simulation = theory.TS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    simulation = simulate_NIST(simulation, observer, analyzer)       







