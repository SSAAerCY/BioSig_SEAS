"""



"""

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_ML.MCMC.retrieve_on_water_VMR import generate_water_VMR,generate_mixing_ratio_file


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
    
    
    # load rayleigh scattering
    s.Rayleigh_enable = utils.to_bool(s.user_input["Atmosphere_Effects"]["Rayleigh"]["enable"])
    if s.Rayleigh_enable:
        s.normalized_rayleigh = s.load_rayleigh_scattering()   
    
    # calculate theoretical transmission spectra
    s.Reference_Transit_Signal = s.load_atmosphere_geometry_model()
    
    # calculate observed transmission spectra
    nu,ref_trans = o.calculate_convolve(s.nu, s.Reference_Transit_Signal)



    
    return nu,ref_trans

def run_mcmc_most_simplistic():
    
    """
    in the simplistic case, assume we want to constrain surface water vapor abundance... for clear atmosphere?
    Therefore only 2 free parameter:
    
        P_0 and VMR_0
    
    """
    
    
    
    return 0
    
def run_modern_earth():


    VMR_0 = 4*10**-2
    VMR_1 = 10**-5.5
    VMR_2 = VMR_1
    VMR_3 = 10**-7
    P_0 = 10**5
    P_1 = P_0*np.e**-2
    P_2 = P_0*np.e**-10
    P_3 = P_0*np.e**-15

    PRange,VMR_Range = generate_water_VMR(VMR_0,VMR_1,VMR_2,VMR_3,P_0,P_1,P_2,P_3)
    
    pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 
                 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 
                 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05]

    f = interp1d(PRange, VMR_Range)
    VMR_interpolate = f(pressures)    
    
    VMR = ["H2O"]
    for i in VMR_interpolate:
        VMR.append(float("%.4g"%i))
    
    MR_name = "temp_RE.txt"
    
    generate_mixing_ratio_file(VMR,MR_name)
    
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename = "Test_Earth"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "earth.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = os.path.join("Retrieval",MR_name)

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
    
    simulation = theory.TS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    nu,ref_trans = simulate_NIST(simulation, observer, analyzer)   
  

    return nu,ref_trans


def run_uniform_model():

    VMR_0 = 2*10**-4
    VMR_1 = 2*10**-4
    VMR_2 = VMR_1
    VMR_3 = 10**-7
    P_0 = 10**5
    P_1 = P_0*np.e**-25
    P_2 = P_0*np.e**-25
    P_3 = P_0*np.e**-25


    PRange,VMR_Range = generate_water_VMR(VMR_0,VMR_1,VMR_2,VMR_3,P_0,P_1,P_2,P_3)
    
    pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 
                 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 
                 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05]

    f = interp1d(PRange, VMR_Range)
    VMR_interpolate = f(pressures)    
    
    VMR = ["H2O"]
    for i in VMR_interpolate:
        VMR.append(float("%.4g"%i))
    
    MR_name = "temp_RE.txt"
    
    generate_mixing_ratio_file(VMR,MR_name)
    
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename = "Test_Earth"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "earth.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = os.path.join("Retrieval",MR_name)

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
    
    simulation = theory.TS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    nu,ref_trans = simulate_NIST(simulation, observer, analyzer)   
  

    return nu,ref_trans



if __name__ == "__main__":
    
    
    VMR_0 = 4*10**-2
    VMR_1 = 10**-5.5
    VMR_2 = VMR_1
    VMR_3 = 10**-7
    P_0 = 10**5
    P_1 = P_0*np.e**-2
    P_2 = P_0*np.e**-10
    P_3 = P_0*np.e**-15


    VMR_0 = 4*10**-3
    VMR_1 = 4*10**-3
    VMR_2 = VMR_1
    VMR_3 = 10**-7
    P_0 = 10**5
    P_1 = P_0*np.e**-25
    P_2 = P_0*np.e**-25
    P_3 = P_0*np.e**-25

    VMR_0 = 2*10**-2
    VMR_1 = 10*10**-5.5
    VMR_2 = VMR_1
    VMR_3 = 10*10**-5.2
    P_0 = 10**5
    P_1 = P_0*np.e**-2
    P_2 = P_0*np.e**-5
    P_3 = P_0*np.e**-23


    PRange,VMR_Range = generate_water_VMR(VMR_0,VMR_1,VMR_2,VMR_3,P_0,P_1,P_2,P_3)
    
    pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 
                 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 
                 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05]

    f = interp1d(PRange, VMR_Range)
    VMR_interpolate = f(pressures)    
    
    VMR = ["H2O"]
    for i in VMR_interpolate:
        VMR.append(float("%.4g"%i))
    
    MR_name = "temp_RE.txt"
    
    generate_mixing_ratio_file(VMR,MR_name)
    
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename = "Test_Earth"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "earth.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = os.path.join("Retrieval",MR_name)

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
    
    simulation = theory.TS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    nu,ref_trans = simulate_NIST(simulation, observer, analyzer)             
    
    
    nu,earth_trans = run_modern_earth()
    nu,uni_trans = run_uniform_model()

    user_input["Plotting"]["Figure"]["Title"] = "Simulated Earth Like exoplanet atmosphere with varying Water Vapor VMR"
    user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    
    
    sim_plot = plotter.Simulation_Plotter(user_input)
    
    
    b = sim_plot.plot_xy(nu,earth_trans,"Modern Earth.")
    c = sim_plot.plot_xy(nu,uni_trans,"Uniform Dist.")
    a = sim_plot.plot_xy(nu,ref_trans,"Test.")
    sim_plot.set_legend([a,b,c])
    
    sim_plot.show_plot()
    
    
    