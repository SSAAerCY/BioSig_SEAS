"""

We want to map out the n and k space for which both the water and cloud cross section are combined to be 1


looking for optical thickness of the atmosphere
"""

import os
import sys
import numpy as np
import random

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze

import SEAS_Utils as utils
import SEAS_Utils.common_utils.data_plotter as plotter
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
from SEAS_Utils.common_utils.timer import simple_timer

import matplotlib.pyplot as plt

def load_atmosphere_properties(s):
    
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
    
    #s.Transit_Signal_CIA = s.load_atmosphere_geometry_model(CIA=True)
    #nu,ref_trans = o.calculate_convolve(s.nu, s.Transit_Signal_CIA)
    #s.Transit_Signal = s.load_atmosphere_geometry_model()

    return s

def calculate_atmosphere(s,o):

    s.Cloudy_Transit_Signal = s.load_atmosphere_geometry_model(Cloud=True,CIA=True)
    nu,clo_trans = o.calculate_convolve(s.nu, s.Cloudy_Transit_Signal)
    
    
    wav = 10000./nu[::-1]
    trans = clo_trans[::-1]
    
    cutoff1_low = 0.95
    cutoff1_high = 1.05
    
    cutoff2_low = 1.35
    cutoff2_high = 1.45
    
    
    cutoff1 = []
    cutoff2 = []
    
    
    
    for i,info in enumerate(wav):
        if info >cutoff1_low and info <cutoff1_high:
            cutoff1.append(trans[i])
            
    for i,info in enumerate(wav):
        if info >cutoff2_low and info <cutoff2_high:
            cutoff2.append(trans[i])      
    
    final_1 = sum(cutoff1)/len(cutoff1)      
    final_104 = sum(cutoff2)/len(cutoff2)  
    
    
    return final_1,final_104
    """
    print final_1,final_104
            
    
    s.user_input["Plotting"]["Figure"]["Title"]   = r"Transit Signal and Atmospheric Window for Simulated Atmosphere of %s"%"_".join(s.normalized_molecules)
    s.user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    s.user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    
    
    sim_plot = plotter.Simulation_Plotter(s.user_input)
    
    #plt_ref_1 = sim_plot.plot_xy(s.nu,s.Transit_Signal,"Molecular")
    plt_ref_1 = sim_plot.plot_xy(nu,ref_trans,"Molecular+CIA")
    plt_ref_2 = sim_plot.plot_xy(nu,clo_trans,"ZnS")
        
    sim_plot.plot_dot(1,final_1)
    sim_plot.plot_dot(1.4,final_104)
    sim_plot.set_legend([plt_ref_1, plt_ref_2])
    sim_plot.show_plot()
    sys.exit()
    """
    


def bigfig():


    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_simulation.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "isothermal_300K.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "H2&More.txt"

    user_input["Planet"]["R_Planet"] = 2.7
    user_input["Planet"]["M_Planet"] = 6.0
    
    
    user_input["Atmosphere_Effects"]["CIA"]["enable"] = "true"

    #user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "Temp_H2&He_Cross_Section.npy"
    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "Temp_H2&More.npy"
    
    simulation = theory.TS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    #Raw_TS = simulation.simulate_CIA()
    
    
    simulation = load_atmosphere_properties(simulation)
    
    Rlist = []
    nlist = []
    radin = np.logspace(-2,1,10)
    imagine = [10**-4]
    cloudn = np.arange(1.3,2.5,0.2)
    cloudn = [1,2,3]
    radin = [0.01,0.1,1,10]
    
    figslope = np.zeros([len(radin),len(cloudn)])
    countrad = 0
    for meanrad in radin:
        countdark = 0
        for dark in cloudn:
            wavelength_list = [1.0,1.4]
            #Zhuchang put code here
            #rats = getratio(wavelength_list,[imagine],[dark],meanrad)
            
            print "hi"
            simulation.cloudn = dark
            simulation.cloudk = imagine[0]
            
            simulation.user_input["Atmosphere_Effects"]["Cloud"]["distribution"]["mean"] = meanrad
            
            rats = calculate_atmosphere(simulation,observer)
            
            slope = (rats[1]-rats[0])/(np.log(wavelength_list[1]-wavelength_list[0]))
            figslope[countrad,countdark]=slope
            Rlist.append(meanrad)
            nlist.append(dark)
            countdark = countdark+1
        countrad = countrad+1
    return figslope,Rlist,nlist

def contour_():
    
    
    figslope,Rlist,nlist = bigfig()
    
    
    plt.xscale('log')
    plt.scatter(Rlist,nlist,c = figslope, cmap = 'seismic')
    plt.colorbar()
    plt.xlabel('k')
    plt.ylabel('n')
    
    plt.show()
    #plt.clim([-.6,1])
    


if __name__ == "__main__":
    
    contour_()
    #bigfig()

    
    #Raw_TS = simulate_CIA(simulation, observer) 











