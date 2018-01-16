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
import numpy as np
from scipy import stats
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
from SEAS_Utils.common_utils.constants import *
from SEAS_Main.observation_effects.noise import Photon_Noise

import SEAS_Aux.atmosphere_processes.TP_profile_generator as TPgen
import SEAS_Aux.atmosphere_processes.mixing_ratio_generator as mix

def simulate_NIST(s,o,a,atmosphere,stuff, sensitivity,separation):

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
    
    # load biosignature molecules
    bio_enable = utils.to_bool(s.user_input["Atmosphere_Effects"]["Bio_Molecule"]["enable"])
    if bio_enable == True:
        data_type     = s.user_input["Atmosphere_Effects"]["Bio_Molecule"]["data_type"]
        bio_molecule  = s.user_input["Atmosphere_Effects"]["Bio_Molecule"]["molecule"]
        bio_abundance = utils.to_float(s.user_input["Atmosphere_Effects"]["Bio_Molecule"]["abundance"])
        s.is_smile      = utils.to_bool(s.user_input["Atmosphere_Effects"]["Bio_Molecule"]["is_smile"])
        s.nu, s.bio_cross_section = s.load_bio_molecule_cross_section(bio_molecule, data_type)
        

        s.bio_normalized_cross_section = np.concatenate([s.normalized_cross_section, s.bio_cross_section], axis=0)

        # modify the molecular abundance after adding biosignatures
        # scale heights untouched still since effect is small
        s.bio_normalized_abundance = []
        for i,abundance in enumerate(s.normalized_abundance):
            s.bio_normalized_abundance.append([])
            for j in abundance:
                s.bio_normalized_abundance[i].append(j*(1-bio_abundance))
            s.bio_normalized_abundance[i].append(bio_abundance)
        s.bio_normalized_molecules = np.concatenate([s.normalized_molecules,[bio_molecule]], axis=0)

    print "load time", s.Timer.elapse()
    
    # load rayleigh scattering
    s.Rayleigh_enable = utils.to_bool(s.user_input["Atmosphere_Effects"]["Rayleigh"]["enable"])
    if s.Rayleigh_enable:
        s.normalized_rayleigh = s.load_rayleigh_scattering()   
    
    
    s.Overlay_enable = utils.to_bool(s.user_input["Atmosphere_Effects"]["Overlay"]["enable"])
    if s.Overlay_enable:
        o_nu, o_xsec = s.load_overlay_effects()
        s.normalized_overlay = s.interpolate_overlay_effects(o_nu,o_xsec)
         
    # calculate theoretical transmission spectra
    s.Reference_Transit_Signal = s.load_atmosphere_geometry_model()
    s.Bio_Transit_Signal = s.load_atmosphere_geometry_model(bio=bio_enable)
    
    # calculate observed transmission spectra
    nu,ref_trans = o.calculate_convolve(s.nu, s.Reference_Transit_Signal)
    nu,bio_trans = o.calculate_convolve(s.nu, s.Bio_Transit_Signal)
    
    # analyze the spectra
    #s.nu_window = a.spectra_window(nu,ref_trans,"T",0.5, 100.,s.min_signal)
    #s.nu_window = a.new_spectra_window(nu,ref_trans,0.8,100)


    noise = Photon_Noise(s.user_input)
    diff = np.array(s.Bio_Transit_Signal)#-np.array(s.Reference_Transit_Signal)
    bin_edges, bin_width, bin_centers = noise.determine_bin()    
    
    bin_means, bin_edges, binnumber = stats.binned_statistic(10000./s.nu[::-1], diff[::-1], bins=bin_edges)

    new_bin_means = []
    for i,info in enumerate(bin_means):
        if float(info) < 0:
            new_bin_means.append(0)
        else:
            new_bin_means.append(info)
    bin_means = np.array(new_bin_means)    

    # what about SNR of atmosphere vs atmosphere ground?

    signal, photon_noise, SNR = noise.calculate_noise(bin_means)


    
    
    s.user_input["Plotting"]["Figure"]["Title"] = "Transit Signal and Atmospheric Window for Simulated %s with traces of %s at %s ppm\n%s"%(atmosphere,bio_molecule,bio_abundance*10**6,stuff)
    s.user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    s.user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    
    s.user_input["Plotting"]["Figure"]["x_scale"] = "linear"
    
    sim_plot = plotter.Simulation_Plotter(s.user_input)
    
    #plt_ref_1 = sim_plot.plot_xy(nu,bio_trans,"with_bio")
    plt_ref_2 = sim_plot.plot_xy(nu,ref_trans,"ref.")
    plt_ref_3 = sim_plot.plot_xy(10000./np.array(bin_centers),bin_means,"bio_binned.", Marker = "o")
    
    #sim_plot.plot_window(s.nu_window,"k", 0.2)
    sim_plot.set_legend([plt_ref_2, plt_ref_3])
    
    if utils.to_bool(s.user_input["Save"]["Plot"]["save"]):
        sim_plot.save_plot()
    else:
        sim_plot.show_plot()
    
    return s


if __name__ == "__main__":
    
    
    Timer = simple_timer()
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    atmo_input = config.Configuration("../../bin_stable/a.Main/atmosphere_prototype.cfg")
    atmosphere_types = atmo_input["Proposed"]
    
    for i, atmosphere in enumerate(atmosphere_types):
        
        name = atmosphere
        print "Simulating '%s' Atmosphere"%atmosphere
        Type                = atmosphere_types[atmosphere]["Type"]
        Molecules           = atmosphere_types[atmosphere]["Molecules"]
        Mixing_Raio         = atmosphere_types[atmosphere]["Ratio"]
        Filler              = atmosphere_types[atmosphere]["Filler"]
        TP_Name             = atmosphere_types[atmosphere]["TP_File"]
        Window_Threshold    = atmosphere_types[atmosphere]["Window_Threshold"]
        Window_Span         = atmosphere_types[atmosphere]["Window_Span"]
        Simple_Filename = "_".join([Type[:3],"_".join(Molecules),name])
        Complex_Filename = "_".join([Type[:3],"_".join(Molecules),str(hash(str(atmosphere_types[atmosphere].values()))%10**8)])  
        MR_Name = "MR_%s.txt"%Complex_Filename
        
        stuff = ""
        for mol,ct in zip(Molecules, Mixing_Raio):
            if stuff == "":
                stuff = "%s:%s"%(mol,ct)
            else:
                stuff = "%s, %s:%s"%(stuff,mol,ct)

        
        ratio_input = {}
        for i,m in enumerate(Molecules):
            ratio_input[m] = {}
            ratio_input[m]["Surface_Ratio"]  = Mixing_Raio[i]
            ratio_input[m]["End_Ratio"]      = None
            ratio_input[m]["Type"]           = "constant"
            ratio_input[m]["Transition"]     = None
            ratio_input[m]["Start_Pressure"] = None
            ratio_input[m]["End_Pressure"]   = None         
        
        if not os.path.isfile(os.path.join(Mixing_Ratio_Data,MR_Name)):
            simulator = mix.mixing_ratio_generator(ratio_input,
                                                   filler=True,
                                                   filler_molecule=Filler,
                                                   pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                                   name=MR_Name,
                                                   overwrite=True)
            simulator.generate()
            simulator.save()        
        
        Filename = Simple_Filename
        Filename1 = Complex_Filename
        
        user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
        user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
        user_input["Simulation_Control"]["TP_Profile_Name"]     = TP_Name
        user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = MR_Name
    

        """
        info = NIST_Smile_List()
        molecule_smiles = info[0]
        Bio_Molecule = random.choice(molecule_smiles)
        """
        Bio_Molecule = "O=CC"
        #Bio_Molecule = molecule_smiles[0]
        user_input["Atmosphere_Effects"]["Bio_Molecule"]["enable"] = True
        user_input["Atmosphere_Effects"]["Bio_Molecule"]["data_type"] = "NIST"
        user_input["Atmosphere_Effects"]["Bio_Molecule"]["molecule"] = Bio_Molecule
        user_input["Atmosphere_Effects"]["Bio_Molecule"]["abundance"] = 100*10**-6
        user_input["Atmosphere_Effects"]["Bio_Molecule"]["is_smile"] = True
        
        user_input["Atmosphere_Effects"]["Overlay"]["enable"] = True

        user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
        user_input["Save"]["Window"]["path"] = "../../output/Simple_Atmosphere_Window"
        user_input["Save"]["Window"]["name"] = "%s_Window_A1000_S100.txt"%Filename1
    
        user_input["Save"]["Plot"] = {}
        user_input["Save"]["Plot"]["save"] = True    
        user_input["Save"]["Plot"]["path"] = "../../result/Carbonyl_Plots"
        user_input["Save"]["Plot"]["name"] = "%s_Plot.png"%atmosphere


        
        user_input["Star"]["R_Star"] = 0.6
        user_input["Star"]["T"] = 4000.
        user_input["Planet"]["R_Planet"] = 1
        #user_input["Planet"]["R_Atmosphere"] = 40*1000
        user_input["Telescope"]["Aperture"] = 6.5
        user_input["Telescope"]["Distance"] = 10*Psec
        user_input["Telescope"]["Duration"] = 100*3600
        user_input["Telescope"]["Quantum_Efficiency"] = 0.3
        user_input["Telescope"]["min_wavelength"] = 1
        user_input["Telescope"]["max_wavelength"] = 25
        user_input["Observation_Effects"]["Noise"]["Multiplier"] = 1.2
        user_input["Observation_Effects"]["bin_exponent"] = 3/2.
        user_input["Observation_Effects"]["bin_width"] = 0.01
    
        #user_input["Plotting"]["Figure"]["y_multiplier"] = 1
        

        sensitivity = 3
        separation = 5

        
        simulation = theory.TS_Simulator(user_input)
        observer   = observe.OS_Simulator(user_input)
        analyzer   = analyze.Spectra_Analyzer(user_input)
        
        simulation = simulate_NIST(simulation, observer, analyzer, atmosphere, stuff, sensitivity, separation)  
  
    