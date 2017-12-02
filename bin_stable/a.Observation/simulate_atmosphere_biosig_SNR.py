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

Simulate an exoplanet atmosphere with/without biosignature, then figure out the 
SNR_bio for the biosignature detection.

"""
import os
import sys
import numpy as np
import random
from scipy import stats

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Main.simulation.observed_spectra_simulator as observe
import SEAS_Main.simulation.spectra_analyzer as analyze
from SEAS_Main.observation_effects.noise import Photon_Noise

import SEAS_Utils as utils
import SEAS_Utils.common_utils.data_plotter as plotter
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.constants import *

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
    s.Reference_Transit_Signal = s.load_atmosphere_geometry_model(result="Height")
    s.Bio_Transit_Signal = s.load_atmosphere_geometry_model(bio=bio_enable,result="Height")
    
    
    #nu,ref_trans = s.nu, s.Reference_Transit_Signal
    #nu,bio_trans = s.nu, s.Bio_Transit_Signal
    # calculate observed transmission spectra
    nu,ref_trans = o.calculate_convolve(s.nu, s.Reference_Transit_Signal)
    nu,bio_trans = o.calculate_convolve(s.nu, s.Bio_Transit_Signal)
    
    # analyze the spectra
    #s.nu_window = a.spectra_window(nu,ref_trans,"T",0.3, 100.,s.min_signal)
    s.nu_window,thres = a.new_spectra_window(nu,ref_trans,0.3,100,"nu",thres=True)

    s.user_input["Star"]["R_Star"] = 0.6
    s.user_input["Star"]["T"] = 4000.
    s.user_input["Planet"]["R_Planet"] = 1
    s.user_input["Planet"]["R_Atmosphere"] = 40*1000
    s.user_input["Telescope"]["Aperture"] = 6.5
    s.user_input["Telescope"]["Distance"] = 10*Psec
    s.user_input["Telescope"]["Duration"] = 100*3600
    s.user_input["Telescope"]["Quantum_Efficiency"] = 0.3
    s.user_input["Telescope"]["min_wavelength"] = 1
    s.user_input["Telescope"]["max_wavelength"] = 25
    s.user_input["Observation_Effects"]["Noise"]["Multiplier"] = 1.2
    s.user_input["Observation_Effects"]["bin_exponent"] = 3/2.
    s.user_input["Observation_Effects"]["bin_width"] = 0.01
    
    
    noise = Photon_Noise(s.user_input)

    diff = np.array(s.Bio_Transit_Signal)-np.array(s.Reference_Transit_Signal)
    difference = np.array(bio_trans)-np.array(ref_trans)


    bin_edges, bin_width, bin_centers = noise.determine_bin()
    
    bin_means1, bin_edges1, binnumber1 = stats.binned_statistic(10000./nu[::-1], difference[::-1], bins=bin_edges)

    bin_means2, bin_edges2, binnumber2 = stats.binned_statistic(10000./s.nu[::-1], diff[::-1], bins=bin_edges)

    new_bin_means = []
    for i,info in enumerate(bin_means1):
        if float(info) < 0:
            new_bin_means.append(0)
        else:
            new_bin_means.append(info)
    bin_means1 = np.array(new_bin_means)
    new_bin_means = []
    for i,info in enumerate(bin_means2):
        if float(info) < 0:
            new_bin_means.append(0)
        else:
            new_bin_means.append(info)
    bin_means2 = np.array(new_bin_means)
    
    signal1, photon_noise1, SNR1 = noise.calculate_noise(bin_means1)
    signal2, photon_noise2, SNR2 = noise.calculate_noise(bin_means2)
    
    print signal1,photon_noise1
    
    s.user_input["Plotting"]["Figure"]["Title"] = "Effective Height for Simulated Earth Atmosphere with traces of %s at %s ppm"%(bio_molecule,bio_abundance*10**6)
    s.user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    s.user_input["Plotting"]["Figure"]["y_label"] = r"Effective Height (m)"    
    s.user_input["Plotting"]["Figure"]["y_multiplier"] = 1
    sim_plot = plotter.Simulation_Plotter(s.user_input)
    
    sim_plot.plot_xy(bin_centers, SNR1)
    sim_plot.plot_xy(bin_centers, SNR2)
    
    #sim_plot.plot_xy(s.nu,diff)
    #sim_plot.plot_xy(nu,difference)
    #sim_plot.plot_xy(10000./np.array(bin_centers),bin_means,"o")
    #sim_plot.plot_bin(bin_centers, bin_means, photon_noise)
    
    #plt_ref_1 = sim_plot.plot_xy(nu,bio_trans,"with_bio")
    #plt_ref_2 = sim_plot.plot_xy(nu,ref_trans,"ref.")
    
    #sim_plot.plot_window(s.nu_window,"k", 0.2)
    #sim_plot.set_legend([plt_ref_1, plt_ref_2])
    
    if utils.to_bool(s.user_input["Save"]["Plot"]["save"]):
        sim_plot.save_plot()
    else:
        sim_plot.show_plot()
    
    return s


if __name__ == "__main__":
    
    
    Timer = simple_timer()
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename = "Test_Earth"
    Filename1 = "Test_Earth_with_biosig"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "earth.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "earth.txt"

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
    user_input["Save"]["Window"]["path"] = "../../output/Simple_Atmosphere_Window"
    user_input["Save"]["Window"]["name"] = "%s_Window_A1000_S100.txt"%Filename1
    
    
    info = NIST_Smile_List()
    molecule_smiles = info[0]
    Bio_Molecule = random.choice(molecule_smiles)
    
    Bio_Molecule = "CSC"
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["enable"] = True
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["data_type"] = "NIST"
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["molecule"] = Bio_Molecule
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["abundance"] = 1*10**-6
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["is_smile"] = True
    
    user_input["Atmosphere_Effects"]["Overlay"]["enable"] = True
    
    user_input["Save"]["Plot"] = {}
    user_input["Save"]["Plot"]["save"] = False    
    user_input["Save"]["Plot"]["path"] = "../../output/Plot_Result"
    user_input["Save"]["Plot"]["name"] = "%s_Plot.png"%Filename1
    
    simulation = theory.TS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    simulation = simulate_NIST(simulation, observer, analyzer)         
    