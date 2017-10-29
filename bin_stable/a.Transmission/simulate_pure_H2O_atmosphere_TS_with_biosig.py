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
    s.Reference_Transit_Signal = s.load_atmosphere_geometry_model()
    s.Bio_Transit_Signal = s.load_atmosphere_geometry_model(bio=bio_enable)
    
    # calculate observed transmission spectra
    nu,ref_trans = o.calculate_convolve(s.nu, s.Reference_Transit_Signal)
    nu,bio_trans = o.calculate_convolve(s.nu, s.Bio_Transit_Signal)
    
    # analyze the spectra
    s.nu_window = a.spectra_window(nu,ref_trans,"T",0.3, 100.,s.min_signal)
    
    
    user_input["Plotting"]["Figure"]["Title"] = "Transit Signal and Atmospheric Window for Simulated Earth Atmosphere with traces of %s at %s ppm"%(bio_molecule,bio_abundance*10**6)
    user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    
    
    sim_plot = plotter.Simulation_Plotter(s.user_input)
    
    plt_ref_1 = sim_plot.plot_xy(nu,bio_trans,"with_bio")
    plt_ref_2 = sim_plot.plot_xy(nu,ref_trans,"ref.")
    
    sim_plot.plot_window(s.nu_window,"k", 0.2)
    sim_plot.set_legend([plt_ref_1, plt_ref_2])
    
    if utils.to_bool(user_input["Save"]["Plot"]["save"]):
        sim_plot.save_plot()
    else:
        sim_plot.show_plot()
    
    return s


if __name__ == "__main__":
    
    
    Timer = simple_timer()
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    Filename = "Test_Water"
    
    user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
    user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
    user_input["Simulation_Control"]["TP_Profile_Name"]     = "TP_iso_H2O_26576258.txt"
    user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = "MR_iso_H2O_26576258.txt"

    user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Filename
    
    
    info = NIST_Smile_List()
    molecule_smiles = info[0]
    Bio_Molecule = random.choice(molecule_smiles)
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["enable"] = True
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["data_type"] = "NIST"
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["molecule"] = Bio_Molecule
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["abundance"] = 1*10**-6
    user_input["Atmosphere_Effects"]["Bio_Molecule"]["is_smile"] = True
    
    user_input["Atmosphere_Effects"]["Overlay"]["enable"] = True
    
    
    simulation = theory.TS_Simulator(user_input)
    observer   = observe.OS_Simulator(user_input)
    analyzer   = analyze.Spectra_Analyzer(user_input)
    
    simulation = simulate_NIST(simulation, observer, analyzer)         
    