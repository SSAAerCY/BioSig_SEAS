"""
Need to have simulation_db to work

"""


import os
import sys
import numpy as np

from openpyxl import load_workbook, Workbook

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils as utils
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data, TP_Profile_Data

from SEAS_Utils.common_utils.timer import simple_timer
import SEAS_Utils.common_utils.configurable as config
import SEAS_Main.simulation.transmission_spectra_simulator as theory
import SEAS_Aux.atmosphere_processes.TP_profile_generator as TPgen
import SEAS_Aux.atmosphere_processes.mixing_ratio_generator as mix

from SEAS_Utils.common_utils.DIRs import HITRAN_Molecule_List

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)




def simulate(s):
    
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
    
    
    # this "H2" toggle is temporary and ...
    # CIA_Enable == True or CIA_Enable.lower() == "true"
    # this is also buggy for cases where there's no biosig gas considered but CIA enabled
    CIA_enable = (utils.to_bool(s.user_input["Atmosphere_Effects"]["CIA"]["enable"])==True and "H2" in s.bio_normalized_molecules)

    if CIA_enable:
        s.CIA_File, s.CIA_Data = s.load_CIA(["H2"])
        s.normalized_CIA = s.interpolate_CIA()
    
    if "H2" not in s.normalized_molecules and "H2" in s.bio_normalized_molecules:
        s.Reference_Transit_Signal = s.load_atmosphere_geometry_model(CIA=False)
    else:
        s.Reference_Transit_Signal = s.load_atmosphere_geometry_model(CIA=CIA_enable)
    
    # load rayleigh scattering
    Rayleigh_enable = utils.to_bool(s.user_input["Atmosphere_Effects"]["Rayleigh"]["enable"])
    if Rayleigh_enable:
        s.normalized_rayleigh = s.load_rayleigh_scattering()    
    
    
    # calculate transmission spectra
    nu,ref_trans = s.calculate_convolve(s.Reference_Transit_Signal)
    
    s.Bio_Transi_Signal = s.load_atmosphere_geometry_model(bio=bio_enable,CIA=CIA_enable)
    nu,bio_trans = s.calculate_convolve(s.Bio_Transi_Signal)
    
    # window calculation?
    s.nu_window = s.spectra_window(nu,ref_trans,"T",0.3,100.)
    
    # determine detection
    s.detection, s.detection_window = s.analyze_spectra_detection(nu,ref_trans,bio_trans)
    
    
    print "calc time", s.Timer.elapse()
    
    return s.detection, s.nu_window, s.detection_window






