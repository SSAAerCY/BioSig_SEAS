"""

construct the schema to make an anoxic atmosphere.


When proposing a new type of atmosphere, should also have function to generate
tp profile and mixing ratio.


atmosphere contains mostly 
    abundance CH4, NH3, H2, N2O 
    and trace of ... $C_NH_N$ C2H2, C2H4, C2H6?
    tiny trace amount of H2O, CO2, O2 and O3


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
from SEAS_Utils.common_utils.DIRs import *
from SEAS_Utils.common_utils.constants import *
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
from SEAS_Utils.common_utils.timer import simple_timer


import SEAS_Aux.atmosphere_processes.TP_profile_generator as TPgen
import SEAS_Aux.atmosphere_processes.mixing_ratio_generator as mix

def simulate_anoxic(s,o,a):

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
    
    if "H2" in s.normalized_molecules:
        s.user_input["Atmosphere_Effects"]["CIA"]["enable"] = True
        s.CIA_File, s.CIA_Data = s.load_CIA(["H2"])
        s.normalized_CIA = s.interpolate_CIA()
        s.Reference_Transit_Signal = s.load_atmosphere_geometry_model(CIA=True,result="Height")
        s.Bio_Transit_Signal = s.load_atmosphere_geometry_model(bio=bio_enable,CIA=True,result="Height")
    else:
        s.Reference_Transit_Signal = s.load_atmosphere_geometry_model(result="Height")
        s.Bio_Transit_Signal = s.load_atmosphere_geometry_model(bio=bio_enable,result="Height")
    
    # calculate observed transmission spectra
    nu,ref_trans = o.calculate_convolve(s.nu, s.Reference_Transit_Signal)
    nu,bio_trans = o.calculate_convolve(s.nu, s.Bio_Transit_Signal)
    
    # analyze the spectra
    #s.nu_window = a.spectra_window(nu,ref_trans,"T",0.3, 100.,s.min_signal)
    
    Window_Threshold = utils.to_float(s.user_input["Observation_Effects"]["Atmosphere_Window"]["threshold"])
    Window_Span = utils.to_float(s.user_input["Observation_Effects"]["Atmosphere_Window"]["span"])
    s.nu_window,thres = a.new_spectra_window(nu,ref_trans,Window_Threshold,Window_Span,"nu",thres=True)
    

    
    
    s.user_input["Plotting"]["Figure"]["Title"] = "Transit Signal and Atmospheric Window for Simulated %s with traces of %s at %s ppm"%(s.user_input["Title_name"],bio_molecule,bio_abundance*10**6)
    s.user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    s.user_input["Plotting"]["Figure"]["y_label"] = r"Transit Signal (ppm)"    

    s.user_input["Plotting"]["Figure"]["Title"] = "Atmosphere Height and Atmospheric Window for Simulated %s with traces of %s at %s ppm"%(s.user_input["Title_name"],bio_molecule,bio_abundance*10**6)
    s.user_input["Plotting"]["Figure"]["x_label"] = r'Wavelength ($\mu m$)'
    s.user_input["Plotting"]["Figure"]["y_label"] = r"Atmosphere Height (m)"   
    
    
    s.user_input["Star"]["R_Star"] = 0.6
    s.user_input["Star"]["T"] = 4000.
    s.user_input["Planet"]["R_Planet"] = 1
    #s.user_input["Planet"]["R_Atmosphere"] = 40*1000
    s.user_input["Telescope"]["Aperture"] = 6.5
    s.user_input["Telescope"]["Distance"] = 10*Psec
    s.user_input["Telescope"]["Duration"] = 100*3600
    s.user_input["Telescope"]["Quantum_Efficiency"] = 0.3
    s.user_input["Telescope"]["min_wavelength"] = 1
    s.user_input["Telescope"]["max_wavelength"] = 25
    s.user_input["Observation_Effects"]["Noise"]["Multiplier"] = 1.2
    s.user_input["Observation_Effects"]["bin_exponent"] = 3/2.
    s.user_input["Observation_Effects"]["bin_width"] = 0.01


    s.user_input["Plotting"]["Figure"]["y_multiplier"] = 1
    
    noise = Photon_Noise(s.user_input)

    diff = np.array(s.Bio_Transit_Signal)-np.array(s.Reference_Transit_Signal)
    print diff

    bin_edges, bin_width, bin_centers = noise.determine_bin()    
    
    bin_means_bio, bin_edges_bio, binnumber_bio = stats.binned_statistic(10000./s.nu[::-1], s.Bio_Transit_Signal[::-1], bins=bin_edges)
    
    bin_means, bin_edges, binnumber = stats.binned_statistic(10000./s.nu[::-1], diff[::-1], bins=bin_edges)

    new_bin_means = []
    for i,info in enumerate(bin_means):
        if float(info) < 0:
            new_bin_means.append(0)
        else:
            new_bin_means.append(info)
    bin_means = np.array(new_bin_means)    


    signal, photon_noise, SNR = noise.calculate_noise(bin_means)

    sim_plot = plotter.Simulation_Plotter(s.user_input)
    
    sim_plot.plot_hline(thres[0][0],[thres[0][1],thres[0][2]],"Threshold")
    sim_plot.plot_hline(thres[1][0],[thres[1][1],thres[1][2]],"Threshold")
    sim_plot.plot_hline(thres[2][0],[thres[2][1],thres[2][2]],"Threshold")
    
    sim_plot.plot_bin(bin_centers, bin_means_bio, photon_noise)
    #sim_plot.plot_xy(bin_centers,bin_means_bio)
    
    
    plt_ref_1 = sim_plot.plot_xy(nu,bio_trans,"with_bio")
    plt_ref_2 = sim_plot.plot_xy(nu,ref_trans,"ref.")
    
    
    
    sim_plot.plot_window(s.nu_window,"k", 0.2)
    #sim_plot.set_legend([plt_ref_1, plt_ref_2])
    
    if utils.to_bool(s.user_input["Save"]["Plot"]["save"]):
        sim_plot.save_plot()
    else:
        sim_plot.show_plot()
    

    return s

def simulation_prep(name, atmosphere,type="File"):
    
    
    if type == "File":
        Filename = name
        Molecules = atmosphere["Molecules"]
        MRFile = atmosphere["Mixing_Ratio_File"]
        TPFile = atmosphere["TP_File"]
        Window_Threshold = 0.3
        Window_Span = 100
        
        return Filename, Filename, Molecules, TPFile, MRFile, Window_Threshold,Window_Span


    elif type == "Half":
        Type                = atmosphere["Type"]
        Molecules           = atmosphere["Molecules"]
        Mixing_Raio         = atmosphere["Ratio"]
        Filler              = atmosphere["Filler"]
        TP_Name             = atmosphere["TP_File"]
        Window_Threshold    = atmosphere["Window_Threshold"]
        Window_Span         = atmosphere["Window_Span"]
        Simple_Filename = "_".join([Type[:3],"_".join(Molecules),name])
        Complex_Filename = "_".join([Type[:3],"_".join(Molecules),str(hash(str(atmosphere.values()))%10**8)])
        
        
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
        
        MR_Name = "MR_%s.txt"%Complex_Filename
        if not os.path.isfile(os.path.join(Mixing_Ratio_Data,MR_Name)):
            simulator = mix.mixing_ratio_generator(ratio_input,
                                                   filler=True,
                                                   filler_molecule=Filler,
                                                   pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                                   name=MR_Name,
                                                   overwrite=True)
            simulator.generate()
            simulator.save()

    else:
        Type                = atmosphere["Type"]
        Surface_Temperature = atmosphere["Surface_Temperature"]
        Molecules           = atmosphere["Molecules"]
        Mixing_Raio         = atmosphere["Ratio"]
        Simple_Filename = "_".join([Type[:3],"_".join(Molecules),name])
        Complex_Filename = "_".join([Type[:3],"_".join(Molecules),str(hash(str(atmosphere.values()))%10**8)])
        Window_Threshold = 0.3
        Window_Span = 100
                
        try:
            Filler = atmosphere["Filler"]
        except:
            Filler = "N2"
        Simple_Filename = "_".join([Type[:3],"_".join(Molecules),name])
        Complex_Filename = "_".join([Type[:3],"_".join(Molecules),str(hash(str(atmosphere.values()))%10**8)])
        
        if Type == "isothermal":
            TP_input = config.Configuration("../../bin_stable/TP_Profile/TP_selection.cfg")["Test Isothermal Atmosphere"]
            TP_Name = "TP_%s.txt"%Complex_Filename
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
        
        MR_Name = "MR_%s.txt"%Complex_Filename
        if not os.path.isfile(os.path.join(Mixing_Ratio_Data,MR_Name)):
            simulator = mix.mixing_ratio_generator(ratio_input,
                                                   filler=True,
                                                   filler_molecule=Filler,
                                                   pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05],
                                                   name=MR_Name,
                                                   overwrite=True)
            simulator.generate()
            simulator.save()

    return Simple_Filename, Complex_Filename, Molecules,TP_Name, MR_Name, Window_Threshold, Window_Span

def setup_simulation():

    Timer = simple_timer()
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")

    
    atmo_input = config.Configuration("../../bin_stable/a.Main/atmosphere_prototype.cfg")
    atmosphere_types = atmo_input["Proposed"]
    
    for i, atmosphere in enumerate(atmosphere_types):
    
        Simple_Filename, Complex_Filename, Molecules, TP_Name, MR_Name, Window_Threshold, Window_Span = simulation_prep(atmosphere, atmosphere_types[atmosphere],"Half")
        print Simple_Filename, Molecules, TP_Name, MR_Name
    
        user_input["Title_name"] = atmosphere
    
        user_input["Simulation_Control"]["DB_DIR"]              = "Simulation_Band"
        user_input["Simulation_Control"]["DB_Name"]             = None#"cross_sec_Example.db"
        user_input["Simulation_Control"]["TP_Profile_Name"]     = TP_Name
        user_input["Simulation_Control"]["Mixing_Ratio_Name"]   = MR_Name
    
        user_input["Save"]["Intermediate_Data"]["cross_section_savename"] = "%s_Cross_Section.npy"%Simple_Filename
        
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
        
        user_input["Observation_Effects"]["Atmosphere_Window"]["enable"] = True
        user_input["Observation_Effects"]["Atmosphere_Window"]["threshold"] = Window_Threshold
        user_input["Observation_Effects"]["Atmosphere_Window"]["span"] = Window_Span
        
        simulation = theory.TS_Simulator(user_input)
        observer   = observe.OS_Simulator(user_input)
        analyzer   = analyze.Spectra_Analyzer(user_input)
        
        
        simulation = simulate_anoxic(simulation, observer, analyzer)      



if __name__ == "__main__":
    
    setup_simulation()
    
    
    
    
    
    
    
    