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

Functions related to calculate and using of transmission spectra

for 0.8, need a full scale conversion of all list into dicts
instead of normalized_xxx, let's have a dict with pressure_layers as keys and relevent data as data

"""
import os
import sys
import numpy as np
import time
from scipy import interpolate
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)

import SEAS_Main.atmosphere_geometry
import SEAS_Main.atmosphere_property
import SEAS_Main.analysis
from SEAS_Main.atmosphere_effects.biosig_molecule import load_NIST_spectra, biosig_interpolate


from SEAS_Aux.calculation.interpolation import interpolate1d
import SEAS_Aux.calculation.astrophysics as calc 
import SEAS_Aux.cross_section.hapi as hp

from SEAS_Utils.common_utils.constants import *
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.DIRs import TP_Profile_Data, Mixing_Ratio_Data, molecule_info, DB_DIR,Intermediate_DIR, HITRAN_CIA
from SEAS_Utils.common_utils.data_loader import two_column_file_loader,multi_column_file_loader, json_loader, molecule_cross_section_loader2
from SEAS_Utils.common_utils.data_saver import check_file_exist, check_path_exist
import SEAS_Utils.common_utils.db_management2 as dbm







class TS_Simulator():
    """
    This is the main class for the planet transmission spectra simulation.
    
    simulate a theoratical spectra of the planet transmission. 
    does not have any observational effects
    
    mixing ratio and TP profiles are pre-calculated data.
    simulating these should not be the task for TS_Simulator 
    
    The simulation will also not determine biosignatures or 
    try to mod the mixing ratio list, etc. All of these should be 
    pre-determined when generating the mixing ratio file.
    """
    
    def __init__(self,user_input):
        
        self.user_input = user_input
        
        self.DB_DIR = os.path.join(DB_DIR,self.user_input["Simulation_Control"]["DB_DIR"])
        
        self.T_grid = [float(x) for x in user_input["Simulation_Control"]["T_Grid"]]
        self.P_grid = [float(x) for x in user_input["Simulation_Control"]["P_Grid"]]

        self.load_spectral_properties()

    def simulate_example(self):
        
        self.Timer = simple_timer(4)
        
        #create a base flat spectra based on planetary parameters
        self.Surface_g, self.Base_TS_Value = self.load_astrophysical_properties()
        
        # normalized pressure directly correspond to atmosphere layers.
        self.normalized_pressure = self.load_atmosphere_pressure_layers()
        
        # load mixing ration files and determine what molecules are added to the simulation
        # acquire the mixing ratio at each pressure (need to be interpolated to normalized pressure)
        self.normalized_molecules, self.MR_Pressure, self.molecule_MR = self.load_mixing_ratio()
        self.normalized_abundance = self.interpolate_mixing_ratio()
        
        # load temperature pressure profile
        self.TP_Pressure, self.TP_Temperature = self.load_TP_profile()        
        self.normalized_temperature = self.interpolate_TP_profile()

        # calculate the scale height for each layer of the atmosphere
        self.normalized_scale_height = self.calculate_scale_height()
   
        # load molecular cross section for main constituent of the atmosphere
        # will check first to see if the molecule is in the database 
        self.cross_db = self.check_molecules()
        self.nu, self.normalized_cross_section = self.load_molecule_cross_section()
        
        
        self.normalized_rayleigh = self.load_rayleigh_scattering()
        
        
        
        print "load time", self.Timer.elapse()
    
        self.Transit_Signal = self.load_atmosphere_geometry_model()
        nu,trans = self.calculate_convolve(self.Transit_Signal)
        
        self.plot_result(nu,trans)
        #return self.Transit_Signal   
    
    def simulate_window(self):
        
        
        self.Timer = simple_timer(4)
        
        #create a base flat spectra based on planetary parameters
        self.Surface_g, self.Base_TS_Value = self.load_astrophysical_properties()
        
        # normalized pressure directly correspond to atmosphere layers.
        self.normalized_pressure = self.load_atmosphere_pressure_layers()
        
        # load mixing ration files and determine what molecules are added to the simulation
        # acquire the mixing ratio at each pressure (need to be interpolated to normalized pressure)
        self.normalized_molecules, self.MR_Pressure, self.molecule_MR = self.load_mixing_ratio()
        self.normalized_abundance = self.interpolate_mixing_ratio()
        
        # load temperature pressure profile
        self.TP_Pressure, self.TP_Temperature = self.load_TP_profile()        
        self.normalized_temperature = self.interpolate_TP_profile()

        # calculate the scale height for each layer of the atmosphere
        self.normalized_scale_height = self.calculate_scale_height()
   
        # load molecular cross section for main constituent of the atmosphere
        # will check first to see if the molecule is in the database 
        self.cross_db = self.check_molecules()
        self.nu, self.normalized_cross_section = self.load_molecule_cross_section()
        
        self.normalized_rayleigh = self.load_rayleigh_scattering()
        
        print "load time", self.Timer.elapse()

        self.effects = self.load_effects()

        # determine atmosphereic effects to add in
        
        self.Transit_Signal = self.load_atmosphere_geometry_model()
        nu,trans = self.calculate_convolve(self.Transit_Signal)
        
        self.absorption = self.load_atmosphere_geometry_model2()
        nu,absorp = self.calculate_convolve(self.absorption)
        print "calc time", self.Timer.elapse()
        
        self.nu_window = self.spectra_window(nu,trans,"T",0.3, 100.)
       
        """
        plt.title("Absorption and Atmospheric Window for Simulated Atmosphere of %s"%"_".join(self.normalized_molecules))
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel("absorption")    

        fig_size = plt.rcParams["figure.figsize"]
        print fig_size
        plt.rcParams["figure.figsize"] = [20,6]

        ax = plt.gca()
        ax.set_xscale('log')
        #ax.set_yscale("log")
        plt.tick_params(axis='x', which='minor')
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))          

        for k in self.nu_window:
            up,down = 10000./k[1],10000./k[0]
            plt.axvspan(up,down,facecolor="k",alpha=0.2)

        #plt.axhline(y=8520, linewidth=2, color = 'k')


        plt.plot(10000./nu, absorp,color="r")
        plt.plot(10000./np.array(self.stuff),np.ones(len(self.stuff))*1000,".")
       
        plt.show()
        """
        
        self.plot_result(nu,trans)
        #self.plot_result(self.nu,self.Transit_Signal)
        
        
        return self.Transit_Signal

    def simulate_CIA(self):
        
        self.Timer = simple_timer(4)
        
        #create a base flat spectra based on planetary parameters
        self.Surface_g, self.Base_TS_Value = self.load_astrophysical_properties()
        
        # normalized pressure directly correspond to atmosphere layers.
        self.normalized_pressure = self.load_atmosphere_pressure_layers()
        
        # load mixing ration files and determine what molecules are added to the simulation
        # acquire the mixing ratio at each pressure (need to be interpolated to normalized pressure)
        self.normalized_molecules, self.MR_Pressure, self.molecule_MR = self.load_mixing_ratio()
        self.normalized_abundance = self.interpolate_mixing_ratio()
        
        # load temperature pressure profile
        self.TP_Pressure, self.TP_Temperature = self.load_TP_profile()        
        self.normalized_temperature = self.interpolate_TP_profile()

        # calculate the scale height for each layer of the atmosphere
        self.normalized_scale_height = self.calculate_scale_height()

        self.normalized_molecules, self.MR_Pressure, self.molecule_MR = self.load_mixing_ratio()
        self.normalized_abundance = self.interpolate_mixing_ratio()
        
        # load temperature pressure profile
        self.TP_Pressure, self.TP_Temperature = self.load_TP_profile()        
        self.normalized_temperature = self.interpolate_TP_profile()

        # calculate the scale height for each layer of the atmosphere
        self.normalized_scale_height = self.calculate_scale_height()        

        self.cross_db = self.check_molecules()
        self.nu, self.normalized_cross_section = self.load_molecule_cross_section()
        
        self.normalized_rayleigh = self.load_rayleigh_scattering()

        self.CIA_File, self.CIA_Data = self.load_CIA(["H2"])
        self.normalized_CIA = self.interpolate_CIA()
        
        self.Transit_Signal_CIA = self.load_atmosphere_geometry_model(CIA=True)
        
        self.Transit_Signal = self.load_atmosphere_geometry_model()
        
        ax = plt.gca()
        ax.set_xscale('log')
        #ax.set_yscale("log")
        plt.tick_params(axis='x', which='minor')
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))           
        

        plt.title("Transit Signal and Atmospheric Window for Simulated Atmosphere of %s"%"_".join(self.normalized_molecules))
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel("Transit Signal (ppm)")  
    
       
        plt1, = plt.plot(10000./self.nu, 1000000*self.Transit_Signal, label="Molecular")
        plt2, = plt.plot(10000./self.nu, 1000000*self.Transit_Signal_CIA, label="Molecular+CIA")
        
        plt.legend(handles=[plt1,plt2])
        plt.show()

    def load_astrophysical_properties(self):
        """
        There is going to be more development here in the future
        For now it's just loading the star and planet radius, then calculate the baseline
        for the simulation
        """
        
        self.R_Star          = float(self.user_input["Star"]["R_Star"])*R_Sun
        self.R_planet        = float(self.user_input["Planet"]["R_Planet"])*R_Earth
        self.M_planet        = float(self.user_input["Planet"]["M_Planet"])*M_Earth
        
        Surface_g = calc.calc_SurfaceG(self.M_planet, self.R_planet)
        Base_TS_Value   = (self.R_planet/self.R_Star)**2
        
        return Surface_g, Base_TS_Value

    def load_spectral_properties(self):
        """
        cross section query should take this in as input
        interpolations also need to consider this
        """
        
        self.numin = float(self.user_input["Spectra"]["Numin"])
        self.numax = float(self.user_input["Spectra"]["Numax"])
        
    def load_atmosphere_pressure_layers(self):
        
        Surface_Pressure = float(self.user_input["Atmosphere"]["Surface_Pressure"])
        P_Cutoff         = float(self.user_input["Atmosphere"]["P_Cut_Off"])

        self.P_grid = []
        P = Surface_Pressure
        while P> P_Cutoff:
            self.P_grid.append(float("%.3g"%P))
            P = P*np.e**-(1./float(self.user_input["Atmosphere"]["Sub_Layers"]))
            
        return self.P_grid
    
    def load_mixing_ratio(self):

        if self.user_input["Simulation_Control"]["Mixing_Ratio"] == "File":
            filename = os.path.join(Mixing_Ratio_Data,self.user_input["Simulation_Control"]["Mixing_Ratio_Name"])
            data = multi_column_file_loader(filename,type="mixed")

            # mixing ratio are in %, so need to divide by 100
            molecules, molecule_MR = [],[]
            for i in data[1:]:
                molecules.append(i[0])
                molecule_MR.append([float(x)/100. for x in i[1:]])

            MR_pressure = [float(x) for x in data[0][1:]]
            
            return molecules, MR_pressure, molecule_MR
            
        else:
            print "Mixing_Ratio currently only takes files"
            sys.exit()
        
    def interpolate_mixing_ratio(self):
        
        molecule_MR = []
        x = np.log(self.MR_Pressure)
        X = np.log(self.normalized_pressure) 
        for y in self.molecule_MR:
            molecule_MR.append(interpolate1d(x,y,X))
        
        normalized_abundance = np.array(molecule_MR).T
        
        return normalized_abundance

    def load_TP_profile(self):

        if self.user_input["Simulation_Control"]["TP_Profile"] == "File":
            filename = os.path.join(TP_Profile_Data,self.user_input["Simulation_Control"]["TP_Profile_Name"])
            return two_column_file_loader(filename)
        else:
            print "Mixing_Ratio currently only takes files"
            sys.exit()
            
    def interpolate_TP_profile(self):
        
        x = np.log(self.TP_Pressure)
        y = self.TP_Temperature
        X = np.log(self.normalized_pressure)
        
        normalized_temperature = interpolate1d(x,y,X)
        return normalized_temperature    

    def load_rayleigh_scattering(self):
        """
        currently not caring about biosignature molecule rayleigh?
        """
        
        Rayleigh_array = []
        for molecule in self.normalized_molecules:
            Rayleigh_array.append(calc.calc_rayleigh(molecule, self.nu))
        
        return Rayleigh_array
            
    def interpolate_rayleigh_scattering(self):
        pass
    
    def load_CIA(self,molecules=None):
        
        if molecules == None:
            molecules = self.normalized_molecules
        
        CIA_Enable = self.user_input["Atmosphere_Effects"]["CIA"]["enable"]
        if  CIA_Enable == True or CIA_Enable.lower() == "true":
            
            import SEAS_Main.atmosphere_effects.collision_induced_absorption as cia 
            
            data_file = cia.select_molecule_cia(molecules)
            data = {}
            for filename in data_file:
                processor = cia.HITRAN_CIA_data_processor(HITRAN_CIA, filename)
                temp,nu,xsec = processor.load()
                data[filename] = {}
                
                data[filename]["nu"] = nu
                for i,T in enumerate(temp):
                    data[filename][T] = {}
                    
                    data[filename][T]["xsec"] = xsec[i]
            return data_file,data
        else:
            return {},data

    def interpolate_CIA(self):
        
        print self.normalized_temperature
        
        normalized_CIA = []
        for i in self.CIA_File:
            x1 = np.array(self.CIA_Data[i]["nu"],dtype=float)
            y1 = np.array(self.CIA_Data[i][300.]["xsec"],dtype=float)
            xsec = biosig_interpolate(x1,y1,self.nu,"C")
                

            current_CIA = []
            for j in range(len(self.normalized_temperature)):
                current_CIA.append(xsec)
                
            normalized_CIA.append(current_CIA)
        
        # the temperature need to be interpolated first
        # then the nu need to be interpolated and augmented
        # need db information on how to slice
        return normalized_CIA

    def calculate_scale_height(self):
        
        # load molecular weight for each molecule in the list
        molecular_weight_list = []
        molecule_reference = json_loader(os.path.join(molecule_info,"Simulation_Molecule.json"))
        for molecule in self.normalized_molecules:
            molecular_weight_list.append(molecule_reference[molecule]["Molecular_Weight"])
    
        normalized_scale_height = []
        for i,abundance_per_layer in enumerate(self.normalized_abundance):
            mean_mw  = calc.calc_MeanMolWeight(abundance_per_layer, molecular_weight_list)
            scale_height = calc.calc_H(self.normalized_temperature[i], mean_mw*mH, self.Surface_g)
            normalized_scale_height.append(scale_height)

        return normalized_scale_height
  
    def check_molecules(self):
        """
        need some rework for this function to accomodate data from different database
        also how do we distinguish "main" molecule and "biosig" molecule?
        or should we have a smart way of classifying them?
        say... data ranking system? data selection preferences?
        
        for now only from simulation db
    
        """
        if self.user_input["Simulation_Control"]["DB_Name"] == None:
            return


        
        kwargs = {"dir"        :self.DB_DIR,
                  "db_name"    :self.user_input["Simulation_Control"]["DB_Name"],
                  "user"       :self.user_input["User"]["username"],
                  "DEBUG"      :False,
                  "REMOVE"     :False,
                  "BACKUP"     :False,
                  "OVERWRITE"  :False}
    
        cross_db = dbm.database(**kwargs)
        cross_db.access_db()

        for molecule in self.normalized_molecules:
            if cross_db.check_table_exist(molecule) == False:
                print "molecule: %s not in database, simulation_terminated"%molecule
                print "also need to check if database exists"
                sys.exit()
        
        return cross_db
                
    def load_molecule_cross_section(self):
        """
        a smarter interpolation is needed here 
        a smarter saver is needed too with dynamic filename
        """
    
        normalized_temperature = self.normalized_temperature
        normalized_pressure = self.normalized_pressure
    
        normalized_cross_section = []
        number_of_layers = len(normalized_pressure)   


        if self.user_input["Save"]["Intermediate_Data"]["cross_section_saved"] == "true":
            savepath = os.path.join(Intermediate_DIR,self.user_input["Save"]["Intermediate_Data"]["cross_section_savepath"])
            savename = os.path.join(savepath,self.user_input["Save"]["Intermediate_Data"]["cross_section_savename"])
            if check_file_exist(savename):
                print "Cross Section Loaded from Save"
                return np.load(savename)      

        print "Cross Section Interpolated from Database"
        for molecule in self.normalized_molecules:
            print molecule,
            if molecule == "N2" or molecule == "He":
                # hard coded for now. need to change
                # doing this means we can't just have N2 or H2 by itself
                processed_cross_section = np.zeros((len(self.normalized_temperature), len(self.normalized_temperature), 12000))
                
            else:
                nu, raw_cross_section_grid = molecule_cross_section_loader2(self.user_input, self.DB_DIR, molecule)
            
            
                print "load:%s"%self.Timer.elapse(),
                processed_cross_section = []
                for layer_cross_section in raw_cross_section_grid:
                    
                    raw_cross_section_by_wn = np.array(layer_cross_section).T
                    pro_cross_section_by_wn = []
    
                    x = self.T_grid
                    X = self.normalized_temperature           
                    for y in raw_cross_section_by_wn:
                        pro_cross_section_by_wn.append(interpolate1d(x,y,X))
    
                    # something is definately wrong here with saving too much data. we don't need a 24x24 grid. just 24x2
                    # also no need to interpolate if isothermal atmosphere
                    processed_cross_section.append(np.array(pro_cross_section_by_wn).T)
                
            normalized_cross_section.append(processed_cross_section)
            print "inte:%s"%self.Timer.elapse()

        if self.user_input["Save"]["Intermediate_Data"]["cross_section_saved"] == "true":
            
            savepath = os.path.join(Intermediate_DIR,self.user_input["Save"]["Intermediate_Data"]["cross_section_savepath"])
            check_path_exist(savepath)
            savename = os.path.join(savepath,self.user_input["Save"]["Intermediate_Data"]["cross_section_savename"])           
            
            np.save(savename, [nu,normalized_cross_section])
            print "\nInterpolated Cross Section Saved!"
        
        
        
        return nu, normalized_cross_section

    def load_bio_molecule_cross_section(self, bio_molecule, data_type):
        """
        currently constrained to one molecule
        """
        
        normalized_bio_molecule_xsec = []
        
        if bio_molecule in self.normalized_molecules:
            print "bio molecule already in molecule list, no changes"
            return [],[]
        
        if bio_molecule == "N2" or bio_molecule == "He":
            print "Invalid Bio Molecule"
            return [],[]
        
        
        # maybe add a checker to see if molecule is in database x?
        
        
        # need to check the interpolation time cost for the nist molecules... assuming very fast for now
        if data_type == "NIST":
            
            # is we chose nist data, it is assuming we're searching with smiles since conflict using formula?
            #x1,y1 = load_NIST_spectra(bio_molecule,["wn","T"],self.is_smile)
            x1,y1 = load_NIST_spectra(bio_molecule,["wn","T"],True)
            
            # trying to subtract the transmission spectra by its baseline
            # um... this.... so much more wavy. 
            y1 = np.array(y1)+(1-(np.mean(y1)+np.median(y1))/2)
            y1new = []
            for i in y1:
                if i > 1:
                    y1new.append(1)
                else:
                    y1new.append(i)
            y1 = y1new        
    
            # these params are a bit wavy ... need a formal analysis for each spectra?
            # and comparision with known spectra
            Pref = 10000.
            Tref = 300.        
            nref = Pref/(BoltK*Tref)
            lref = 0.05        
            
            # interpolation
            yinterp = biosig_interpolate(x1,y1,self.nu,"T")
            sigma = -np.log(yinterp)/(nref*lref)*10000  # multiply by a factor of 10000 due to unit conversion
        
        
            # need to know how to scale cross sections
            X = self.normalized_pressure        
            Y = self.normalized_temperature  
            normalized_bio_molecule_xsec = []
            for i in range(len(X)):
                normalized_bio_molecule_xsec.append([])
                for j in range(len(Y)):
                    normalized_bio_molecule_xsec[i].append(sigma)
            
            
        
            return self.nu, [normalized_bio_molecule_xsec]
        
        

        elif data_type == "HITRAN_Lines":
            # need to check if molecule is in database
            
            if self.user_input["Save"]["Intermediate_Data"]["bio_cross_section_saved"] == "true":
                savepath = os.path.join(Intermediate_DIR,self.user_input["Save"]["Intermediate_Data"]["bio_cross_section_savepath"])
                savename = os.path.join(savepath,self.user_input["Save"]["Intermediate_Data"]["bio_cross_section_savename"])
                if check_file_exist(savename):
                    print "Bio Cross Section Loaded from Save"
                    return np.load(savename)      
            
            nu, raw_cross_section_grid = molecule_cross_section_loader2(self.user_input, self.DB_DIR, bio_molecule)

            print "load:%s"%self.Timer.elapse(),
            processed_cross_section = []
            for layer_cross_section in raw_cross_section_grid:
                
                raw_cross_section_by_wn = np.array(layer_cross_section).T
                pro_cross_section_by_wn = []

                x = self.T_grid
                X = self.normalized_temperature           
                for y in raw_cross_section_by_wn:
                    pro_cross_section_by_wn.append(interpolate1d(x,y,X))

                # something is definately wrong here with saving too much data. we don't need a 24x24 grid. just 24x2
                # also no need to interpolate if isothermal atmosphere
                processed_cross_section.append(np.array(pro_cross_section_by_wn).T)

            normalized_bio_molecule_xsec.append(processed_cross_section)
            print "inte:%s"%self.Timer.elapse()
            
            
            if self.user_input["Save"]["Intermediate_Data"]["bio_cross_section_saved"] == "true":
                
                savepath = os.path.join(Intermediate_DIR,self.user_input["Save"]["Intermediate_Data"]["bio_cross_section_savepath"])
                check_path_exist(savepath)
                savename = os.path.join(savepath,self.user_input["Save"]["Intermediate_Data"]["bio_cross_section_savename"])           
                
                np.save(savename, [nu,normalized_bio_molecule_xsec])
                print "\nInterpolated Bio Cross Section Saved!"            
            
            
            return nu, normalized_bio_molecule_xsec
                     
    def update_mixing_ratio(self):
        pass

    def load_overlay_effects(self):
        """
        effects that are not officially included and still under test
        
        """
        filename = "../../input/absorption_data/HITRAN_Cross_Section/O3/O3_300.0_0.0_29164.0-40798.0_04.xsc"
        nu,coef = [],[]
        with open(filename,"r") as f:
            result = f.read().split("\n")
            for i in result[1:]:
                if i != "":
                    for j in i.split():
                        coef.append(j)
            
            numin = 29164.
            numax = 40798.
            
            wavmin = 10000./numax
            wavmax = 10000./numin
            
            npts = 5818
            
            wav = np.linspace(wavmin,wavmax,npts)
            
            nu = 10000./wav[::-1]
            
        return nu,coef  
    
    def interpolate_overlay_effects(self,x1,y1):
        
        xsec = biosig_interpolate(x1,y1,self.nu,"C")
        
        
        normalized_xsec = []
        X = self.normalized_pressure        
        Y = self.normalized_temperature  
        normalized_bio_molecule_xsec = []
        for i in range(len(X)):
            normalized_bio_molecule_xsec.append([])
            for j in range(len(Y)):
                normalized_bio_molecule_xsec[i].append(xsec)
                
        return [normalized_bio_molecule_xsec]

    def load_atmosphere_geometry_model(self, bio=False, CIA=False, Rayleigh=True, result="Trans"):
        # convert all the toggle into self variables.

        TotalBeams = len(self.normalized_pressure)
        
        normalized_pressure         = self.normalized_pressure
        normalized_temperature      = self.normalized_temperature
        normalized_molecules        = self.normalized_molecules
        normalized_abundance        = self.normalized_abundance
        normalized_cross_section    = self.normalized_cross_section
        normalized_scale_height     = self.normalized_scale_height     


        if Rayleigh:
            normalized_rayleigh      = self.normalized_rayleigh        

        if CIA:
            normalized_CIA           = self.normalized_CIA
            normalized_abundance_ref = {}
            for _,mol in enumerate(normalized_molecules):
                normalized_abundance_ref[mol] = np.array(normalized_abundance).T[_]     

        
        if bio:
            normalized_cross_section = self.bio_normalized_cross_section
            normalized_abundance     = self.bio_normalized_abundance
            normalized_molecules     = self.bio_normalized_molecules
            
            # assuming no rayleigh due to biosig molecules... could be wrong
            normalized_rayleigh.append(np.zeros(len(normalized_rayleigh[0])))

        if self.Overlay_enable:
            normalized_overlay = self.normalized_overlay



        Total_Tau = np.zeros(len(self.nu))
        Total_Transit_Signal = np.ones(len(self.nu))*self.Base_TS_Value
        base_layer = self.R_planet
        
        for i in range(TotalBeams):
            
            if i == 0:
                prev_layer = base_layer
                base_layer += normalized_scale_height[i]
                # skip the bottom layer? this is the beam that "touch" the surface
                # need to think about this when rounding
                continue
    
            # opacity per beam
            BeamTau = []
            prev_pathl = 0
            target_layer = base_layer        
            for j in range(TotalBeams-i):
                
                target_layer += normalized_scale_height[j+i]
                pathl = np.sin(np.arccos(base_layer/target_layer))*target_layer - prev_pathl
                prev_pathl += pathl   
                
                # opacity per chunk of the beam, this can be thought as the test tube case
                ChunkTau = []        
                for m, molecule in enumerate(normalized_molecules):        
                    
                    #weird how abundance and cross section are wired differently
                    molecular_ratio = normalized_abundance[j+i][m]
                    number_density = (normalized_pressure[j+i]/(BoltK*normalized_temperature[j+i]))*molecular_ratio
                    
                    rayleigh = normalized_rayleigh[m]*molecular_ratio
                    sigma = normalized_cross_section[m][j+i][j+i]
                    
                    # adding ozone cross section augmentation
                    if molecule == "O3" and self.Overlay_enable:
                        overlay_sigma = normalized_overlay[0][j+i][j+i]
                        sigma = sigma+overlay_sigma
                    
                    #rayleigh_scat, etc should be pre calculated
                    effects = sigma+rayleigh#+CIA+cloud
        
                    ChunkTau_Per_Molecule = number_density*(effects)*pathl*2*0.0001
        
                    if ChunkTau == []:
                        ChunkTau = ChunkTau_Per_Molecule
                    else:
                        ChunkTau += ChunkTau_Per_Molecule   

                if CIA:
                    for k,CIA_data in enumerate(normalized_CIA):
                        
                        molecule1, molecule2 = self.CIA_File[k].replace("_","-").split("-")[:2]
                        
                        
                        molecular_ratio1 = normalized_abundance_ref[molecule1][j+i]
                        molecular_ratio2 = normalized_abundance_ref[molecule2][j+i]
                        
                        number_density1 = (normalized_pressure[j+i]/(BoltK*normalized_temperature[j+i]))*molecular_ratio1
                        number_density2 = (normalized_pressure[j+i]/(BoltK*normalized_temperature[j+i]))*molecular_ratio2
                        
                        CIA_sigma = CIA_data[j+i]
                        
                        Chunk_CIA_Tau = number_density1*number_density2*CIA_sigma*pathl*2*100*100**-3*100**-3
    
                        ChunkTau += Chunk_CIA_Tau      
                

                
                        
                if BeamTau == []:
                    BeamTau = ChunkTau
                else:
                    BeamTau += ChunkTau                       
            
            if result == "Trans":
                
                BeamTrans = calc.calc_transmittance(BeamTau)  
                RingArea = (base_layer**2-prev_layer**2)/self.R_Star**2
                
                Ring_Transit_Signal = (1-BeamTrans)*RingArea
                Total_Transit_Signal += Ring_Transit_Signal
            
            elif result == "Absorp":
                
                Total_Tau += BeamTau
            
            # update to the next beam up
            prev_layer = base_layer
            base_layer += normalized_scale_height[i]        
        
        
        self.min_signal = (self.R_planet/self.R_Star)**2
        self.max_signal = ((self.R_planet+sum(normalized_scale_height))/self.R_Star)**2
        
        
        if result == "Trans":
            Raw_Transit_Signal = Total_Transit_Signal
        elif result == "Absorp":
            Raw_Transit_Signal = Total_Tau
        
        return Raw_Transit_Signal


    def calculate_convolve(self,ydata):
        #if self.user_input["Observation"]["Convolve"] == "true":
        amount = float(self.user_input["Observation"]["Convolve_Amount"])
        nu,Transit_Signal,i1,i2,slit = hp.convolveSpectrum(self.nu,ydata,SlitFunction=hp.SLIT_RECTANGULAR,Resolution=amount,AF_wing=20.0)
        
        return nu,Transit_Signal
    
    def spectra_window(self, nu, coef, type="A",threshold=200.,span=100.):
        
        if type == "A":
        
            window = []
            
            start = True
            win = [0,0]
            
            self.stuff = []
            for i,n in enumerate(nu):
                
                if coef[i] < threshold and start == True:
                    win[0] = n
                    self.stuff.append(n)
                    start = False
                if coef[i] > threshold and start == False:
                    
                    self.stuff.append(n)
                    if n>=win[0]+span:
                        win[1] = n
                        start = True
                        window.append(np.array(win))
                    else:
                        win[0] = n
                        start = True
        
        elif type == "T":
            
            window = []
            start = True
            win = [0,0]
                        
            Min = self.min_signal
            Max = max(coef)#self.max_signal
            if threshold > 1:
                threshold = 1000/threshold
            
            
            threshold = Min+(Max-Min)*threshold
            self.threshold = threshold
            
            self.stuff = []
            for i,n in enumerate(nu):
                
                if coef[i] < threshold and start == True:
                    win[0] = n
                    self.stuff.append(n)
                    start = False
                if coef[i] > threshold and start == False:
                    
                    self.stuff.append(n)
                    if n>=win[0]+span:
                        win[1] = n
                        start = True
                        window.append(np.array(win))
                    else:
                        win[0] = n
                        start = True

        if self.user_input["Save"]["Window"]["save"] == "true":
            with open(os.path.join(self.user_input["Save"]["Window"]["path"],
                                   self.user_input["Save"]["Window"]["name"]),"w") as f:
                for i in window:
                    f.write("%s-%s\n"%(i[0],i[1]))
        
        return window

    def plot_result(self,nu,coef):
        
        def close_event():
            plt.close() #timer calls this function after 3 seconds and closes the window 
        
        fig = plt.figure(figsize=(16,6))
        timer = fig.canvas.new_timer(interval = 30000) #creating a timer object and setting an interval of 3000 milliseconds
        timer.add_callback(close_event)

        
        #plt.title("Transit Signal and Atmospheric Window for Simulated Atmosphere of %s"%"_".join(self.normalized_molecules))
        plt.title("Transit Signal and Molecule Detectable Window for Simulated exo-Earth Atmosphere")
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel("Transit Signal (ppm)")    


        ax = plt.gca()
        ax.set_xscale('log')
        #ax.set_yscale("log")
        plt.tick_params(axis='x', which='minor')
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))          

        for k in self.nu_window:
            up,down = 10000./k[1],10000./k[0]
            plt.axvspan(up,down,facecolor="k",alpha=0.2)

        #plt.axhline(y=self.threshold*10**6, linewidth=2, color = 'k')


        plt.plot(10000./nu, 1000000*coef,color="r")
        
        
        #save_dir  = self.user_input["Save"]["Plot"]["path"]
        #save_name = self.user_input["Save"]["Plot"]["name"] 
        
        #plt.savefig(os.path.join(save_dir,save_name))

        #timer.start()
        plt.show()

    def analyze_spectra_detection(self,nu,trans,bio_trans,method="max"):
        """
        How to implement area under curve?
        """
        
        noise_level = 10
        comp = 2
        detection = False
        Detected = []
        
        for i in self.nu_window:
            detected = False
            reference =  trans[list(nu).index(i[0]):list(nu).index(i[1])]
            signal = bio_trans[list(nu).index(i[0]):list(nu).index(i[1])]
            
            
            # above certain ppm
            difference = max(signal-reference)*10**6
            # above certain comparision
            comparison = max((signal-self.min_signal)/(reference-self.min_signal))
        
            if difference > 3*noise_level:
                detection = True
                detected = True
            if comparison > comp:
                detection = True
                detected = True
                
            
            Detected.append(detected)
                
        
        return detection, Detected
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        








