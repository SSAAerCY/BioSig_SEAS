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

"""
import os
import sys
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)

import SEAS_Main.atmosphere_geometry
import SEAS_Main.atmosphere_property
import SEAS_Main.effects

from SEAS_Aux.calculation.interpolation import interpolate1d
import SEAS_Aux.calculation.astrophysics as calc 

from SEAS_Utils.common_utils.constants import *
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.DIRs import TP_Profile_Data, Mixing_Ratio_Data, molecule_info, DB_DIR,Intermediate_DIR
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
    
    
    def example_simulate(self):
        
        
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
        print "load time", self.Timer.elapse()

        self.effects = self.load_effects()


        # determine atmosphereic effects to add in
        self.Final_TS = self.load_atmosphere_geometry_model()
        print "calc time", self.Timer.elapse()
        
        
        
        plt.title("Transit Signal for Simulated Atmosphere")
        plt.xlabel(r'Wavelength ($\mu m$)')
        plt.ylabel("Transit Signal (ppm)")    
        
        
        plt.plot(10000./self.nu, self.Final_TS)
        
        ax = plt.gca()
        ax.set_xscale('log')
        plt.tick_params(axis='x', which='minor')
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))          
        
        plt.show()
        return self.Final_TS

    def CIA_simulate(self):
        pass

        

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
        pass
    
    def interpolate_rayleigh_scattering(self):
        pass
    
    def load_CIA(self):
        pass
    
    def interpolate_CIA(self):
        pass




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

    def load_atmosphere_geometry_model(self):
        
        
        
        TotalBeams = len(self.normalized_pressure)
        
        normalized_pressure         = self.normalized_pressure
        normalized_temperature      = self.normalized_temperature
        normalized_molecules        = self.normalized_molecules
        normalized_abundance        = self.normalized_abundance
        normalized_cross_section    = self.normalized_cross_section
        normalized_scale_height     = self.normalized_scale_height     
        
        
        Total_Transit_Signal = np.ones(len(self.nu))*self.Base_TS_Value
        base_layer = self.R_planet


        Rayleigh_array = []
        for molecule in normalized_molecules:
            Rayleigh_array.append(calc.calc_rayleigh(molecule, self.nu))


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
                for m,molecule in enumerate(normalized_molecules):        
                    
                    #weird how abundance and cross section are wired differently
                    molecular_ratio = normalized_abundance[j+i][m]
                    number_density = (normalized_pressure[j+i]/(BoltK*normalized_temperature[j+i]))*molecular_ratio
                    
                    rayleigh = Rayleigh_array[m]
                    sigma = normalized_cross_section[m][j+i][j+i]
                    
                    #rayleigh_scat, etc should be pre calculated
                    effects = sigma+rayleigh#+CIA+cloud
        
                    ChunkTau_Per_Molecule = number_density*(effects)*pathl*2*0.0001
        
                    if ChunkTau == []:
                        ChunkTau = ChunkTau_Per_Molecule
                    else:
                        ChunkTau += ChunkTau_Per_Molecule   
                        
                if BeamTau == []:
                    BeamTau = ChunkTau
                else:
                    BeamTau += ChunkTau                       
            
            
            BeamTrans = calc.calc_transmittance(BeamTau)  
            RingArea = (base_layer**2-prev_layer**2)/self.R_Star**2
            
            Ring_Transit_Signal = (1-BeamTrans)*RingArea
            Total_Transit_Signal += Ring_Transit_Signal
            
            # update to the next beam up
            prev_layer = base_layer
            base_layer += normalized_scale_height[i]        
            
            
        Raw_Transit_Signal = Total_Transit_Signal
        return Raw_Transit_Signal

    def load_effects(self):
        """
        setup the effect template for what goes into each "test tube"
        """
        pass    
    
    def load_other_effects(self):
        pass

if __name__ == "__main__":
    
    import SEAS_Utils.common_utils.configurable as config
    user_input = config.Configuration("../../bin/dev/user_input_dev.cfg")
    simulation = TS_Simulator(user_input)
    Raw_TS = simulation.simulate()








