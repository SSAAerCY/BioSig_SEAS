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

Reflection Spectra Simulator

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
from SEAS_Main.atmosphere_effects.biosig_molecule import load_NIST_spectra, biosig_interpolate
from transmission_spectra_simulator import TS_Simulator

from SEAS_Aux.calculation.interpolation import interpolate1d
import SEAS_Aux.calculation.astrophysics as calc 
import SEAS_Aux.cross_section.hapi as hp

from SEAS_Utils.common_utils.constants import *
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.DIRs import TP_Profile_Data, Mixing_Ratio_Data, molecule_info, DB_DIR,Intermediate_DIR, HITRAN_CIA
from SEAS_Utils.common_utils.data_loader import two_column_file_loader,multi_column_file_loader, json_loader, molecule_cross_section_loader2
from SEAS_Utils.common_utils.data_saver import check_file_exist, check_path_exist
import SEAS_Utils.common_utils.db_management2 as dbm


class RS_Simulator(TS_Simulator):
    
    def __init__(self, user_input):
        super(RS_Simulator, self).__init__()
    
    
    def load_atmosphere_geometry_model(self):
        
        normalized_pressure         = self.normalized_pressure
        normalized_temperature      = self.normalized_temperature
        normalized_molecules        = self.normalized_molecules
        normalized_abundance        = self.normalized_abundance
        normalized_cross_section    = self.normalized_cross_section
        normalized_scale_height     = self.normalized_scale_height     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    