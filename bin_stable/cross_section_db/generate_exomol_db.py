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
create database for exomol cross section

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


import SEAS_Aux.cross_section.cross_section_database_generator as csdb
from SEAS_Utils.common_utils.DIRs import Exomol_Xsec,Exomol_DB
import SEAS_Utils.common_utils.configurable as config
import matplotlib.pyplot as plt


def exomol_db():
    """
    Simplistic database with one molecule
    with same resolution for all wavelength
    """    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    molecule        = "31P-1H3"
    outpath         = Exomol_DB
    db_name         = "Exomol.db"
    user            = "Azariven"
    
    db_info = [outpath, db_name, user]
    
    numin       = 400
    numax       = 9999
    step        = 1
    method      = "Doppler"
    source      = "Exomol"
    
    T_grid      = user_input["Simulation_Control"]["T_Grid"]
    P_grid      = user_input["Simulation_Control"]["P_Grid"]
    
    db_init = csdb.xsec_generator(db_info,numin,numax,step,method,source)
    
    d_path = os.path.join(Exomol_Xsec,molecule)
    print d_path
    
    db_init.exomol_inject(d_path,"PH3",P_grid,T_grid)


def exomol_db2():
    """
    Simplistic database with one molecule
    with same resolution for all wavelength
    """    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    molecule        = "31P-1H3"
    outpath         = Exomol_DB
    db_name         = "Exomol_PH3.db"
    user            = "Azariven"
    
    db_info = [outpath, db_name, user]
    
    numin       = 400
    numax       = 9999
    step        = 1
    method      = "Doppler"
    source      = "Exomol"
    
    T_grid      = user_input["Simulation_Control"]["T_Grid"]
    P_grid      = user_input["Simulation_Control"]["P_Grid"]
    
    db_init = csdb.xsec_generator(db_info,numin,numax,step,method,source)
    
    d_path = os.path.join(Exomol_Xsec,molecule)
    print d_path
    
    db_init.exomol_inject2(d_path,"PH3",P_grid,T_grid)
    
if __name__ == "__main__":


    exomol_db2()
    
    
    
    
    