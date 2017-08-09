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
create an example database using small portion of the database

we should be creating database with an info table that help describe what the 
database contains

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


import SEAS_Aux.cross_section.cross_section_database_generator as csdb
from SEAS_Utils.common_utils.DIRs import Example_DB
from SEAS_Utils.common_utils.DIRs import HITRAN_Lines
import SEAS_Utils.common_utils.configurable as config
import matplotlib.pyplot as plt


def example_db():
    """
    Simplistic database with one molecule
    with same resolution for all wavelength
    """    
    data = config.Configuration("Example_db.cfg")

    molecule_data   = data["Molecule"]
    outpath         = Example_DB
    db_name         = "cross_sec_Example.db"
    user            = "Azariven"
    
    db_info = [outpath, db_name, user]
    
    numin       = 3000
    numax       = 4000
    step        = 10
    method      = "Voigt"
    source      = "HITRAN_LL"
    
    T_grid      = [250,300,350]
    P_grid      = [100000,10000,1000]
    gamma       = "gamma_self"
    cross       = True
    
    db_init = csdb.xsec_generator(db_info,numin,numax,step,method,source)
    
    for molecule in molecule_data:
        component = np.array(molecule_data[molecule]["component"],dtype=np.int)
        d_path = os.path.join(HITRAN_Lines,molecule)
        db_init.single_molecule_inject(d_path,molecule,component,gamma,cross,P_grid,T_grid)
    
    
if __name__ == "__main__":


    example_db()
    
    
    
    
    