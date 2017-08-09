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
from SEAS_Utils.common_utils.DIRs import Demo_DB
from SEAS_Utils.common_utils.DIRs import HITRAN_Lines
import SEAS_Utils.common_utils.configurable as config


def demo_db():
    
    data = config.Configuration("Demo_db.cfg")

    molecule_data   = data["Molecule"]
    outpath         = Demo_DB
    db_name         = "cross_sec_Demo.db"
    user            = "Azariven"
    
    db_info = [outpath, db_name, user]
    
    numin       = 400
    numax       = 30000
    step        = 10
    method      = "Voigt"
    source      = "HITRAN_LL"
    
    T_grid      = [200,300,400,500,600,700,800]
    P_grid      = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152]
    gamma       = "gamma_self"
    cross       = True
    
    db_init = csdb.xsec_generator(db_info,numin,numax,step,method,source)
    
    for molecule in molecule_data:
        component = np.array(molecule_data[molecule]["component"],dtype=np.int)
        d_path = os.path.join(HITRAN_Lines,molecule)
        db_init.single_molecule_inject(d_path,molecule,component,gamma,cross,P_grid,T_grid)
    
    
if __name__ == "__main__":

    demo_db()

    
    
    
    
    