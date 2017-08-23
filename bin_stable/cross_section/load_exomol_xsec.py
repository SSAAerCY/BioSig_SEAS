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

Load cross sections form exomol datafiles

"""

import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.configurable as config
import SEAS_Utils.common_utils.db_management2 as dbm
from SEAS_Utils.common_utils.DIRs import Exomol_Xsec, Exomol_DB
from SEAS_Utils.common_utils.data_loader import two_column_file_loader,exomol_cross_section_loader
from SEAS_Utils.common_utils.data_plotter import Plotter

def simple_load():
    
    ph3_data = Exomol_Xsec+"/31P-1H3/31P-1H3_400-9999_300K_1.000000.sigma"
    nu, xsec = two_column_file_loader(ph3_data)
    
    P = 0.0306
    T = 275
    pathl = 1
    n = 10**25
    
    tau = n*xsec*pathl
    trans = np.e**(-tau)
    
    plotter = Plotter() 
    plotter.plot_xy(10000./nu,trans)

def db_load():


    inputs = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")


    molecule        = "PH3"
    
    kwargs = {"dir"        :Exomol_DB,
              "db_name"    :"Exomol_PH3.db",
              "user"       :"Azariven",
              "DEBUG"      :False,
              "REMOVE"     :True,
              "BACKUP"     :False,
              "OVERWRITE"  :True}

    cross_db = dbm.database(**kwargs)
    cross_db.access_db()

    result = cross_db.c.execute("SELECT nu, coef FROM {} WHERE P={} AND T={} AND nu>={} AND nu<={} ORDER BY nu".format(molecule,13500.0,300,400,9999))
    fetch = np.array(result.fetchall()).T 
    
    nu = fetch[0]
    coef = fetch[1]
    plotter = Plotter() 
    plotter.plot_xy(10000./nu,coef)

def db_load_bulk():

    inputs = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    molecule        = "PH3"
    
    kwargs = {"dir"        :Exomol_DB,
              "db_name"    :"Exomol_PH3.db",
              "user"       :"Azariven",
              "DEBUG"      :False,
              "REMOVE"     :True,
              "BACKUP"     :False,
              "OVERWRITE"  :True}

    cross_db = dbm.database(**kwargs)
    cross_db.access_db()
    
    nu, data = exomol_cross_section_loader(inputs, cross_db, molecule)
    
    print len(nu)
    

if __name__ == "__main__":
    
    db_load_bulk()
    
    
    
    
    
    