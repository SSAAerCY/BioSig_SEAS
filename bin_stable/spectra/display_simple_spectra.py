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
This is an example of showing spectra from a molecule

Temperature Pressure Grid to sample from 
    T_grid          = [100,150,200,250,275,300,325,350,400]
    P_grid          = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 
                       0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05]
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.db_management2 as dbm
from SEAS_Utils.common_utils.DIRs import Simulation_DB
#import SEAS_Main.simulation.astrophysics as astro


if __name__ == "__main__":
    
    kwargs = {"dir"        :"/Users/mac/Workspace/BioSig2/SEAS_Utils/common_utils/../../input/database/Simulation_Band",
              "db_name"    :"cross_section_Simulation.db",
              "user"       :"azariven",
              "DEBUG"      :False,
              "REMOVE"     :True,
              "BACKUP"     :False,
              "OVERWRITE"  :True}

    cross_db = dbm.database(**kwargs)
    cross_db.access_db()
    
    molecule = "CS"
    P = 0.0306
    T = 275
    numin = 400
    numax = 30000
    
    result = cross_db.c.execute("SELECT nu, coef FROM {} WHERE P={} AND T={} AND nu>={} AND nu<{} ORDER BY nu".format(molecule,P,T,numin,numax))
    nu,xsec = np.array(result.fetchall()).T
    
    
    pathl = 1
    n = 10**25
    
    tau = n*xsec*pathl
    
    trans = np.e**(-tau)
    plt.plot(nu,trans)
    plt.show()
    

    
    
    
    
    



