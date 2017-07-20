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
Misc. Data Loading functions


Should I put the hitran data reader here as well?
i mean, after all the user should want control what data line list data
it's reading in from...

"""

import os
import sys
import json
import numpy as np
import SEAS_Utils.common_utils.db_management2 as dbm


def two_column_file_loader(path,spliter="",type="float"):
    """
    load data from files that contains only two column
    """
    
    with open(path) as data_file:   
        if spliter == "":
            xdata,ydata = (np.array([x.split() for x in data_file.read().split("\n")])).T
        else:
            xdata,ydata = (np.array([x.split(spliter) for x in data_file.read().split("\n")])).T
            
        if type == "float":
            return np.array(xdata,dtype=np.float),np.array(ydata,dtype=np.float)



def multi_column_file_loader(path,spliter="",type="float"):
    """
    load data from files that contains multiple columns
    """
    
    with open(path) as data_file:   
        if spliter == "":
            data = np.array([x.split() for x in data_file.read().split("\n")]).T
        else:
            data = np.array([x.split(spliter) for x in data_file.read().split("\n")]).T
        

        #if data[-1] == []:
        #    data = data[:-1]
        
        if type == "float":
            return np.array(data,dtype=np.float)
        elif type == "mixed":
            return data
    

def json_loader(path):
    
    with open(path) as data_file:
        data = json.load(data_file)

    return data



def database_loader(input):
    "TBD,"
    pass


def cross_section_loader(inputs, cross_db, molecule):

    Pgrid = inputs["Simulation_Control"]["P_Grid"]
    Tgrid = inputs["Simulation_Control"]["T_Grid"]
    numin = inputs["Spectra"]["Numin"]
    numax = inputs["Spectra"]["Numax"]


    data = [[0 for y in range(len(Pgrid))] for x in range(len(Tgrid))]
    for i,T in enumerate(Tgrid):
        for j,P in enumerate(Pgrid):
            
            
            result = cross_db.c.execute("SELECT nu, coef FROM {} WHERE P={} AND T={} AND nu>={} AND nu<{} ORDER BY nu".format(molecule,P,T,numin,numax))
            fetch = np.array(result.fetchall()).T
            data[i][j] = fetch[1]
            
    nu = fetch[0]

    
    return nu, data

def cross_section_loader2(inputs, outpath, molecule):

    T_grid = inputs["Simulation_Control"]["T_Grid"]
    P_grid = inputs["Simulation_Control"]["P_Grid"]
    numin = inputs["Spectra"]["Numin"]
    numax = inputs["Spectra"]["Numax"]

    kwargs = {"dir"        :outpath,
              "db_name"    :"cross_section_Simulation_%s.db"%molecule,
              "user"       :"azariven",
              "DEBUG"      :False,
              "REMOVE"     :True,
              "BACKUP"     :False,
              "OVERWRITE"  :True}
    
    cross_db = dbm.database(**kwargs)  
    cross_db.access_db()  

    data = [[0 for y in range(len(T_grid))] for x in range(len(P_grid))]
    for j,P in enumerate(P_grid):
        for i,T in enumerate(T_grid):
            table_name = "T%sP%s"%(i,j)
            result = cross_db.c.execute("SELECT * FROM {} ORDER BY nu".format(table_name))
            fetch = np.array(result.fetchall()).T
            data[j][i] = fetch[1]
            
    nu = fetch[0]

    
    return nu, data




def HITRAN_Line_List_reference():
    
    from SEAS_Utils.common_utils.DIRs import HITRAN_Molecule_List
    
    molecule = open(HITRAN_Molecule_List).read().split("\n")
    
    component = []
    for i in range(len(molecule)):
        component.append([i+1,1,1])
        
    return molecule,component





