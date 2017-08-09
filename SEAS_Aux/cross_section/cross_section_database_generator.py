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
Cross Section Database Generator

uses the cross_section_calculator to calculate bulk cross sections and store in 
sqlite database

"""

import SEAS_Utils.common_utils.db_management2 as dbm
import SEAS_Aux.cross_section.cross_section_calculator as csc
from SEAS_Utils.common_utils.timer import simple_timer
import lines2xsec as l2x




class xsec_generator():
    
    def __init__(self,db_info,numin,numax,step,method,source):
        
        self.db_info    = db_info
        self.numin      = numin
        self.numax      = numax
        self.step       = step
        self.method     = method
        self.source     = source

        self.kwargs = {"dir"        :db_info[0],
                       "db_name"    :db_info[1],
                       "user"       :db_info[2],
                       "DEBUG"      :False,
                       "REMOVE"     :True,
                       "BACKUP"     :False,
                       "OVERWRITE"  :True}

        self.cross_db = dbm.database(**self.kwargs)
        
        if self.cross_db.is_db():
            self.access_db()
        else:
            self.create_db()
    
    def create_db(self):
        
        self.cross_db.create_db()
    
    def access_db(self):
        
        self.cross_db.access_db()
    
    def single_molecule_inject(self, d_path, molecule, component, gamma, cross, P_grid, T_grid):
        
        
        
        Timer = simple_timer()
        table_name = molecule
        columns = [("nu","real"),("coef","real"),("P","real"),("T","real"),("Resolution","real"),("Method","text"),("Source","text")]
        
        if self.cross_db.check_table_exist(table_name) == False:
            print "table created for {}".format(molecule)
            self.cross_db.create_table(table_name, *columns)
        else:
            self.cross_db.delete_table(table_name)
    
        # convert to csc in the future
        molecule_data = l2x.read_data(d_path,molecule,self.numin,self.numax,-1,-1)
        
    
        for P in P_grid:
            for T in T_grid:
                # huge room for multiprocessing here....
                print "Generating cross section for {},{} at P: {}pa and T: {}K".format(molecule,component,P,T),
                
                nu,coef = l2x.absorption_Voigt_calculation(molecule_data, component, gamma, P, T, self.numin, self.numax, self.step, cross)
                
                data = [x+(P,T,self.step,self.method,self.source) for x in zip(nu,coef)]        
                print "Generation_Time:", Timer.elapse(),
        
                for i in data:
                    self.cross_db._insert_data_single_dev(table_name, i)
                    
                self.cross_db.conn.commit()                        
                print "Injection_Time:", Timer.elapse()
        print "Total_Time:",Timer.total()
    
    def bulk_inject(self):
        
        pass



    def test_select(self):
        
        pass
    
    def dev_db(self):
        
        pass





