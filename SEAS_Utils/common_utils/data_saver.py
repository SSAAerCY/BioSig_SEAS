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
Misc. Data Saving functions

Code related to saving generated intermediate data
"""


import os,sys
import shutil
import numpy as np


def check_file_exist(file_path):

    return os.path.isfile(file_path)

def check_path_exist(save_path, create=True, overwrite=False, Verbose=False):
    """
    check if path exists. If path exist and overwrite is True, will remove the existing path and content
    then replace with an empty path. If path doesn't exist and create is True, will create an empty path
    
    """
    
    if os.path.isdir(save_path):
        
        if Verbose:
            print "Path exist, ",
        if overwrite:
            shutil.rmtree(save_path)
            os.makedirs(save_path)
            if Verbose:
                print "overwrite existing path with clean path"
        else:
            if Verbose:
                print "no overwrite"
        return True
    else:
        if Verbose:
            print "Path not exist,",
        if create:
            os.makedirs(save_path)
            if Verbose:
                print "path created"
        else:
            if Verbose:
                print "path not created"        
        return False

def save_txt(savepath, savename, data, extension=".txt", overwrite=True, check=False):

    if check:
        check_path_exist(savepath)

    save = os.path.join(savepath, savename)
    
    with open(save, "w") as f:
        f.write(data)
        f.close()
        
    
    
def save_npy(savepath, savename, data, overwrite=True, check=False):
    
    if check:
        check_path_exist(savepath)
        
    save = os.path.join(savepath, savename)

    np.save(save, data)
    
    
def save_excel(savepath, savename, data, overwrite=True):
    
    check_path_exist(savepath)
    save = os.path.join(savepath, savename)    



class Saver():
    
    
    def __init__(self, savepath, savename, data, overwrite=True):
        """
        overwrite should be on the file level
        """
        
        self.savepath = savepath
        self.savename = savename
        self.data = data
        self.overwrite = overwrite
        
        self.save = os.path.join(savepath, savename)
    
    
        self.check_data()
    
    def check_data(self):
        pass
    

    def check_file_exist(self):
    
        return os.path.isfile(self.save)
    
    def to_text(self):
        
        pass
    
    def to_npy(self):
        
        np.save(self.save, self.data)

    def to_excel(self):
        
        from openpyxl import Workbook, load_workbook
        
        pass




    