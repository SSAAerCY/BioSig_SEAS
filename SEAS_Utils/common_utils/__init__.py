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
This code is inherited the util/__init__.py in TSIG code written by
Martin Owens, et al. 
"""

from .constants import *
import os,sys
import shutil

def to_bool(x):
    try:
        if x.lower() in ['true', 'yes']:
            return True
        elif x.lower() in ['false', 'no']:
            return False
    except AttributeError:
        pass
        try:
            return bool(int(x))
        except (ValueError, TypeError):
            pass
    raise ValueError("Unknown boolean specifier: %s" % x)

def to_float(x):
    if x is None:
        return None
    return float(x)

def to_int(x):
    if x is None:
        return None
    return int(x)

def isfile(path):

    return os.path.isfile(path)

def isdir(path):
    
    return os.path.isdir(path)


def check_path(save_path, overwrite=False):
    """
    check if path exists. If path exist and overwrite is True, will remove the existing path and content
    then replace with an empty path. If path doesn't exist and create is True, will create an empty path
    
    """
    
    if overwrite:
        if os.path.isdir(save_path):
            shutil.rmtree(save_path)
            os.makedirs(save_path)
            print "Overwrite Existing Path with clean Path"

        else:
            os.makedirs(save_path)
            print "Created Path for Simulation Saves"
            
    else:
        if os.path.isdir(save_path):
            print "Path Exist but No Overwrite Permission"

        else:
            os.makedirs(save_path)
            print "Created Path for Simulation Saves"        
        

        

