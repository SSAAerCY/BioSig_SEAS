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


def check_file_exist(file_path):

    return os.path.isfile(file_path)

def check_path_exist(save_path, create=True, overwrite=False):
    """
    check if path exists. If path exist and overwrite is True, will remove the existing path and content
    then replace with an empty path. If path doesn't exist and create is True, will create an empty path
    
    """
    
    if os.path.isdir(save_path):
        print "Path exist, ",
        if overwrite:
            shutil.rmtree(save_path)
            os.makedirs(save_path)
            print "overwrite existing path with clean path"
        else:
            print "no overwrite"
        return True
    else:
        print "Path not exist, "
        if create:
            os.makedirs(save_path)
            print "path created"
        else:
            print "path not created"        
        return False
    
    
    
    