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
Select which CIA will be chosen based on a given list of molecules
"""

import os
import sys

from SEAS_Utils.common_utils.DIRs import HITRAN_CIA



class cia_selector():
    
    def __init__(self,molecules):
        
        self.molecules = molecules
    
    def select_fuzzy(self):
        """
        the fuzzy search will search the data for molecules match position 1 and 2 
        """
        
        filename = []
        for file in os.listdir(HITRAN_CIA):
            info = file.replace("_","-").split("-")
            if "cia" in file and len(info) == 3:
                if info[0] in self.molecules and info[1] in self.molecules:
                    filename.append(os.path.join(HITRAN_CIA,file))
        
        return filename
    
    def select_exact(self):
        """
        the exact search will search the data for molecules and positions
        """
        filename = []
        for molecule in self.molecules:
            a = molecule[0]
            b = molecule[1]
            for file in os.listdir(HITRAN_CIA):
                info = file.replace("_","-").split("-")
                if "cia" in file and len(info) == 3:
                    if info[0] in self.molecules and info[1] in self.molecules:
                        filename.append(os.path.join(HITRAN_CIA,file))            
        
        
    