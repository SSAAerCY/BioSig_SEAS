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
Download the 2016 HITRAN Line List Database

"""


import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


import SEAS_Aux.data_downloader.web_downloader as wd



if __name__ == "__main__":
    
    data_output_path = "../../input/absorption_data/HITRAN_Line_List"

    ID = wd.get_HITRAN_ID()
    numin = 0
    numax = 50000
    
    for i in ID:
        wd.HITRAN_Line_List_downloader(data_output_path,i,numin,numax,True,False)


    
    
    
    