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
Data Plotting Tools

More Modes is going to come

"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_Utils.common_utils.data_saver import check_path_exist, check_file_exist

class Plotter():
    
    def __init__(self, x, y,
                 line_style = "-", color = None,
                 title = "A plot", xlabel = "xlabel", ylabel = "ylabel",
                 show = True, time = None,
                 save = False, save_path = "", save_name = "Temp_Plot.png",
                 overwrite = False
                 ):
        """
        x and y are double arrays here... or in another plotter?
        or ... simply make this powerful?
        """
        
        self.Figure = plt.Figure()
        
        self.x = y
        self.y = y
        
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title = title
        
        self.line_style = line_style
        self.color = color
        
        

        self.show = show
        self.time = time
        
        
        self.save = save
        
        
        check_path_exist(save_path)
        
        self.savename = os.path.join(save_path,save_name)        
    
    

    
    
    def set_log(self):
        pass
    

    def show_plot(self):
        
        plt.show()  
    
