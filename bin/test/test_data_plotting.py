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
Test if data plotter is functioning correctly
 

"""


import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.data_plotter as plotter
import SEAS_Utils.common_utils.configurable as config


if __name__ == "__main__":
    
    
    x = [1,2,3]
    y = [2,4,5]
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    
    user_input["Plotting"]["Figure"]["Title"] = "Testing Plot"


    plotter = plotter.Simulation_Plotter(user_input)
    
    
    plt_ref = plotter.plot_xy(x,y)
    
    plotter.set_legend([plt_ref])
    
    
    plotter.show_plot()







    
    
    