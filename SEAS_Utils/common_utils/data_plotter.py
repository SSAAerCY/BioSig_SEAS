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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
ml = MultipleLocator(10)

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils as utils
from SEAS_Utils.common_utils.data_saver import check_path_exist, check_file_exist


class Plotter():
    """
    A more genralized plotter for multiple purpose needs
    """
    
    def __init__(self):
        pass


class Simulation_Plotter():
    """
    Plotter for the transit simulation that requires an user input cfg
    """
    '''
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
    '''
    
    def __init__(self, user_input):
        
        self.user_input = user_input
        
        fig_w = utils.to_int(self.user_input["Plotting"]["Figure"]["figsize_w"])
        fig_h = utils.to_int(self.user_input["Plotting"]["Figure"]["figsize_h"])
        
        self.fig = plt.figure(figsize=(fig_w, fig_h))
        self.ax = plt.gca()
        
        plt.tick_params(axis='x', which='minor')
        self.ax.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))   
        
        self.x_unit  = self.user_input["Plotting"]["Figure"]["x_unit"]
        self.x_scale = self.user_input["Plotting"]["Figure"]["x_scale"]
        self.x_label = self.user_input["Plotting"]["Figure"]["x_label"]
        self.x_multi = utils.to_float(self.user_input["Plotting"]["Figure"]["x_multiplier"])
        
        self.y_unit  = self.user_input["Plotting"]["Figure"]["y_unit"]
        self.y_scale = self.user_input["Plotting"]["Figure"]["y_scale"]
        self.y_label = self.user_input["Plotting"]["Figure"]["y_label"]
        self.y_multi = utils.to_float(self.user_input["Plotting"]["Figure"]["y_multiplier"])
        
        self.Title   = self.user_input["Plotting"]["Figure"]["Title"]
        
        self.set_scale()
        
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
        plt.title(self.Title)
        
    def plot_xy(self, x, y, Label="Data", Color = ""):
        
        if self.x_unit == "um":
            x = 10000./np.array(x)
        elif self.x_unit == "nm":
            x = 10000000./np.array(x)
        else:
            x = np.array(x)
        
        y = np.array(y)
    
        if Color == "":
            plt_ref, = plt.plot(x*self.x_multi,y*self.y_multi,label=Label)            
        else:
            plt_ref, = plt.plot(x*self.x_multi,y*self.y_multi,label=Label,color=Color)

        return plt_ref
    
    def plot_xy_list(self, xlist, ylist, ):
        pass
    
    def plot_window(self, window, Color="k", Alpha=0.2):
        
        for k in window:
            if self.x_unit == "um":
                up,down = 10000./k[1],10000./k[0]
            elif self.x_unit == "nm":
                up,down = 10000000./k[1],10000000./k[0]          
            else:
                up,down = k[0],k[1]
            
            plt.axvspan(up,down,facecolor=Color,alpha=Alpha)
        
    def plot_line(self):
        pass
    
    def set_scale(self):
        
        if self.x_scale == "log":
            self.ax.set_xscale('log')
        if self.y_scale == "log":
            self.ax.set_yscale('log')
    
    def set_legend(self, legends):
        plt.legend(handles=legends)

    def save_plot(self):
        
        save_dir  = self.user_input["Save"]["Plot"]["path"]
        save_name = self.user_input["Save"]["Plot"]["name"] 
        plt.savefig(os.path.join(save_dir,save_name))

    def show_plot(self):
        
        plt.show()  
    
