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

"""


import os
import sys
from Tkinter import *
DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_Aux.Visualizer.SEAS_Visualizer import SEAS_Main_GUI


def test_visualizer():
    
    
    root = Tk()
    root.wm_title("Spectra Search for All Small Molecules V")
    GUI = SEAS_Main_GUI(root)
    root.mainloop() 


if __name__ == "__main__":

    test_visualizer()








