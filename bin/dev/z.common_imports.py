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

Here is a list of useful common imports for SEAS
This code should not be executed

"""
import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

#plotting
import SEAS_Utils.common_utils.data_plotter as plt
import matplotlib.pyplot as plt

#timer
from SEAS_Utils.common_utils.timer import simple_timer

#dbm
import SEAS_Utils.common_utils.db_management2 as dbm

#config
import SEAS_Utils.common_utils.configurable as config

#DIR
from SEAS_Utils.common_utils.DIRs import Simulation_DB

#constants
from SEAS_Utils.common_utils.constants import *

# data
from SEAS_Utils.common_utils.data_loader import two_column_file_loader
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
import SEAS_Utils.common_utils.jdx_Reader as jdx

# save data
from SEAS_Utils.common_utils.data_saver import check_path_exist, check_file_exist



