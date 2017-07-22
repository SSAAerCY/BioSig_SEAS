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

This is an demonstration code for loading cross sections from database


not urgently need to be implemented since need fixing with molecule db imports. placeholder for now

"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.cross_section.cross_section_calculator as csc
import SEAS_Utils.common_utils.data_plotter as plt
from SEAS_Utils.common_utils.DIRs import Temp_DIR, HITRAN_Water_Lines
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.data_loader import molecule_cross_section_loader






