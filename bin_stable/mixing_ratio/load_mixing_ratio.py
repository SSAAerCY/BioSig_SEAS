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
This is an example simulation for beginner users to try out the code
"""

import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_Utils.common_utils.data_loader import multi_column_file_loader
from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data


if __name__ == "__main__":
    filename = os.path.join(Mixing_Ratio_Data,"test_earth.txt")
    data = multi_column_file_loader(filename,type="mixed"))
    print data









