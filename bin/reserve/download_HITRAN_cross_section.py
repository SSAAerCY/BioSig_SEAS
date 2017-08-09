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
Download the HITRAN Absorption Cross Section Database

This code is a bit tricky and takes some time to make it working
basically hitran website is now java scripted, so data scraping is a bit challenging

However, the data file is still downloadable, but you need to know the correct prior informations
which should be stored in a table (manually created)

Take Ozone for example,
http://hitran.org/data/xsec/O3_220.0_0.0_29164.0-40798.0_04.xsc

the link is divided into

<molecule>_<temperature>_<pressure>_<numin>-<numax>_<resolution>.xsc


"""


import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


import SEAS_Aux.data_downloader.web_downloader as wd



if __name__ == "__main__":
    pass