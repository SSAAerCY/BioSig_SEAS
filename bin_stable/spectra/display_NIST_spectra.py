"""

Code attempting to figure out the state of each spectra from NIST downloaded

  
"""


import os
import sys
from openpyxl import Workbook, load_workbook
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils as utils
from SEAS_Utils.common_utils.DIRs import molecule_info, NIST_Spectra

import SEAS_Utils.common_utils.jdx_Reader as jdx


def display_NIST_spectra(CAS,file):
    
    name = os.path.join(NIST_Spectra,CAS,"%s_%s.jdx"%(CAS,file))

    data = jdx.JdxFile(name)
    
    x,y = data.wl(), data.trans()
    
    plt.plot(x,y)
    plt.show()




if __name__ == "__main__":
    
    CAS = "C6320963"
    file = 0
    
    display_NIST_spectra(CAS,file)








