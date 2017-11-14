"""
comparing between hitran and nist absorption cross sections

"""


import os
import sys
import numpy as np
import random

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.data_loader import NIST_Smile_List, HITRAN_to_NIST


def test_HITRAN_to_NIST():
    
    
    molecule = "H2O"
    molecule = "HCOOH"
    
    molecule = HITRAN_to_NIST(molecule,"Smiles")
    
    print molecule


if __name__ == "__main__":
    
    test_HITRAN_to_NIST()









