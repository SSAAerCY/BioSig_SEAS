import os
import sys
import numpy as np
import random

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.data_loader import NIST_Smile_List


def test_NIST():
    
    
    molecule_type = "NIST"
    info = NIST_Smile_List()
    molecule_smiles = info[0]
    formula = info[3]
    IUPAC_name = info[4]
    
    for i,mo in enumerate(molecule_smiles):
        if "CSC" in mo:
            print formula[i]
        


if __name__ == "__main__":
    
    test_NIST() 
    