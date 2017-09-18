"""

This code is meant to test simulation of ATMOS in TS simulation


"""

import os
import sys

import numpy as np
from scipy.special import wofz

import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Aux.ATMOS_1.ATMOS_Xsec import ATMOS_1_Simulator
from SEAS_Main.atmosphere_effects.biosig_molecule import biosig_interpolate, load_NIST_spectra



if __name__ == "__main__":
    
    #smiles = 'NCC(O)(CC)'
    smiles = "C(C)NCC(O)"
    window = "earth2_atmosphere"
    
    pressure = 10e6
    temperature = 300

    Simulator = ATMOS_1_Simulator(smiles, window, pressure, temperature)
    
    nu  = np.arange(400,10000,1)
    xsec = Simulator.get_cross_section(nu)

    plt.plot(10000./nu, xsec)
    plt.xlabel("Wavenumber (cm^-1)")
    plt.ylabel("Intensity (proxy)")
    plt.title('Simulated ATMOS Molecule Cross Section for: %s'%smiles)
    plt.show()

