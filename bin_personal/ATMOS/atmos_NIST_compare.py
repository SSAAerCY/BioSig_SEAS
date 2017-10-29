"""

compare atmos result with nist result in TS simulation

"""

import os
import sys

import numpy as np
from scipy.special import wofz

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from matplotlib import ticker
ml = MultipleLocator(10)

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils.configurable as config
from SEAS_Aux.ATMOS_1.ATMOS_Xsec import ATMOS_1_Simulator
from SEAS_Main.atmosphere_effects.biosig_molecule import biosig_interpolate, load_NIST_spectra


def simple_compare():

    #smiles = 'NCC(O)(CC)'
    smiles = "SCCC(C)C"
    window = "earth2_atmosphere"
    
    pressure = 10e6
    temperature = 300

    Simulator = ATMOS_1_Simulator(smiles, window, pressure, temperature)
    nu  = np.arange(400,30000,1)
    xsec = Simulator.get_cross_section(nu)
    leg1, = plt.plot(10000./nu, xsec, label = "ATMOS")

    x1,y1 = load_NIST_spectra(smiles,["wn","A"],True)
    y2 = biosig_interpolate(x1,y1,nu,"A")
    leg2, = plt.plot(10000./nu, np.array(y2), label = "NIST")
    
    legends = [leg1, leg2]

    plt.legend(handles=legends)
    plt.xlabel("Wavenumber (cm^-1)")
    plt.ylabel("Absorption (proxy)")
    plt.title('Simulated ATMOS vs NIST Molecular Absorption for: %s'%smiles)
    plt.show()


if __name__ == "__main__":
    
    molecule_input = config.Configuration("atmos_data.cfg")

    window = "earth2_atmosphere"
    pressure = 10e6
    temperature = 300
    
    for smiles in molecule_input:
        print 
        Simulator = ATMOS_1_Simulator(smiles, window, pressure, temperature)
        nu  = np.arange(400,30000,1)
        xsec = Simulator.get_cross_section(nu)
        leg1, = plt.plot(10000./nu, xsec, label = "ATMOS")
    
        x1,y1 = load_NIST_spectra(smiles,["wn","A"],True)
        y2 = biosig_interpolate(x1,y1,nu,"A")
        leg2, = plt.plot(10000./nu, np.array(y2), label = "NIST")
        
        legends = [leg1, leg2]
    
        plt.legend(handles=legends)
        plt.xlabel("Wavenumber (cm^-1)")
        plt.ylabel("Absorption (proxy)")
        plt.title('Simulated ATMOS vs NIST Molecular Absorption for: %s'%smiles)
        plt.savefig("result/ATMOS_NIST_Compare_%s.png"%smiles)

