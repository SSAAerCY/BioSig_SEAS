"""

This code is meant to test simulation of ATMOS in TS simulation


"""

import os
import sys

import matplotlib.pyplot as plt


DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Aux.ATMOS_1.ATMOS_Xsec import ATMOS_1_Simulator





if __name__ == "__main__":
    
    plotted_molecule = 'NCC(O)(CC)'
    Simulator = ATMOS_1_Simulator(plotted_molecule)
    
    functional_dictionary = Simulator.load_functional_dictionary()
    
    plotables = Simulator.load_plotables()
    
    molecules = Simulator.load_molecules(functional_dictionary)

    # Finds all the strong features whose average frequency is within a window
    atmosphere = Simulator.load_atmosphere_windows("earth2_atmosphere")
    
    Simulator.calculate_detection(atmosphere, molecules, plotables)
    
    example = molecules[plotted_molecule]
    
    xs, ys = zip(*example.average_points())
    
    filtered_list = list(filter(lambda x: 1600 < x[0] < 1900, example.average_points())) #window = range(3250, 3450)
    
    markerline, stemlines, baseline = plt.stem(xs, ys, '-')
    plt.setp(baseline, 'color', 'r', 'linewidth', 1)
    plt.xlabel("Wavenumber (cm^-1)")
    plt.ylabel("Intensity (proxy)")
    plt.title('Molecule Plotted : %s'%plotted_molecule)
    plt.show()

