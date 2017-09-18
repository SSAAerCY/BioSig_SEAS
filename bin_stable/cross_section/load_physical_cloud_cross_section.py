
import os
import sys
import numpy as np
import time
from scipy import interpolate
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_Main.atmosphere_effects.cloud import File_Cloud_Simulator




if __name__ == "__main__":
    
    legends = []
    for radius in [1,2,3]:
        index  = 1.33
        lambd  = np.arange(400,30000,10)
        lambd = 10000./lambd
        #radius = 3
        
        Simulator = File_Cloud_Simulator(lambd, radius)
        
        n = Simulator.calc_cloud_number_density(particle_radius=radius*10**-4)
        wav, xsec = Simulator.get_cloud_cross_section("CROSS_ZnS.npy")
        
        
        lag, = plt.plot(wav, xsec, label="R=%s um"%radius)
        legends.append(lag)
        
    ax = plt.gca()
    #ax.set_yscale('log')
    plt.legend(handles = legends)
    plt.title("ZnS Cloud Particle Cross Section with Varying Radius")
    plt.xlabel("Wavelength (um)")
    plt.ylabel("Intensity (cm^2/particle)")

    
    plt.show()
    