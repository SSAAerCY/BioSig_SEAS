"""



"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyfits

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.DIRs import Stellar_Spectra



def display_stellar_spectra():

    folder = "star_data"
    name = "albedo_grid_G5V.dat"
    
    filename = os.path.join(Stellar_Spectra,folder,name)

    f = open(filename)
    data = []
    for i in f:
        data.append(i.split())
    
    sdata = np.array(data).T
    
    front = 0
    back = 10000
    
    x = sdata[0][front:back]
    y = sdata[1][front:back]
    
    plt.plot(x,y)
    plt.show()
        
def varying_stellar_spectra():
    
    orbit = 5.  #AU
    
    R_E = 6400000
    R_S = 6.95*10**8
    R_P = 1*R_E
    P_Surface = np.pi*R_P**2
    S_Surface = np.pi*R_S**2
    
    
    filename = os.path.join(Stellar_Spectra,"star_data","albedo_grid_G5V.dat")
    #filename = os.path.join(Stellar_Spectra,"star_data","SOLARSPECTRUM.dat")
    f = open(filename)
    data = []
    for i in f:
        data.append(i.split())
    
    sdata = np.array(data).T
    
    x = np.array(sdata[0][:20000],dtype="float")
    
    # y is flux density at 1AU
    y = np.array(sdata[1][:20000],dtype="float")
    
    
    #assuming that 1% of light received by the planet reflected to observer
    modifier = (1./orbit)**2*P_Surface*(0.01)
    y = y*modifier
    
    
    
    plt.plot(x,y)
    plt.show()


if __name__ == "__main__":
    
    #display_stellar_spectra()
    varying_stellar_spectra()
    
    
    
    
    
    