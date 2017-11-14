"""

calculate the number of photons per bin and associated photon noise


"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from scipy.integrate import simps
from numpy import trapz


DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.data_loader import multi_column_file_loader
from SEAS_Main.atmosphere_effects.biosig_molecule import biosig_interpolate


def calc_photon_energy(wavelength):
    
    h = 6.626e-34
    c = 2.99e8
    
    photon_energy = h*c/wavelength
    
    return photon_energy

def calc_photon_number(flux, photon_energy):
    
    photon_number = flux/photon_energy
    
    return photon_number
    
  
def calc_photon_noise_snr(photon_number):
    
    return photon_number/np.sqrt(photon_number)  


def calc_star_flux():
    
    pi = np.pi
    r = 0.1*6.95e8
    sigma = 5.67e-8
    T = 3000
    
    
    r = 6.95e8
    T = 5770
    
    
    L = 4*pi*r**2*sigma*T**4
    
    d = 1.5e11
    
    E = L/(4*pi*d**2)
    
    print L, E

def display_stellar_spectra():
    
    
    filename = "../../input/absorption_data/Stellar_Spectra/star_data/albedo_grid_K0V.dat"
    
    data = multi_column_file_loader(filename)
    x = np.array(data[0],dtype="float")
    y = np.array(data[1],dtype="float")
    
    nu  = np.arange(400,30000,1)
    
    wav = np.arange(0.3,4,0.01)
    normalized_flux = biosig_interpolate(x,y,wav,"C")
    
    
    area_under_curve = trapz(normalized_flux, dx=0.05)
    
    print area_under_curve
    
    
    #plt.plot(x,y)
    plt.plot(wav,normalized_flux)
    plt.show()



    
if __name__ == "__main__":
    
    
    wave = 800*1e-9
    flux = 1e-8
    
    energy = calc_photon_energy(wave)
    
    photon_number = calc_photon_number(flux,energy)
    
    print photon_number
    
    photon_noise_snr = calc_photon_noise_snr(photon_number)
    
    print photon_noise_snr
    
    calc_star_flux()
    display_stellar_spectra()
    
    
    
    
    