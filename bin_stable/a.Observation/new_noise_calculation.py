"""

This module is to calculate the new noise model


"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.constants import *
import SEAS_Utils.common_utils.configurable as config

def blackbody_lam(wav, T):
    """ Blackbody as a function of wavelength (m) and temperature (K).
    """
    a = 2*HPlanck*CLight**2
    b = HPlanck*CLight/(wav*BoltK*T)
    intensity = a/((wav**5)*(np.exp(b)-1.0))
    return intensity


def calculate_noise(user_input,bin_centers,bin_width):
    """
    bin_centers and bin_width are taken with unit of microns
    
    """

    R_Star = user_input["Star"]["R_Star"]
    R_planet = user_input["Planet"]["R_Planet"] 
    R_obs = user_input["Telescope"]["Aperture"]
    R_atmosphere = user_input["Planet"]["R_Atmosphere"]
    Distance = user_input["Telescope"]["Distance"]
    Duration = user_input["Telescope"]["Duration"]
    Quantum  = user_input["Telescope"]["Quantum_Efficiency"]
    T_Star = user_input["Star"]["T"]
    Noise_M = user_input["Observation_Effects"]["Noise"]["Multiplier"]
    
    # calculate number of photons 
    B_Body = blackbody_lam(bin_centers*10**-6, T_Star) #J/s/m^2/m/Steradian
    Bin_width = bin_width*10**-6 
    A_Star = np.pi*R_Star**2
    Psi_Tele = np.pi*R_obs**2/Distance**2
    E_Total = B_Body*Bin_width*A_Star*Psi_Tele*Duration
    num_photon = (E_Total*bin_centers*10**-6)/(HPlanck*c)*Quantum

    # calculate photon noise
    signal = 2.*R_planet*R_atmosphere/R_Star**2.
    photon_noise = Noise_M/np.sqrt(num_photon)
    snr = signal/photon_noise

    return signal,photon_noise,snr


def determine_bin(user_input):

    lambda_init = user_input["Telescope"]["min_wavelength"]
    lambda_max = user_input["Telescope"]["max_wavelength"]
    bin_width_init = user_input["Observation_Effects"]["bin_width"]
    bin_exponent = user_input["Observation_Effects"]["bin_exponent"]
    
    bin_width,bin_centers,bin_edges = [],[],[]
    
    i=0
    lambda_current = lambda_init
    while True:
        
        new_bin_width = bin_width_init*(lambda_current/lambda_init)**(bin_exponent)
        #new_bin_width = bin_width_init
        lambda_center = lambda_current+0.5*new_bin_width
        lambda_current += new_bin_width
        
        if lambda_center > lambda_max:
            print i
            break
        
        bin_edges.append(lambda_current)
        bin_centers.append(lambda_center)
        bin_width.append(new_bin_width)
    
        
        i+=1
        
    return np.array(bin_edges), np.array(bin_width), np.array(bin_centers)
        
        
        
if __name__ == "__main__":
    
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")

    user_input["Star"]["R_Star"] = 0.6*R_Sun
    user_input["Star"]["T"] = 4000.
    user_input["Planet"]["R_Planet"] = 1*R_Earth
    user_input["Planet"]["R_Atmosphere"] = 40*1000
    user_input["Telescope"]["Aperture"] = 6.5
    user_input["Telescope"]["Quantum_Efficiency"] = 1
    user_input["Telescope"]["Distance"] = 3*Psec
    user_input["Telescope"]["Duration"] = 100*3600
    user_input["Telescope"]["min_wavelength"] = 1
    user_input["Telescope"]["max_wavelength"] = 25
    user_input["Observation_Effects"]["Noise"]["Multiplier"] = 1.2
    user_input["Observation_Effects"]["bin_exponent"] = 3/2.
    user_input["Observation_Effects"]["bin_width"] = 0.01
    
    
    bin_edges, bin_width, bin_centers = determine_bin(user_input)
    signal,photon_noise,snr = calculate_noise(user_input,bin_centers,bin_width)
    
    print snr
    """
    Mdwarf has more FarUV (magnetic field), less MUV (cooler photosphere)
    Mdwarf has 100 - 1000x less bioactive UV => minimum threshold needed
    But... Mdwarf flares! 10x more uv than the sun. (1 every 10-40 days) ~ 20% mdwarfs
    """



