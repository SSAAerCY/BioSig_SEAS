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


def calculate_noise(user_input,lam=10**-6,bin=0.1*10**-6):

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
    B_Body = blackbody_lam(lam, T_Star) #J/s/m^2/m/Steradian
    Bin_width = bin 
    A_Star = np.pi*R_Star**2
    Psi_Tele = np.pi*R_obs**2/Distance**2
    E_Total = B_Body*Bin_width*A_Star*Psi_Tele*Duration
    num_photon = (E_Total*lam)/(HPlanck*c)*Quantum
    print num_photon,
    # calculate photon noise
    photon_noise = Noise_M/np.sqrt(num_photon)
    print photon_noise, 
    snr = (2.*R_planet*R_atmosphere/R_Star**2.)/(photon_noise)#R_Star**2/(2*R_planet*R_atmosphere)*photon_noise#*10**6
    print R_Star**2/(2*R_planet*R_atmosphere),transmission_noise,
    sys.exit()
    return transmission_noise


def determine_bin(user_input):

    lambda_init = user_input["Telescope"]["min_wavelength"]
    lambda_max = user_input["Telescope"]["max_wavelength"]
    bin_width_init = user_input["Observation_Effects"]["bin_width"]
    bin_exponent = user_input["Observation_Effects"]["bin_exponent"]
    
    bin_width,bin_centers,error_bar = [],[],[]
    
    i=0
    lambda_current = lambda_init
    while True:
        
        #new_bin_width = bin_width_init*(lambda_current/lambda_init)**(bin_exponent)
        new_bin_width = bin_width_init
        lambda_center = lambda_current+0.5*new_bin_width
        lambda_current += new_bin_width
        
        if lambda_center > lambda_max:
            print i
            break
        
        bin_centers.append(lambda_center)
        bin_width.append(new_bin_width)
        noise = calculate_noise(user_input,lambda_center*10**-6,new_bin_width*10**-6)
        error_bar.append(noise)
        
        print lambda_center, new_bin_width, noise*10**6
        
        i+=1
        
    return bin_centers, error_bar
        
        
        
if __name__ == "__main__":
    
    
    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")

    user_input["Star"]["R_Star"] = 0.6*R_Sun
    user_input["Star"]["T"] = 4000.
    user_input["Planet"]["R_Planet"] = 1*R_Earth
    user_input["Planet"]["R_Atmosphere"] = 200*1000
    user_input["Telescope"]["Aperture"] = 6.5
    user_input["Telescope"]["Quantum_Efficiency"] = 1
    user_input["Telescope"]["Distance"] = 3*Psec
    user_input["Telescope"]["Duration"] = 200*3600
    user_input["Telescope"]["min_wavelength"] = 1
    user_input["Telescope"]["max_wavelength"] = 25
    user_input["Observation_Effects"]["Noise"]["Multiplier"] = 1.2
    user_input["Observation_Effects"]["bin_exponent"] = 3/2.
    user_input["Observation_Effects"]["bin_width"] = 0.01
    
    
    #user_input["Planet"]["R_Atmosphere"] = 20000*1000
    
    
    bin_centers, error_bar = determine_bin(user_input)
    
    
    """
    Mdwarf has more FarUV (magnetic field), less MUV (cooler photosphere)
    Mdwarf has 100 - 1000x less bioactive UV => minimum threshold needed
    But... Mdwarf flares! 10x more uv than the sun. (1 every 10-40 days) ~ 20% mdwarfs
    """



