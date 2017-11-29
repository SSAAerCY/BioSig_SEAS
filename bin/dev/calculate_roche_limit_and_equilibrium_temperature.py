"""

I want to understand the regime in which disintegrating exoplanets will form. 
Why are some planets Hot jupiters... why are some disintegrating...

I believe this has to do with the 


roche limit vs hill sphere?


what's the albedo of the planet?


"""

import numpy as np
import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.constants import *

def calc_roche_limit(R_Planet,M_Star,M_Planet):
    
    return 1.442*R_Planet*(M_Star/M_Planet)**(1/3.)

def calc_roche_limit2(R_Star, Rho_Star, Rho_Planet):

    return 1.442*R_Star*(Rho_Star/Rho_Star)**(1/3.)


def calc_roche_limit3(M_star, Rho_Planet):
    
    return 0.8947*(M_star/Rho_Planet)**(1./3)


def equilibrium_temperature(T_star, R_star, Dist, albedo):
    
    return T_star*(1-albedo)**0.25*(R_star/2/Dist)**0.5
    

def thermal_escape_rate():
    
    pass
    
def equilibrium_temperature_at_roche():


    R_star = R_Sun*0.26
    M_star = M_Sun*0.2
    T_star = 3100.
    
    R_planet = R_Earth
    M_planet = M_Earth
    albedo = 0.3
    
    Rho_star = M_star/(4*np.pi/3*R_star**3)
    Rho_planet = M_planet/(4*np.pi/3*R_planet**3)
    
    R_roche3 = calc_roche_limit3(M_star,Rho_planet)
    print R_roche3, R_roche3/(R_Sun*0.26), R_roche3/AU  
    
    Dist = R_roche3
    
    Teq = equilibrium_temperature(T_star, R_star, Dist, albedo)
    
    print Teq



if __name__ == "__main__":
    
    equilibrium_temperature_at_roche()

    
    
    
    
    