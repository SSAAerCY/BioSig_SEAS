"""
Code related with using the pandexo web tool

learning place to figure out what each component is doing and what needs to be calculated

"""
import os
import sys
import numpy as np

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.common_utils as utils
import SEAS_Utils.common_utils.configurable as config
from SEAS_Utils.common_utils.constants import *




def calc_period(M_star, a):
    
    G = 6.67*10**-11
    
    return np.sqrt(4*np.pi**2*a**3/G/M_star)

def calc_transit_duration(M_star, R_star, a, R_planet):
    

    P = calc_period(M_star, a)
    b = a*np.cos(i)/R_star
    
    param = np.sqrt((R_star+R_planet)**2-(b*R_star)**2)/a
    
    T_duration = P/np.pi*np.arcsin(param)
    
    return T_duration

def ratio_in_transit(P, T_duration):
    
    return T_duration/P
    


def calc_stellar_params():
    
    pass

def M_Dwarf_Mass_Radius_Relation():
    """
    This also depends on metallicity of the star.
    """
    
    pass    




if __name__ == "__main__":

    M_sun   = 2*10**30
    R_sun   = 695000000.
    R_earth = 6400000.
    a_earth = 1.5*10**11
    i       = (90./180)*np.pi


    user_input = config.Configuration("../../bin_stable/a.Main/user_input_dev.cfg")
    stellar_input = user_input = config.Configuration("../../bin_stable/a.Main/stellar_info.cfg")

    star = "M0V"
    

    mass = utils.to_float(stellar_input[star]["Mass"])*M_sun
    radi = utils.to_float(stellar_input[star]["Radi"])*R_sun
    a    = utils.to_float(stellar_input[star]["Habz"])*a_earth
    
    T_dur = calc_transit_duration(mass, radi, a, R_earth)
    P = calc_period(mass,a)
    ratio = ratio_in_transit(P,T_dur)
    
    print T_dur, T_dur/3600, 200/(T_dur/3600.), 200/(T_dur/3600.)*(P/86400)/365., ratio
    
    
    