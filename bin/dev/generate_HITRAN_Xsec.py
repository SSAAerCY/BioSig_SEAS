"""

The purpose of this is to generate new HITRAN cross sections from the 2016 data


We need extremely hot gases... 

we also need really hot mie scattering 


"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Aux.cross_section.cross_section_calculator as csc
#import SEAS_Utils.common_utils.data_plotter as plt
from SEAS_Utils.common_utils.DIRs import Temp_DIR, HITRAN_Lines
from SEAS_Utils.common_utils.timer import simple_timer

import matplotlib.pyplot as plt
import numpy as np

def calculate_extremely_hot_cross_section():

    Timer = simple_timer()
    molecule    = "HBr"
    
    d_path      = os.path.join(HITRAN_Lines,molecule)
    r_path      = Temp_DIR
    component   = [1,1,1]
    numin       = 400
    numax       = 30000
    step        = 1
    P           = 10.
    T           = 1300.
    legends = []
    
    result_data = []
    
    if os.path.isfile("temp.npy"):
        nu,result_data = np.load("temp.npy")
        
        for i,data in enumerate(result_data[::-1]):
            a, = plt.plot(10000./nu,data,label=str([300,800,1300,1800,2300][::-1][i]))
            legends.append(a)
    else:
        for T in [300,800,1300,1800,2300]:
        
            
            calc = csc.cross_section_calculator(d_path,r_path,molecule,component,numin,numax,step,P,T)
    
            data = calc.read_data()
            
            nu, coef = calc.personal_calculator(data)
            
            k = 1.38*10**-23
            n = P/(k*T)
            l = 10
            a, = plt.plot(10000./nu,coef,label=str(T))
            legends.append(a)
            #plt.plot(10000./nu,np.exp(-n*coef*l))
            
            result_data.append(coef)
        
        np.save("temp.npy",[nu,result_data])
    
    
    ax = plt.gca()
    
    plt.xscale("log")
    plt.yscale("log")
    
    plt.legend(legends)
    
    plt.title("Water Vapor at different Temperature")
    plt.xlabel("Wavelegnth (cm^-1)")
    plt.ylabel("cm^2/molecule")    
    plt.show()
    
        
if __name__ == "__main__":
    calculate_extremely_hot_cross_section()
    
    