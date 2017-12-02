"""
trying to model or parametrize the VMR of water vapor

it should depend on the TP profile? (or augmented?)

ocean coverage?

exponential decay?

"""


import matplotlib.pyplot as plt
import numpy as np

from scipy import signal
from scipy.fftpack import fft, fftshift
from scipy.interpolate import interp1d

pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05]


def calculate_VMR():
    """
    
    
    All Parameters: 
        Surface VMR: VMR_0
        VMR at turn point 1: VMR_1
        VMR at turn point 2: VMR_2
        Terminator VMR: VMR_3
        Surface Pressure: P_0
        Pressure at turn point 1: P_1
        Pressure at turn point 2: P_2
        Terminator Pressure: P_3
    
    
    Dependent VMR:
        Between P_1 and P_2, the VMR is assumed to be constant, thus
        VMR_2 = VMR_1
        
        Terminator is only used as a end point, usually assumed to have reached the top of the atmosphere
        Further are truncated. So assume top of the atmosphere to be 15 scale height up, or around 120km
        which is the top of the thermosphere for earth. it's around 0.01 pa here or 0.0001 mb
    
    
    Free Parameter:
        VMR_0
        VMR_1
        P_0
        P_1
        P_2
    
    We can probably constrain the free parameters with better cloud model or biochemistry code?
    
    can we generalize this to all molecules?
    
    
    """
    
    VMR_0 = 10**-2
    
    VMR_1 = 10**-5.5
    
    VMR_2 = VMR_1
    
    VMR_3 = 10**-7
    
    P_0 = 10**5
    
    P_1 = P_0*np.e**-2
    
    P_2 = P_0*np.e**-10
    
    P_3 = P_0*np.e**-15
    
    
    # Initial decrease from sea level
    Press = np.linspace(0,np.log(P_1/P_0),101)
    PRange1 = P_0*np.e**Press
    
    a = np.log10(VMR_0/VMR_1)
    b = np.log(P_0/P_1)
    k = a/b*np.log(PRange1/P_0)
        
    VMR_Range_1 = VMR_0*10**(k)
    
    # well mixed region
    Press2 = np.linspace(np.log(P_1/P_0),np.log(P_2/P_0),101)
    PRange2 = P_0*np.e**Press2
    
    VMR_Range_2 = VMR_1*np.ones(101)
    
    
    
    # final decrease to zero at top of the atmosphere
    # H2O are broken down into H and O due to UV radiation from the star
    # photolysis: maximum wavelength needed for photolysis is ~419nm
    # Quiet M dwarfs may lack photolysis of water?
    # active M dwarfs will be active in UV though.
    Press3 = np.linspace(np.log(P_2/P_0),np.log(P_3/P_0),101)
    PRange3 = P_0*np.e**Press3
    
    m = np.log10(VMR_2/VMR_3)
    n = np.log(P_2/P_3)
    
    k = m/n*np.log(PRange3/P_2)
    
    VMR_Range_3 = VMR_2*10**(k)

    Press4 = np.linspace(np.log(P_3/P_0),np.log(10**5*np.e**-25/P_0),11)
    PRange4 = P_0*np.e**Press4
    VMR_Range4 = np.ones(11)*VMR_3
    
    PRange = np.concatenate((PRange1[:-1],PRange2[:-1],PRange3[:-1]))
    VMR_Range = np.concatenate((VMR_Range_1[:-1],VMR_Range_2[:-1],VMR_Range_3[:-1]))
    
    """
    f = interp1d(PRange, VMR_Range)
    VMR_interpolate = f(pressures)
    """
    
    plt.plot(VMR_Range, -np.log(PRange/P_0))
    #plt.plot(VMR_interpolate, -np.log(np.array(pressures)/P_0) )
    
    plt.xscale("log")
    
    
    plt.xlabel("Volume Mixing Ratio of Water Vapor")
    plt.ylabel("Scale Height from Ground")
    plt.title("Parameterized Water Vapor VMR")
    plt.show()
    
    
    


if __name__ == "__main__":
    calculate_VMR()
    
    
    
    
