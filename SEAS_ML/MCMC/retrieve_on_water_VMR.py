"""

a way to generate parameters for water VMR

"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from SEAS_Utils.common_utils.DIRs import Mixing_Ratio_Data
from SEAS_Utils.common_utils.data_loader import multi_column_file_loader


def generate_water_VMR(VMR_0 = 10**-2,VMR_1 = 10**-5.5,VMR_2 = 10**-5.5, VMR_3 = 10**-10,
                       P_0 = 10**5, P_1 = 10**5*np.e**-2, P_2 = 10**5*np.e**-10, P_3 = 10**5*np.e**-15):


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
    
    PRange = np.concatenate((PRange1[:-1],PRange2[:-1],PRange3[:-1], PRange4))
    VMR_Range = np.concatenate((VMR_Range_1[:-1],VMR_Range_2[:-1],VMR_Range_3[:-1],VMR_Range4))
    

    
    return PRange,VMR_Range

def generate_mixing_ratio_file(VMR,outname):
    
    
    ref_file = os.path.join(Mixing_Ratio_Data,"earth.txt")
    out_file = os.path.join(Mixing_Ratio_Data,"Retrieval",outname)
    

    data = multi_column_file_loader(ref_file,type="mixed")
    
    old_stdout = sys.stdout
    sys.stdout = open(out_file, 'w')
    
    for i in range(len(data[0])):
        
        sumAll = 0.0
        
        for j in range(len(data)):
            
            if i == 0:
                print data[j][i],
                continue
            
            
            if j == 0:
                print data[j][i],
            
            elif j == 1:
                try:
                    # multiply by 100 because converting from vmr to percent
                    print VMR[i]*100,
                    sumAll+=float(VMR[i])
                except:
                    print 0.0,
            
            elif j == len(data)-1:
                print "%.4g"%(100-sumAll),
            
            else:
                print data[j][i],
                sumAll+=float(data[j][i])

        if i == len(data[0])-1:
            break
        else:
            print
    sys.stdout = old_stdout

def plot_VMR(VMR_Range,PRange):


    plt.plot(VMR_Range*100, -np.log(PRange/P_0))
    
    plt.xscale("log")
    
    plt.xlabel("Volume Mixing Ratio of Water Vapor (Percent)")
    plt.ylabel("Scale Height from Ground")
    plt.title("Parameterized Water Vapor VMR")
    plt.show()    
    


if __name__ == "__main__":
    
    
    VMR_0 = 3*10**-3
    VMR_1 = 10**-5.5
    VMR_2 = VMR_1
    VMR_3 = 10**-5.2
    P_0 = 10**5
    P_1 = P_0*np.e**-2
    P_2 = P_0*np.e**-5
    P_3 = P_0*np.e**-23

    VMR_0 = 4*10**-2
    VMR_1 = 10**-5.5
    VMR_2 = VMR_1
    VMR_3 = 10**-7
    P_0 = 10**5
    P_1 = P_0*np.e**-2
    P_2 = P_0*np.e**-10
    P_3 = P_0*np.e**-15    
    
    VMR_0 = 2*10**-2
    VMR_1 = 10*10**-5.5
    VMR_2 = VMR_1
    VMR_3 = 10*10**-5.2
    P_0 = 10**5
    P_1 = P_0*np.e**-2
    P_2 = P_0*np.e**-5
    P_3 = P_0*np.e**-23

    
    
    PRange,VMR_Range = generate_water_VMR(VMR_0,VMR_1,VMR_2,VMR_3,P_0,P_1,P_2,P_3)
    
    pressures = [100000.0, 36800.0, 13500.0, 4980.0, 1830.0, 674.0, 248.0, 
                 91.2, 33.5, 12.3, 4.54, 1.67, 0.614, 0.226, 0.0832, 0.0306, 
                 0.0113, 0.00414, 0.00152, 0.00056, 0.000206, 7.58e-05, 2.79e-05, 1.03e-05]

    
    f = interp1d(PRange, VMR_Range)
    VMR_interpolate = f(pressures)    
    
    VMR = ["H2O"]
    for i in VMR_interpolate:
        VMR.append(float("%.4g"%i))
    
    plot_VMR(VMR_Range,PRange)
    
    name = "temp_RE.txt"
    
    #generate_mixing_ratio_file(VMR,name)
