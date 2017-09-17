
import os
import sys
import numpy as np
import time
from scipy import interpolate
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_Main.atmosphere_effects.cloud import Physical_Cloud_Simulator



def calc_cloud_number_density():
    
    index  = 1.33
    lambd  = np.arange(400,30000,10)
    lambd = 10000./lambd
    radius = 1

    
    c = Physical_Cloud_Simulator(lambd,radius)
    number_density = c.calc_cloud_number_density()
    
    print number_density
    
    
def read_cloud_xsec():
    
    filename = "CROSS_ZnS.npy"
    
    
    data = np.load(filename)
    

    
    radius = data[0]
    wavelength = data[1]
    result = data[2]
    
    print wavelength
    
    xsec_list = []
    
    #for i,radius in enumerate(result):
    i = 1
    radius = result[1]
    for j,xsec in enumerate(radius):
        print result[i][j][0],
        xsec_list.append(result[i][j][0])
    print 
    #break
    
    plt.plot(wavelength,xsec_list)
    plt.show()
    
    
    
    
    
    
if __name__ == "__main__":
    #calc_cloud_number_density()
    
    read_cloud_xsec()
    
    
    
    