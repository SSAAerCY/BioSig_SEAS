import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from scipy.io import netcdf
from SEAS_Utils.common_utils.DIRs import HITRAN_Mineral
from SEAS_Main.atmosphere_effects.cloud import Physical_Cloud_Simulator_2

def load_particulate(filename,output="wave"):
    
    info = netcdf.netcdf_file(filename, 'r')
    
    text = info.variables["text"].data
    ri   = info.variables["ri"].data
    rn   = info.variables["rn"].data
    wave = info.variables["wavelength"].data
    wcm  = info.variables["wcm"].data
    lenx = info.variables["nlines"].data
    
    info.close()

    new_x = []
    new_rn = []
    new_ri = []
    
    if output=="wave":
        for i,dat in enumerate(wave):
            if dat >= 25:
                break
            new_x.append(wave[i])
            new_rn.append(rn[i])
            new_ri.append(ri[i])
            
    if output=="wcm":
        for i,dat in enumerate(wcm):
            print dat
            if dat <= 500:
                break
            new_x.append(wcm[i])
            new_rn.append(rn[i])
            new_ri.append(ri[i])
                
    return new_x,new_rn,new_ri

def plot_particulate(particulate,data):
    
    x,n,i = data
    
    fig, ax1 = plt.subplots()
    ax1.set_title("Index of Refraction (n,i) for %s"%particulate)
    ax1.plot(x,n,label="Real",color="b")
    ax1.set_xlabel("wavelength (um)")
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('Index of Refraction n', color='b')
    ax1.tick_params('y', colors='b')
    
    ax2 = ax1.twinx()
    ax2.plot(x,i,label="Imag",color="r")
    ax2.set_ylabel('Index of Refraction i', color='r')
    ax2.tick_params('y', colors='r')
    
    fig.tight_layout()
    plt.show
    
def calculate_cloud_xsec(info):
    
    x,n,i = info
    
    cloud = Physical_Cloud_Simulator_2()
    
   
if __name__ == "__main__":
    
    particulate = "mgsio3"
    
    filename = os.path.join(HITRAN_Mineral,"jager_mgsio3.nc")
    
    info = load_particulate(filename,"wcm")
    
    #plot_particulate(particulate,info)
    