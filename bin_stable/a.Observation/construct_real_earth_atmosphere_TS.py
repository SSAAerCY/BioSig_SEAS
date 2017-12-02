"""

Construct a Earth Transmission Spectra from Real TS data from the ACE experiment


In order to estimate the transmission spectra for the entire globe, 
we can assume it's the average of transmission spectra from summer and winter.
For a given time, the entire globe will have half in winter and half in summer, 
thus we can represent the global transmission spectra as the sum of summer and winter.
For such instance, we will represent North as Winter ans South as Summer.

Since the ACE data are collected for: 
    Arctic Winter
    Arctic Summer
    Mid Latitude Winter
    Mid Latitude Summer
    Tropics

we made the simple assumption to slice the earth atmosphere cross section into 8 equal size regions, with:
    Arctic Winter          : 90N-67.5N    (67.5NW - 67.5NE,                  45 degree total)
    Mid Latitude Winter x2 : 67.5N-22.5N  (67.5NW - 22.5NW, 67.5NE - 22.5NE, 90 degree total) 
    Tropics             x2 : 22.5N-22.5S  (22.5NW - 22.5SW, 22.5NE - 22.5SE, 90 degree total)
    Mid Latitude Summer x2 : 22.5S-67.5S  (22.5SW - 67.5SW, 22.5SE - 67.5SE, 90 degree total)    
    Arctic Summer          : 67.5S-90S    (67.5SW - 67.5SE,                  45 degree total)

Thus the total transmission spectra for a given "beam" of the atmosphere is

    T_layer = 1./8 (AW+AS+2*(MW+T+MS))




"""


import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.data_saver import save_txt


def generate_beam():
    
    beam = []
    start = 0
    
    while True:
        start+=4
        end = start+4
        beam.append("_%03d-%03dkm_trim.npy"%(start,end))        
        if end == 124:
            break
    
    
    
    return beam,len(beam)
    


def display_one_layer():

    inputfilepath = "../../input/absorption_data/ACE_Earth_Spectra/"

    Timer = simple_timer()

    location = {"as":"ArcticSummer", 
                "aw":"ArcticWinter",
                "mls":"MidLatitudeSummer",
                "mlw":"MidLatitudeWinter",
                "tro":"Tropics"}

    name = "_004-008km_trim.npy"

    total = np.zeros(12000)
    
    legends = []
    for i,folder in enumerate(location.keys()):
        header = "".join([location[folder],name])
        
        xdata,ydata = np.load(os.path.join(inputfilepath,folder,header))
        
        legs, = plt.plot(xdata,ydata+i+1, label=location[folder])
        legends.append(legs)
        
        if folder == "aw" or folder == "as":
            total += ydata
        else:
            total += ydata*2
    
    legs, = plt.plot(xdata,total/8.,color="k", label = "Global Average")  
    legends.append(legs)
    
    plt.legend(handles=legends)
    
    
    plt.title("ACE Atmosphere Data at 004-008km ")
    plt.xlabel("wavenumber (cm^-1)")
    plt.ylabel("Transmission Spectra")
    
    plt.show()



def main():
        
    
    
    inputfilepath = "../../input/absorption_data/ACE_Earth_Spectra/"

    Timer = simple_timer()

    location = {"as":"ArcticSummer", 
                "aw":"ArcticWinter",
                "mls":"MidLatitudeSummer",
                "mlw":"MidLatitudeWinter",
                "tro":"Tropics"}

    all_beam, beam_num = generate_beam()

    total = np.zeros(12000)    
    
    for name in all_beam:
        
        beam_total = np.zeros(12000) 
        for i,folder in enumerate(location.keys()):
            header = "".join([location[folder],name])
            xdata,ydata = np.load(os.path.join(inputfilepath,folder,header))
            
            if folder == "aw" or folder == "as":
                beam_total += ydata
            else:
                beam_total += ydata*2
        total += (beam_total/8.)
        print Timer.elapse()
    
    average = -np.log(total/beam_num)*0.001+0.01
    
    data = np.array([10000./xdata[::-1],average[::-1]]).T
    
    
    print data
    #plt.xscale("log")
    
    save_txt("../../input/absorption_data/ACE_Earth_Spectra","Test_Sim_1",data)
    
    
    
    
    
    
    plt.plot(10000./xdata,average)
    plt.title("ACE Atmosphere Data, Earth Global TS reconstruction ")
    plt.xlabel("wavenumber (cm^-1)")
    plt.ylabel("(R_P/R_S)^2")
    
    plt.show()   
    

if __name__ == "__main__":
    
    main()
    
            
            
            