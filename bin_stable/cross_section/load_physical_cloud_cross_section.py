
import os
import sys
import numpy as np
import time
from scipy import interpolate
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))


from SEAS_Main.atmosphere_effects.cloud import File_Cloud_Simulator, Physical_Cloud_Simulator
from SEAS_Utils.common_utils.data_loader import two_column_file_loader


def test_1():

    legends = []
    for radius in [1,2,3]:
        index  = 1.33
        lambd  = np.arange(400,30000,10)
        lambd = 10000./lambd
        #radius = 3
        
        Simulator = File_Cloud_Simulator(lambd, radius)
        
        n = Simulator.calc_cloud_number_density(particle_radius=radius*10**-4)
        wav, xsec = Simulator.get_cloud_cross_section("CROSS_ZnS.npy")
        
        
        lag, = plt.plot(wav, xsec, label="R=%s um"%radius)
        legends.append(lag)
        
    ax = plt.gca()
    #ax.set_yscale('log')
    plt.legend(handles = legends)
    plt.title("ZnS Cloud Particle Cross Section with Varying Radius")
    plt.xlabel("Wavelength (um)")
    plt.ylabel("Intensity (cm^2/particle)")

    
    plt.show()

def test_2():
    
    nfile = "../../bin_stable/a.Transmission/cloudy/ZnS_n.txt"
    kfile = "../../bin_stable/a.Transmission/cloudy/ZnS_k.txt"
    
    wavelength, n = two_column_file_loader(nfile)
    wavelength, k = two_column_file_loader(kfile)    


    
    legends = []
    for i in [100000,10000,1000,100,10,1]:

        
        P = i
        T = 300
        BoltK = 1.38*10**-23
        R = 6.022*10**23
        
        air_density = P/(BoltK*T)/R*29/100**3 #g/cm^3
        particle_mixing_ratio = 4.62*10**-6 # %
        particle_density = air_density*particle_mixing_ratio # g of zns in /cm^3
        total_unit_mass = particle_density*1
        ZnS_density = 4.09 #g/cm^3
        
        total_unit_mass = particle_density*1
        mean = 0.5#np.log(2.0) #arb
        stdev = 0.1 #arb
        radrange = []
        while True:
            particle_radius = float("%.2f"%np.random.lognormal(mean, stdev)) #cm
            unit_particle_mass = ZnS_density*4/3*np.pi*(particle_radius*10**-4)**3 #g
        
            total_unit_mass -= unit_particle_mass
            if total_unit_mass <0:
                break
            
            radrange.append(particle_radius)
        
        print len(radrange)
    
        Simulator = Physical_Cloud_Simulator(wavelength[:294],radrange)
        ZnS_abs,ZnS_sca,ZnS_qext,ZnS_x = Simulator.spect(n[0:294], k[0:294])
        ZnS_cross = Simulator.GetSigma(ZnS_qext,unit="cm")
        
        Total_xsec = np.zeros(294)
        for xsec in ZnS_cross.T:
            Total_xsec += xsec     
        
        l = 100
        

        
        #plt.plot(wavelength[:294],np.e**(-Total_xsec*l))
        legs, = plt.plot(wavelength[:294],Total_xsec, label="%spa"%P)
        legends.append(legs)
        
    
    plt.legend(handles=legends)
        
    plt.xlabel("Wavelength, (um)")
    plt.ylabel("Transmission")
    #plt.title("Transmission Spectra from ZnS Clouds with radius=50nm and pathlength=1m")        
    plt.title("ZnS clouds cross section with radius=50nm ") 
    plt.show()
    
def test_3():
    
    nfile = "../../bin_stable/a.Transmission/cloudy/ZnS_n.txt"
    kfile = "../../bin_stable/a.Transmission/cloudy/ZnS_k.txt"
    
    wavelength, n = two_column_file_loader(nfile)
    wavelength, k = two_column_file_loader(kfile)    


    
    legends = []
    #for i in [100000,10000,1000,100,10,1]:
    for i in [50,500,1000,2000,3000]:
        
        P = 100000
        T = 300
        BoltK = 1.38*10**-23
        R = 6.022*10**23
        
        air_density = P/(BoltK*T)/R*29/100**3 #g/cm^3
        particle_mixing_ratio = 4.62*10**-6 # %
        particle_density = air_density*particle_mixing_ratio # g of zns in /cm^3
        total_unit_mass = particle_density*1
        ZnS_density = 4.09 #g/cm^3
        
        total_unit_mass = particle_density*1*i/50
        mean = i/1000.
        stdev = 0.1 #arb
        radrange = []
        for j in range(100):
            particle_radius = float("%.2f"%np.random.lognormal(mean, stdev)) #cm
            radrange.append(particle_radius)
        print "jere"
    
        Simulator = Physical_Cloud_Simulator(wavelength[:294],radrange)
        ZnS_abs,ZnS_sca,ZnS_qext,ZnS_x = Simulator.spect(n[0:294], k[0:294])
        ZnS_cross = Simulator.GetSigma(ZnS_qext,unit="cm")
        
        Total_xsec = np.zeros(294)
        for xsec in ZnS_cross.T:
            Total_xsec += xsec     
        
        l = 100
        

        
        #plt.plot(wavelength[:294],np.e**(-Total_xsec*l))
        legs, = plt.plot(wavelength[:294],Total_xsec, label="%snm"%i)
        legends.append(legs)
        
    
    plt.legend(handles=legends)
    plt.yscale("log")
    plt.xlabel("Wavelength, (um)")
    plt.ylabel("Cross Section")
    #plt.title("Transmission Spectra from ZnS Clouds with radius=50nm and pathlength=1m")        
    plt.title("ZnS clouds cross section with various radius ") 
    plt.show()
if __name__ == "__main__":
    
    test_3()
    