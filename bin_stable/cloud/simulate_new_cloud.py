"""

Purpose of this code:

1. Fix the mie scattering code to be a more generalized purpose.
    1. a generalized way of loading particulates
        maybe a referencer function and a table for possible particulates?

2. Think about VMR for these particles...
 



"""



import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

from scipy.io import netcdf
from SEAS_Utils.common_utils.DIRs import Aerosol_Data
from SEAS_Main.atmosphere_effects.cloud import Physical_Cloud_Simulator_2
from SEAS_Utils.common_utils.data_loader import Particulate_Info_Loader

import numpy as np
import math
import cmath
import scipy.special as special

Pi = math.pi

def mie_abcd(m,x):
    nmax=round(2+x+(4*x**(1./3.)))
    i = 1.0j
    n = np.arange(1,nmax+1,1)
    nu = (n+0.5)
    z = np.multiply(m,x)
    m2 = np.multiply(m,m)
    sqx = np.sqrt(0.5*(math.pi)/x)
    sqz = np.sqrt(0.5*(math.pi)/z)
    bx = (np.multiply((special.jv(nu,x)),sqx))
    bz = np.multiply((special.jv(nu,z)),sqz)
    yx = np.multiply((special.yv(nu,x)),sqx)
    hx = bx+(i*yx)
    sinx = np.array((cmath.sin(x))/x)
    b1x = np.append(sinx,bx[0:int(nmax-1)])
    sinz = np.true_divide(cmath.sin(z),z)
    b1z = np.append(sinz,bz[0:int(nmax-1)])
    cosx = np.array((cmath.cos(x))/x)
    y1x = np.append(-cosx,yx[0:int(nmax-1)])
    h1x = b1x+(i*y1x)
    ax = (np.multiply(x,b1x))-(np.multiply(n,bx))
    az = (np.multiply(z,b1z))-(np.multiply(n,bz))
    ahx = (np.multiply(x,h1x))-(np.multiply(n,hx))
    m2bz = np.multiply(m2,bz)
    antop = (np.multiply(m2bz,ax))-np.multiply(bx,az)
    anbot = (m2*bz*ahx)-(hx*az)
    an = np.true_divide(antop,anbot)
    bn = (bz*ax-bx*az)/(bz*ahx-hx*az)
    cn = (bx*ahx-hx*ax)/(bz*ahx-hx*az)
    dn = m*(bx*ahx-hx*ax)/(m2*bz*ahx-hx*az)
    return np.array([an, bn, cn, dn])

def Mie(m,x):
    if np.any(x)==0: #To avoid a singularity at x=0
        result=0
    else: # This is the normal situation. 
        nmax=round(2+x+(4*x**(1./3.)))
        n1=int(nmax-1);
        n = np.arange(1,nmax+1,1)
        cn=2*n+1
        c1n=np.true_divide((np.multiply(n,(n+2))),(n+1))
        c2n=np.true_divide((np.true_divide(cn,n)),(n+1))
        x2=x*x
        f=mie_abcd(m,x)
        anp=(f[0,:]).real
        anpp=(f[0,:]).imag
        bnp=(f[1,:]).real
        bnpp=(f[1,:]).imag
        g1=np.empty([4,int(nmax)])   
        g1[0,0:int(n1)]=anp[1:int(nmax)] 
        g1[1,0:int(n1)]=anpp[1:int(nmax)]
        g1[2,0:n1]=bnp[1:int(nmax)]
        g1[3,0:n1]=bnpp[1:int(nmax)]
        dn=np.multiply(cn,(anp+bnp))
        q=sum(dn);
        qext=2*q/x2;
        en=np.multiply(cn,(np.multiply(anp,anp)+np.multiply(anpp,anpp)+np.multiply(bnp,bnp)+np.multiply(bnpp,bnpp)))
        q=sum(en);
        qsca=2*q/x2;
        qabs=qext-qsca;
        fn=np.multiply((f[0,:]-f[1,:]),cn)
        gn=(-1)**n;
        f[2,:]=np.multiply(fn,gn)
        q=sum(f[2,:])
        qb=q*(1/q)/x2
        asy1=np.multiply(c1n,(np.multiply(anp,g1[0,:])+np.multiply(anpp,g1[1,:])+np.multiply(bnp,g1[2,:])+np.multiply(bnpp,g1[3,:])))
        asy2=np.multiply(c2n,(np.multiply(anp,bnp)+np.multiply(anpp,bnpp)))
        asy=4/x2*sum(asy1+asy2)/qsca;
        qratio=qb/qsca;
        return qext,qsca,qabs,asy


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
    plt.show()
    
def calculate_cloud_xsec(info, sample=100):
    
    lam,real,imag = info

    # log normal with mean of 0 means radius mean is 1
    # log normal with mean of 1 means radius mean of 2.7 (or e)
    mean = 0
    stdev = 0.2
        
    radius = []
    for i in range(sample):
        particle_radius = float("%.2f"%np.random.lognormal(mean, stdev))*10**-6
        radius.append(particle_radius)

    sigmas = np.zeros(len(lam))
    
    for rad in radius:
        sigma = []
        for i, l in enumerate(lam):
            ref = (real[i])+(imag[i]*1.0j)
            x = np.true_divide((2.0*Pi*rad),float(l))  
            qext,qsca,qabs,g = Mie(ref,x)
            sigma_um=(Pi*(rad**2)*(qext))
            sigma_cm = np.true_divide(sigma_um,10**(8)) 
            sigma.append(sigma_cm)
        sigmas += sigma

    return lam, sigmas/sample


def plot_cross_section(particulate,lam,sigmas):
        
    plt.plot(lam,sigmas/10)
    particulate = "mgsio3"
    
    plt.xlabel("wavelength (um)")
    plt.ylabel("cross section cm^2/molecule")
    plt.title("Particulate Mie Scattering Cross Section for %s"%particulate)
    
    plt.show()
    
    
    
    
    
    
   
if __name__ == "__main__":
    
    particulate = "mgsio3"
    
    
    partical, source, filename = Particulate_Info_Loader(particulate)
    
    
    filepath = os.path.join(Aerosol_Data,source,filename)
    
    info = load_particulate(filepath,"wave")
    
    lam,sigmas = calculate_cloud_xsec(info,100)
    
    plot_cross_section(particulate,lam,sigmas)
    
    
    #plot_particulate(particulate,info)
    
    
    
    
    