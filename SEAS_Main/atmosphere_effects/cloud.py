#!/usr/bin/env python
#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Module used to simulate cloud code

"""


import os
import sys
import math
import cmath
import numpy as np
import scipy.special as special
import matplotlib.pyplot as plt



class Grey_Cloud():
    
    def __init__(self):
        pass





class Cloud_Simulator():
    
    def __init__(self,lambd,radius):
        
        self.lambd  = lambd
        self.radius = radius
    
    
    def mie_abcd(self,m,x):
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
        sinz = (cmath.sin(z))/z
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
    
    def Mie(self,m,x):
        
        if np.any(x)==0: #To avoid a singularity at x=0
            return 0,0,0
        
        
        nmax=round(2+x+(4*x**(1./3.)))
        n1=int(nmax-1);
        n = np.arange(1,nmax+1,1)
        cn=2*n+1
        c1n=np.true_divide((np.multiply(n,(n+2))),(n+1))
        c2n=np.true_divide((np.true_divide(cn,n)),(n+1))
        x2=x*x
        
        f=self.mie_abcd(m,x)
        
        anp=(f[0,:]).real
        anpp=(f[0,:]).imag
        bnp=(f[1,:]).real
        bnpp=(f[1,:]).imag
        g1=np.empty([4,int(nmax)])   # displaced numbers used fo
        g1[0,0:int(n1)]=anp[1:int(nmax)] # asymmetry parameter, p. 120
        g1[1,0:int(n1)]=anpp[1:int(nmax)]
        g1[2,0:n1]=bnp[1:int(nmax)]
        g1[3,0:n1]=bnpp[1:int(nmax)]
        
        
        dn=np.multiply(cn,(anp+bnp))
        q=sum(dn);
        qext=2*q/x2;
        en=np.multiply(cn,(np.multiply(anp,anp)+
                           np.multiply(anpp,anpp)+
                           np.multiply(bnp,bnp)+
                           np.multiply(bnpp,bnpp)))
        q=sum(en);
        qsca=2*q/x2;
        qabs=qext-qsca;
        fn=np.multiply((f[0,:]-f[1,:]),cn)
        gn=(-1)**n;
        f[2,:]=np.multiply(fn,gn)
        q=sum(f[2,:])
        qb=q*(1/q)/x2
        
        
        asy1=np.multiply(c1n,(np.multiply(anp,g1[0,:])+
                              np.multiply(anpp,g1[1,:])+
                              np.multiply(bnp,g1[2,:])+
                              np.multiply(bnpp,g1[3,:])))
        asy2=np.multiply(c2n,(np.multiply(anp,bnp)+
                              np.multiply(anpp,bnpp)))
        
        asy=4/x2*sum(asy1+asy2)/qsca;
        qratio=qb/qsca;
        
        
        return qext, qsca, qabs
    
    
    def spect(self, index):
        #rads can be a particle distribution
        countlam = 0
        index_of_refrac_medium = 1.0
        mat_abs  = np.empty([len(self.lambd),len(self.radius)])
        mat_sca  = np.empty([len(self.lambd),len(self.radius)])
        mat_qext = np.empty([len(self.lambd),len(self.radius)])
        for l in self.lambd:
            countrad = 0
            for r in self.radius:
                x =(2.0*np.pi*index_of_refrac_medium*r)/l
                qext,qsca,qabs = self.Mie(x, index)
                mat_abs[countlam,countrad] = qabs
                mat_sca[countlam,countrad] = qsca
                mat_qext[countlam,countrad] = qext
                countrad=countrad+1
            countlam=countlam+1
        
        return mat_abs,mat_sca,mat_qext,x    
    
    
    def GetSigma(self,mat_sca):
        
        mat_sigma = np.empty([len(self.lambd),len(self.radius)])
        countlam=0
        for l in self.lambd:
            countrad = 0
            for r in self.radius:
                sigma=np.pi*(r**2)*(mat_sca[countlam,countrad])
                mat_sigma[countlam,countrad] = sigma
                countrad=countrad+1
            countlam=countlam+1 
        return mat_sigma
    

    def plot(self):
        "temporary placeholder"
        
                
        plt.xlabel('Wavelength (um)')
        plt.ylabel('Sigma (cm2)')
        plt.title('Water (m = 1.33), monodisperse system')
        plt.show()
