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
Noise Module. Handles all noise associated works

"""
import numpy as np


def convolve_spectra():
    pass




def simple_noise():
    
    
    pass




class Noise():
    
    def __init__(self):
        pass
    
    def get_noise(self):
        pass
    
    
class Poisson_Noise(Noise):
    
    def __init__(self, poisson, shape):
        
        poissonNoise = np.random.poisson(poisson, shape).astype(float)
    
        return poissonNoise
    
class Shot_Noise(Poisson_Noise):
    """
    shot noise is most frequently observed with small currents or low light intensities that have been amplified.
    
    SNR = sqrt(N), where N is average number of event
    
    """    
    def __init__(self):
        pass
        
    


class Gaussian_Noise(Noise):

    def __init__(self, multiplier=1, length=10):
    
        self.length = length
        self.multiplier = multiplier
    
    def get_noise(self):
        return np.random.randn(self.length)*self.multiplier
        
        """
        mu, sigma = 8560, 20 # mean and standard deviation
        s = np.random.normal(mu, sigma, 1000)
        
        import matplotlib.pyplot as plt
        count, bins, ignored = plt.hist(s, 250, normed=True)
        print count, bins
        
        
        func = 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2) )
        
        plt.plot(bins, func, linewidth=2, color='r')
        plt.show()
        """
        

    
    
class Uniform_Noise(Noise):

    def __init__(self):
        
        
        return np.random.random()
    
    
class Laplace_Noise(Noise):

    def __init__(self):
        pass
    

class Lorentz_Noise(Noise):

    def __init__(self):
        pass


class Perlin_Noise(Noise):
    """
    reference https://pypi.python.org/pypi/noise/
    """
    
    def __init__(self):
        pass



class Telescope_Noise(Noise):

    def __init__(self):
        pass


    def add_jitter(self):
        pass


class JWST_Noise(Telescope_Noise):

    def __init__(self):
        pass

    def get_MIRI_noise(self):
        pass
    
    def get_NIRCam_noise(self):
        pass
    
    def get_NIRISS_noise(self):
        pass
    
    def get_NIRSpec_noise(self):
        pass






    
    