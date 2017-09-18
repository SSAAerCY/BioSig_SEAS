import numpy as np
from scipy.special import wofz
import pylab

def gaussian(x, mu, sig):
    return np.exp(-((x - mu)/sig)**2/2)/(sig*(np.sqrt(2.*np.pi)))

def G(x, alpha):
    """ Return Gaussian line shape at x with HWHM alpha """
    return np.sqrt(np.log(2) / np.pi) / alpha* np.exp(-(x / alpha)**2 * np.log(2))

def L(x, gamma):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return gamma / np.pi / (x**2 + gamma**2)

def V(x, alpha, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)


alpha, gamma = 10,10
x = np.linspace(-10,10,1000)
pylab.plot(x, G(x, alpha), ls=':', c='b', label='Gaussian')
pylab.plot(x, L(x, gamma), ls='--', c='k', label='Lorentzian')
pylab.plot(x, V(x, alpha, gamma), c='r', label='Voigt')

pylab.legend()


pylab.show()