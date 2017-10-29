
from __future__ import division
import os
import sys
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

"""
%matplotlib inline
%precision 4
"""
plt.style.use('ggplot')

from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as stats
from functools import partial

np.random.seed(1234)


def analytic():
    n = 100
    h = 61
    p = h/n
    rv = st.binom(n, p)
    mu = rv.mean()
    
    a, b = 10, 10
    prior = st.beta(a, b)
    post = st.beta(h+a, n-h+b)
    ci = post.interval(0.95)
    
    thetas = np.linspace(0, 1, 200)
    plt.figure(figsize=(9,6 ))
    plt.style.use('ggplot')
    plt.plot(thetas, prior.pdf(thetas), label='Prior', c='blue')
    plt.plot(thetas, post.pdf(thetas), label='Posterior', c='red')
    plt.plot(thetas, n*st.binom(n, thetas).pmf(h), label='Likelihood', c='green')
    plt.axvline((h+a-1)/(n+a+b-2), c='red', linestyle='dashed', alpha=0.4, label='MAP')
    plt.axvline(mu/n, c='green', linestyle='dashed', alpha=0.4, label='MLE')
    plt.xlim([0, 1])
    plt.axhline(0.3, ci[0], ci[1], c='black', linewidth=2, label='95% CI')
    plt.xlabel(r'$\theta$', fontsize=14)
    plt.ylabel('Density', fontsize=16)
    plt.legend()
    plt.show()


def numerical():
    """
    One advantage of this is that the prior does not have to be conjugate 
    (although the example below uses the same beta prior for ease of comaprsion), 
    and so we are not restricted in our choice of an approproirate prior distribution. 
    For example, the prior can be a mixture distribution or estimated empirically from data. 
    The disadvantage, of course, is that this is computationally very expenisve when 
    we need to esitmate multiple parameters, since the number of grid points grows 
    as O(n^d), wher nn defines the grid resolution and dd is the size of θθ.
    
    """
    
    n = 100
    h = 61
    p = h/n
    rv = st.binom(n, p)
    mu = rv.mean()
    
    a, b = 10, 10
    prior = st.beta(a, b)
    post = st.beta(h+a, n-h+b)
    ci = post.interval(0.95)

    thetas = np.linspace(0, 1, 200)
    prior = st.beta(a, b)
    
    post = prior.pdf(thetas) * st.binom(n, thetas).pmf(h)
    post /= (post.sum() / len(thetas))
    
    plt.figure(figsize=(9,6))
    plt.plot(thetas, prior.pdf(thetas), label='Prior', c='blue')
    plt.plot(thetas, n*st.binom(n, thetas).pmf(h), label='Likelihood', c='green')
    plt.plot(thetas, post, label='Posterior', c='red')
    plt.xlim([0, 1])
    plt.xlabel(r'$\theta$', fontsize=14)
    plt.ylabel('Density', fontsize=16)
    plt.legend();
    plt.show()



def test_pymc():
    pass



if __name__ == "__main__":
    numerical()
    test_pymc()