"""

understanding log normal distribution

"""

import numpy as np
import matplotlib.pyplot as plt


mu = 1.
sigma = 0.5
s = np.random.lognormal(mu, sigma)

print s


"""
count, bins, ignored = plt.hist(s, 100, normed=True, align='mid')

plt.axis("tight")
plt.show()
"""