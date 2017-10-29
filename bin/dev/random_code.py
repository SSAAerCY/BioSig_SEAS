from scipy.optimize import curve_fit
from numpy import *
import matplotlib.pyplot as plt

# Create a function
# ==> First encounter with *whitespace* in Python <==
def gaussian(x, a, b, c):
    val = a * exp(-(x - b)**2 / c**2)
    return val

# Generate fake data.
# Note: functions in random package, array arithmetic (exp)
n = 100
x = random.uniform(-10., 10., n)
y = exp(-(x - 3.)**2 / 4) * 10. + random.normal(0., 2., n)
e = random.uniform(0.1, 1., n)
# Note: these error bars don't reflect the distribution from which
# they were drawn! Chi^2 of the fit will be poor.

# Fit
popt, pcov = curve_fit(gaussian, x, y, sigma=e)

# Print results
print("Scale =  %.3f +/- %.3f" % (popt[0], sqrt(pcov[0, 0])))
print("Offset = %.3f +/- %.3f" % (popt[1], sqrt(pcov[1, 1])))
print("Sigma =  %.3f +/- %.3f" % (popt[2], sqrt(pcov[2, 2])))

# Plot data
plt.errorbar(x, y, yerr=e, linewidth=1, color='black', fmt=".")

# Plot model
xm = linspace(-10., 10., 100)  # 100 evenly spaced points
plt.plot(xm, gaussian(xm, popt[0], popt[1], popt[2]))



plt.show()