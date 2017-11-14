from scipy import stats
import numpy as np


ydata = np.random.rand(10000)*10
xdata = np.random.rand(10000)*(4400-750)+750
bin_means, bin_edges, binnumber = stats.binned_statistic(xdata, ydata, bins=50, range=(750, 4400))



print bin_means
print bin_edges
print binnumber

bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2
print bin_centers

