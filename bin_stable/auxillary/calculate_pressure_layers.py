"""
Create an array of pressures used in the simulation

"""



import numpy as np

P_surface = 100000
P_Cutoff = 0.00001

p_grid = []
P = P_surface
while P> P_Cutoff:
    p_grid.append(float("%.3g"%P))
    P = P*np.e**-1

print p_grid