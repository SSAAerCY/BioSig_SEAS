import numpy as np


data = [(1,2),(3,4)]


np.save("temp.npy",data)


stuff = np.load("temp.npy")

print data
print stuff

print stuff == data