
import numpy as np


data = np.genfromtxt('test_data.txt')
data[:,0] /= 1000.0 #change from nm to um

data[:,1:] *= 1.0e9 #change from 1/nm^3 to 1/um^3  (*1.0e9)

f = open('new_data.txt','w')
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
       f.write("  %.6e" %data[i,j])
    f.write("\n")
f.close()
