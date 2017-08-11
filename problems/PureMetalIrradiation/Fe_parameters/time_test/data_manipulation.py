import numpy as np

size = 80 #30 #
filename ='spatial_defect_cluster_production_'+str(size)+'.txt'
shift = 25
change_num_size = 70 #25 #60

data = np.genfromtxt(filename)
cols = data.shape[1]
vcluster = range(1,int(cols/2)+1)
icluster = range(int(cols/2)+1,cols)

#data[:,icluster] = 1.1*data[:,vcluster]
bins = data.shape[0]
for vsize in range(int(cols/2)+1-change_num_size,int(cols/2)+1):
    line = np.zeros(bins)
    line[:(bins-shift)] = data[shift:,vsize]
    data[:,vsize] = line

for isize in range(cols-change_num_size,cols):
    line = np.zeros(bins)
    line[:(bins-shift)] = data[shift:,isize]
    data[:,isize] = line

for i in range(1,cols):
    print data[data[:,i].argmax(),0]
np.savetxt(filename.split('.')[0]+'_test.txt',data,delimiter=' ')
