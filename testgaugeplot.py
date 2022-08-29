import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler

# for i in range(12,13):
#     gauge=np.loadtxt('_output/SL_'+ "{:04d}".format(i) +'/gauge05832.txt',skiprows = 3)
#     plt.plot(gauge[:,1],gauge[:,5])
#     y=(gauge[0,5]-gauge[0,2]) - (gauge[1,5]-gauge[1,2])
#     plt.scatter(i, y)
# plt.show()

# data=np.load('/mnt/data/nragu/Tsunami/ML_Tohoku/ML_Slab/_riku1024coded/_data/riku.npy')
# scaler0 = StandardScaler().fit(data[:,0,:])
# standardizeddata0 = scaler0.transform(data[:,0,:])

# scaler1 = StandardScaler().fit(data[:,1,:])
# standardizeddata1 = scaler1.transform(data[:,1,:])

# plt.plot(standardizeddata0[400])
# plt.plot(standardizeddata1[400])

# data256=np.load('/mnt/data/nragu/Tsunami/ML_Tohoku/EGU/_wave/_self/_data/riku.npy')
data1024=np.load('/mnt/data/nragu/Tsunami/ML_Tohoku/EGU/_wave/_riku/_data/riku.npy')
# plt.plot(data256[400,1,:])
plt.plot(data1024[400,1,:])

footprint=np.loadtxt('/mnt/data/nragu/Tsunami/ML_Tohoku/ML_Footprint/fgmax3_flooded.csv', delimiter=',')
plt.plot(footprint[:,400])