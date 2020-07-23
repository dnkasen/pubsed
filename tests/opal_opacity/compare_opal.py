import matplotlib.pyplot as plt
import numpy as np

logT,log_opal = np.loadtxt('opal_X0.9_Y0.1_Z0.0.dat',unpack=True,usecols=[0,8],skiprows=1)
print log_opal

T = 10**logT
opal = 10**log_opal
plt.plot(T,opal,'o')


T,op = np.loadtxt('mean_opacities.dat',unpack=True,usecols=[0,5],skiprows=1)
plt.plot(T,op,color='k')

plt.legend(['OPAL','Sedona'])
plt.xlim(4e3,5e8)
plt.ylim(0.04,4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Temperature')
plt.ylabel('Rosseland Mean opacity')
plt.show()
