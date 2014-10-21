import numpy as np
import pylab as py

data = np.loadtxt('output_optical.spec')
x = data[:,0]
y = data[:,1]
x = x/3600.0/24.0
py.plot(x,y,'o',color='black',markersize=8,markeredgewidth=2,markerfacecolor='none')

data = np.loadtxt('output_gamma.spec')
x = data[:,0]
y = data[:,1]
x = x/3600.0/24.0
py.plot(x,y,'-',color='black')

data = np.loadtxt('lucy_lc.dat')
x = data[:,0]
y = data[:,1]
py.plot(x,y,color='red',linewidth=2)


data = np.loadtxt('lucy_gr.dat')
x = data[:,0]
y = data[:,1]
py.plot(x,y,color='red',linewidth=2)


py.xlim(0,60)
py.ion()
py.show()
j = raw_input()
