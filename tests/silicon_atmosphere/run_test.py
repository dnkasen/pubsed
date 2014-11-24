import os
import pylab as py

x,f,c = py.loadtxt('out_optical_1.spec',unpack=1,skiprows=1)
lam = 3e10/x*1e8
f = f/lam**2
py.plot(lam,f/max(f))

ls,fs,es = py.loadtxt('fort.11',unpack=1)
fs = fs
py.plot(ls,fs/max(fs))

py.ion()
py.show()
j = raw_input()
