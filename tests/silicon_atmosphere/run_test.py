import os
import pylab as py

x,f,c = py.loadtxt('out_optical_1.spec',unpack=1,skiprows=1)
lam = 3e10/x*1e8
ff = f/lam**2
b = py.interp(8000.0,lam,ff)
print b
ff = ff/b*1.2
py.plot(lam,ff)

ls,fs,es = py.loadtxt('fort.11',unpack=1)
b = py.interp(8000,ls,fs)
fs = fs/b

py.plot(ls,fs)
py.xlim(2000,10000)

py.ion()
py.show()
j = raw_input()
