#!/usr/bin/env python
 
import bisect
import optparse
import numpy as np
import pylab as py
import h5py 



parser = optparse.OptionParser()
parser.add_option("-s",dest="ofile")

(opts, args) = parser.parse_args()

fin = h5py.File(args[0],'r')
py.ion()


Z = py.array(fin['Z'])
A = py.array(fin['A'])

rho = py.array(fin['rho'])
rho = py.swapaxes(rho,0,1)
py.matshow(py.log10(rho))
py.colorbar()
py.title('log10 density')
py.show()
j = raw_input("press>")
py.clf()

temp = py.array(fin['temp'])
temp = py.swapaxes(temp,0,1)
py.matshow(temp)
py.colorbar()
py.title('temperature')
py.show()
j = raw_input("press>")
py.clf()

comp = py.array(fin['comp'])
nelm = (comp.shape)[2]
for i in range(nelm):
	carr = py.swapaxes(comp[:,:,i],0,1)
	py.matshow(carr)
	py.colorbar()
	py.title('mass fraction Z = ' + str(Z[i]) + '.' + str(A[i]))
	py.show()
	j = raw_input("press>")
	py.clf()

