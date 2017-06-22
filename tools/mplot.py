#!/usr/bin/env python
 
import optparse
import pylab as py
import h5py 



parser = optparse.OptionParser()
parser.add_option("-s",dest="ofile")

(opts, args) = parser.parse_args()

py.ion()

# read data
fin = h5py.File(args[0],'r')
Z = py.array(fin['Z'])
A = py.array(fin['A'])
rho = py.array(fin['rho'])
temp = py.array(fin['temp'])
comp = py.array(fin['comp'])
dr   = py.array(fin['dr'])
vx   = py.array(fin['vx'])
vz   = py.array(fin['vz'])

nx = rho.shape[0]
nz = rho.shape[1]
mass  = 0
cmass = py.zeros(len(Z))
ke    = 0
for i in range(nx):
	for j in range(nz):
		x = dr[0]*(i+0.5)
		z = dr[1]*(j+0.5)
		vol = 2*3.14159*dr[0]*dr[1]*x
		mass += rho[i,j]*vol
		ke   += 0.5*(vx[i,j]**2 + vz[i,j]**2)*rho[i,j]*vol
		for k in range(len(Z)):
			cmass[k] += rho[i,j]*vol*comp[i,j,k]
m_sun = 1.99e33
print '-----------------------------------------------------------'
print 'dimensions; (nx,nz) = (',nx,',',nz,')'
print ('zonesizes:  (dx,dz) = ({0:.4e}, {1:.4e})'.format(dr[0],dr[1]))
print('mass = {0:.4e} grams ({1:.4e} m_sun)'.format(mass,mass/m_sun))
print('ke   = {0:.4e}  ergs'.format(ke))
for k in range(len(Z)):
	print 'elem ' +str(Z[k]) + '.' + str(A[k]),
	print(': mass = {0:.4e} grams ({1:.4e} m_sun)'.format(cmass[k],cmass[k]/m_sun))
print '-----------------------------------------------------------'


rho = py.swapaxes(rho,0,1)
py.matshow(py.log10(rho))
py.colorbar()
py.title('log10 density')
py.show()
j = raw_input("press>")
py.clf()

temp = py.swapaxes(temp,0,1)
py.matshow(temp)
py.colorbar()
py.title('temperature')
py.show()
j = raw_input("press>")
py.clf()

nelm = (comp.shape)[2]
for i in range(nelm):
	carr = py.swapaxes(comp[:,:,i],0,1)
	py.matshow(carr)
	py.colorbar()
	py.title('mass fraction Z = ' + str(Z[i]) + '.' + str(A[i]))
	py.show()
	j = raw_input("press>")
	py.clf()

