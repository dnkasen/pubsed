#!/usr/bin/env python
#
# Python script for plotting up sedona plot files
#

import optparse
import matplotlib.pyplot as py
import numpy as np
import h5py 

m_sun = 1.99e33
py.ion()

parser = optparse.OptionParser()
parser.add_option("-s",dest="ofile")
(opts, args) = parser.parse_args()


# read data
fin  = h5py.File(args[0],'r')
Z    = np.array(fin['Z'])
A    = np.array(fin['A'])
rho  = np.array(fin['rho'],dtype='d')
temp = np.array(fin['temp'],dtype='d')
comp = np.array(fin['comp'],dtype='d')

ndims = len(rho.shape)

#--------------------------------
#-------- 1D spherical model ----
#--------------------------------
if (ndims == 1):
	print 'hi'

#--------------------------------
#-------- 2D clyndrical model ----
#--------------------------------
if (ndims == 2):
	
	# dimensions
	nx = rho.shape[0]
	nz = rho.shape[1]

	# other quantities
	vx   = np.array(fin['vx'],dtype='d')
	vz   = np.array(fin['vz'],dtype='d')
	dr   = np.array(fin['dr'],dtype='d')

	# calculate integrated quantities
	mass  = 0
	cmass = np.zeros(len(Z))
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

	print '-----------------------------------------------------------'
	print 'dimensions; (nx,nz) = (',nx,',',nz,')'
	print ('zonesizes:  (dx,dz) = ({0:.4e}, {1:.4e})'.format(dr[0],dr[1]))
	print('mass = {0:.4e} grams ({1:.4e} m_sun)'.format(mass,mass/m_sun))
	print('ke   = {0:.4e}  ergs'.format(ke))
	for k in range(len(Z)):
		print 'elem ' +str(Z[k]) + '.' + str(A[k]),
		print(': mass = {0:.4e} grams ({1:.4e} m_sun)'.format(cmass[k],cmass[k]/m_sun))
	print '-----------------------------------------------------------'

	rho = np.swapaxes(rho,0,1)
	py.matshow(np.log10(rho))
	py.colorbar()
	py.title('log10 density')
	py.show()
	j = raw_input("press>")
	py.clf()

	temp = np.swapaxes(temp,0,1)
	py.matshow(temp)
	py.colorbar()
	py.title('temperature')
	py.show()
	j = raw_input("press>")
	py.clf()

	vx = np.swapaxes(vx,0,1)
	py.matshow(vx)
	py.colorbar()
	py.title('x=velocity')
	py.show()
	j = raw_input("press>")
	py.clf()

	vz = np.swapaxes(vz,0,1)
	py.matshow(vz)
	py.colorbar()
	py.title('z-velocity')
	py.show()
	j = raw_input("press>")
	py.clf()
	nelm = (comp.shape)[2]

	for i in range(nelm):
		carr = np.swapaxes(comp[:,:,i],0,1)
		py.matshow(carr)
		py.colorbar()
		py.title('mass fraction Z = ' + str(Z[i]) + '.' + str(A[i]))
		py.show()
		j = raw_input("press to continue >")
		py.clf()


#--------------------------------
#-------- 3D cartesian model ----
#--------------------------------
if (ndims == 3):

	# dimensions
	nx = rho.shape[0]
	ny = rho.shape[1]
	nz = rho.shape[2]

	# other quantities
	vx   = np.array(fin['vx'],dtype='d')
	vy   = np.array(fin['vy'],dtype='d')
	vz   = np.array(fin['vz'],dtype='d')
	dr   = np.array(fin['dr'],dtype='d')
	rmin = np.array(fin['rmin'],dtype='d')

	vol = dr[0]*dr[1]*dr[2]
	# calculate integrated quantities
	mass  = 0
	cmass = np.zeros(len(Z))
	ke    = 0
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				mass += rho[i,j,k]*vol
				ke   += 0.5*(vx[i,j,k]**2 + vz[i,j,k]**2)*rho[i,j,k]*vol
				for l in range(len(Z)):
					cmass[l] += rho[i,j,k]*vol*comp[i,j,k,l]

	print '-----------------------------------------------------------'
	print  'dimensions; (nx,ny,nz) = (',nx,',',ny,',',nz,')'
	print ('zonesizes:  (dx,dy,dz) = ({0:.4e}, {1:.4e},{2:.4e})'.format(dr[0],dr[1],dr[2]))
	print ('min_pos:    (x0,y0,z0) = ({0:.4e}, {1:.4e},{2:.4e})'.format(rmin[0],rmin[1],rmin[2]))
	print('mass = {0:.4e} grams ({1:.4e} m_sun)'.format(mass,mass/m_sun))
	for k in range(len(Z)):
		print 'elem ' +str(Z[k]) + '.' + str(A[k]),
		print(': mass = {0:.4e} grams ({1:.4e} m_sun)'.format(cmass[k],cmass[k]/m_sun))
	print('ke   = {0:.4e}  ergs'.format(ke))
	print '-----------------------------------------------------------'

	py.matshow(np.log10(rho[:,:,nx/2]))
	py.colorbar()
	py.title('log10 density x-y plane')
	py.show()
	j = raw_input("press to continue>")
	py.clf()

	py.matshow(np.log10(rho[:,nx/2,:]))
	py.colorbar()
	py.title('log10 density x-z plane')
	py.show()
	j = raw_input("press to continue>")
	py.clf()

	py.matshow(np.log10(rho[nx/2,:,:]))
	py.colorbar()
	py.title('log10 density y-z plane')
	py.show()
	j = raw_input("press to continue>")
	py.clf()

	#from mayavi import mlab 
	#mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(1, 1, 1))
	#mlab.contour3d(np.log10(rho),contours=20,colormap='jet' )
	#py.show()
	#j = raw_input("press to continue>")


