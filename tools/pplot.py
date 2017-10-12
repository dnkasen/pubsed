#!/usr/bin/env python
#
# Python script for plotting up plot files
#

import optparse
import matplotlib.pyplot as py
import numpy as np
import h5py

parser = optparse.OptionParser()
parser.add_option("-l",dest="log")
parser.add_option("-t",dest="type")

parser.add_option("--xr",dest="xrange")
parser.add_option("--yr",dest="yrange")
parser.add_option("--line",dest="linestyle")

(opts, args) = parser.parse_args()

if (opts.linestyle): line = opts.linestyle
else: line = ''

if (opts.type):
    type = opts.type
else:
    type = "rho,T_gas,T_rad,velr,velz,e_nuc_dep,e_nuc_emit"

nfiles = len(args) 

# get dimensions
fin = h5py.File(args[0],'r')
rho = np.array(fin['rho'])
ndims  = len(rho.shape)

py.ion()

##################### 1D spherical
if (ndims == 1):

  for i in range(nfiles):
    print args[i]

  if (opts.log == 'y'): py.yscale('log')
  if (opts.log == 'x'): py.xscale('log')
  if (opts.log == 'xy' or opts.log == 'yx'): 
    py.yscale('log')
    py.xscale('log')

  if (opts.xrange):
    xx = opts.xrange.split(',')
    x1 = float(xx[0])
    x2 = float(xx[1])
    py.xlim((x1,x2))

  if (opts.yrange):
    xx = opts.yrange.split(',')
    x1 = float(xx[0])
    x2 = float(xx[1])
    py.ylim((x1,x2))

  py.show()
  j = raw_input('press any key to exit > ')

##################### 2D cylndrical
if (ndims == 2):

  for thistype in type.split(','):

    data = np.array(fin[thistype])

    title = thistype
    if (opts.log):
      data = np.log10(data)
      title = 'log10 ' + title

    data = np.swapaxes(data,0,1)
    py.matshow(data)
    py.colorbar()
    py.title(title)
    py.show()
    j = raw_input("press any key to continue >")
    py.clf()


##################### 3D cartesian
if (ndims == 3):

  for thistype in type.split(','):

    data = np.array(fin[thistype])

    title = thistype
    if (opts.log):
      data = np.log10(data)
      title = 'log10 ' + title

    imid = data.shape[2]/2
    arr = data[:,:,imid]
    title = title + ' (x-y plane)'

    py.matshow(arr)
    py.colorbar()
    py.title(title)
    py.show()
    j = raw_input("press any key to continue >")
    py.clf()

