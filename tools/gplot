#!/usr/bin/env python
#
# Python script for plotting up grid files
#

import optparse
import pylab as py
import h5py

parser = optparse.OptionParser()
parser.add_option("-l",dest="log")
parser.add_option("-z",dest="zone")
parser.add_option("-t",dest="type")
parser.add_option("-n",action="store_true",dest="nuplot")

parser.add_option("--xr",dest="xrange")
parser.add_option("--yr",dest="yrange")
parser.add_option("--line",dest="linestyle")

(opts, args) = parser.parse_args()

if (opts.linestyle): line = opts.linestyle
else: line = ''

if (opts.zone):
    zname = 'zonedata/' + str(opts.zone)
else:
    zname = 'zonedata/0'

if (opts.type):
    type = opts.type
else:
    type = 'opacity'

print 'zone = ' + zname + '; type = ' + type

n = len(args) 
for i in range(n):
    fin = h5py.File(args[i],'r')
    nu = py.array(fin['nu'])
    lam = 3e10/nu*1e8

    if (type == 'Snu'):
      opac = py.array(fin[zname + '/opacity'])
      emis = py.array(fin[zname + '/emissivity'])
      vals = emis/opac
    else:
      vals = py.array(fin[zname + '/' + type])
    if (opts.nuplot): py.plot(nu,vals,line)
    else: py.plot(lam,vals,line)
    fin.close()

if (opts.nuplot): py.xlabel('frequency (Hz)',size=15)
else: py.xlabel('wavelength (angstroms)',size=15)

if (type == 'opacity'): py.ylabel('opacity (cm^2/g)',size=15)
if (type == 'emissivity'): py.ylabel('emissivity',size=15)
if (type == 'Jnu'): py.ylabel('radiation field Jnu (ergs/s/Hz/cm^2/str)',size=15)

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

print vals

py.ion()
py.show()
j = raw_input('press > ')
