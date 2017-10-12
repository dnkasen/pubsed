#!/usr/bin/env python
 
import bisect
import optparse
import numpy as np
import pylab as py
import h5py 

#cm = py.get_cmap('rainbow')

parser = optparse.OptionParser()
parser.add_option("-x",dest="cols")
parser.add_option("-s",dest="ofile")
parser.add_option("-l",dest="log")
parser.add_option("-n",action="store_true", dest="norm")
parser.add_option("-t",dest="plottime",default=20.0)
parser.add_option("-m",action="store_true", dest="max")
parser.add_option("--xr",dest="xrange")
parser.add_option("--yr",dest="yrange")
parser.add_option("--yt",dest="ytitle")
parser.add_option("--xt",dest="xtitle")
parser.add_option("--title",dest="title")
parser.add_option("--math",dest="math")
parser.add_option("--line",dest="linestyle")
parser.add_option("--legend",dest="legend")
parser.add_option("--offset",dest="offset")
parser.add_option("--shift",dest="shift")
parser.add_option("--skip",dest="skip")
parser.add_option("--velplot",dest="velplot")
parser.add_option("--nu",action="store_true",dest="nu")
parser.add_option("--secs",action="store_true",dest="secs")
parser.add_option("--bol",action="store_true",dest="plotbol")


(opts, args) = parser.parse_args()
n = len(args)

if (opts.linestyle): line = opts.linestyle
else: line = ''

if (opts.skip): skip = int(opts.skip)
else: skip = 0


for i in range(n):

  if (".h5" in args[i]):
    fin = h5py.File(args[i],'r')
    nu    = py.array(fin['nu'])
    times = py.array(fin['time'])
    Lnu   = py.array(fin['Lnu'])
  else:
    nu,Lnu,count = py.loadtxt(args[i],skiprows=1,unpack=1)
    times = []

  # calculate bolometric LC
  nt  = len(times)
  nnu = len(nu)
  bol = np.zeros(nt)
  for it in range(nt):
    bol[it] = py.trapz(Lnu[it,:],x=nu)

  if (len(times) == 0):
    if (not opts.nu):
      Lnu = Lnu*nu**2/2.99e10/1e8
      nu  = 2.99e10/nu*1e8
      py.plot(nu,Lnu)
  else:
    if (not opts.secs): 
      times = times/3600.0/24.0
    if (opts.plotbol):
      py.plot(times,bol)
    else:
      py.xlabel('frequency (Hz)')
      if (not opts.nu):
        Lnu = Lnu*nu**2/2.99e10/1e8
        nu  = 2.99e10/nu*1e8
        py.xlabel('wavelength (angstroms)')
        
      for thist in (opts.plottime).split(','):
        indt = bisect.bisect(times,float(thist))
        print 'time = ',thist,indt
        py.plot(nu,Lnu[indt,:])
        
      py.ylabel('specific luminosity')
      py.title('spectrum at time t = ' + str(opts.plottime))

######################


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

if (not opts.plotbol): py.xlim(1000,3e4)

if (opts.yrange):
   xx = opts.yrange.split(',')
   x1 = float(xx[0])
   x2 = float(xx[1])
   py.ylim((x1,x2))

if (opts.xtitle):  py.xlabel(opts.xtitle)
if (opts.ytitle):  py.ylabel(opts.ytitle)
if (opts.title):   py.title(opts.title)

if (opts.legend):
   if (opts.legend == '0'): names = args
   else: names = opts.legend.split(',')
   py.legend(names)

if opts.ofile: py.savefig(opts.ofile)

py.ion()
py.show()
j = raw_input('press > ')
