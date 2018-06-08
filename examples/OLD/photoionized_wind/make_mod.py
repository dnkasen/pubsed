import numpy as np
import constants as pc

m_sun = 1.99e33
sec_to_year = 1.0/3600.0/24.0/365.0
pi = 3.14159

nx    = 100
mdot  = 1e-3*(m_sun*sec_to_year)
rin   = 1e15
rout  = 1e16
vel   = 1e8
temp  = 1.0e4
texp  = 0.0
dr = (rout-rin)/(nx*1.0)

fout  = open("wind.mod","w")

fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(rin) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
	r = rin + (i+1)*dr
	rho = mdot/(4.0*pi*r**2*vel)
	line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,vel,rho,temp)
	fout.write(line)

fout.close()

