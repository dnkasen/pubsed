import numpy as np

r_in  = 1e14
r_out = 1e15
nx    = 128
M     = 0.5*1.99e33

v    =   0
temp = 1e4

dr = (r_out - r_in)/(nx*1.0)
r = np.zeros(nx)
for i in range(nx):
	r[i] = r_in + dr*(i+1.0)
print len(r)

rho  = M/(4.0*3.1415*r_in**2*(r_out - r_in))*(r/r_in)**(-2.0)
fout = open("powerlaw.mod","w")

fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(r_in) + "\t" + str(0) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
    line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r[i],v,rho[i],temp)
    fout.write(line)

fout.close()

