# sod set up from 
# http://cococubed.asu.edu/code_pages/exact_riemann.shtml

nx   = 300
rho  = 1
r0   = 1e3
rsh  = 5
dr = (rsh)/(1.0*nx)
v    = 0

k   = 1.380658e-16     # boltzmann constant (ergs/K)
m_p = 1.67262158e-24   # mass of proton (g)


fout = open("sod_model.mod","w")
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(r0) + "\t" + str(0.0) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
	r = (i+1)*dr

	if (r < 2): 
		rho = 10.0
		p   = 100.0
	else: 
		rho = 1.0
		p   = 1.0

	T = p/k/(rho/m_p)
	line = "%16.10e %15.8e %16.8e %16.8e 1.0\n" % (r0 + r,v,rho,T)
	fout.write(line)

fout.close()

