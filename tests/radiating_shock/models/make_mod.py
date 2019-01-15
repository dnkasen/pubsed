# sod set up from 
# http://cococubed.asu.edu/code_pages/exact_riemann.shtml

nx   = 300
rho  = 7.78e-10
r0   = 8e11
rsh  = 7e10
dr = (rsh)/(1.0*nx)
v    = 0

k   = 1.380658e-16     # boltzmann constant (ergs/K)
m_p = 1.67262158e-24   # mass of proton (g)


fout = open("ensman_model.mod","w")
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(r0) + "\t" + str(0.0) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
	r = (i+1)*dr
	T = 10
	line = "%16.10e %15.8e %16.8e %16.8e 1.0\n" % (r0 + r,v,rho,T)
	fout.write(line)

fout.close()

