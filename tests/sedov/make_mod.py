# sod set up from 
# http://cococubed.asu.edu/research_pages/sedov.shtml#

k   = 1.380658e-16     # boltzmann constant (ergs/K)
m_p = 1.67262158e-24   # mass of proton (g)

nx   = 300
rho  = 1
rsh  = 1.5
dr = (rsh)/(1.0*nx)
v    = 0
eblast = 0.851072


zin = 5
rin = (zin)*dr
vol = 4.0*3.14159/3.0*rin**3

fout = open("sedov_model.mod","w")
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(0.0) + " 1 \n")
fout.write("1.1\n")

for i in range(nx):
	r = (i+1)*dr
	rho = 1
	v   = 0.0
	if (i < zin): egas = eblast/vol
	else: egas =  eblast/vol*1e-10
	T   = egas/(k*rho/m_p)*(1.4 - 1.0)

	line = "%16.10e %15.8e %16.8e %16.8e 1.0\n" % (r,v,rho,T)
	fout.write(line)

fout.close()

