nx    = 200
vel   = 0.0
texp  = 0.0
temp  = 1e4;


r_pulse = 5.0
r0      = 1.0e2
rho     = 1.0
dr      = r0/(1.0*nx)

fout = open("constant.mod","w")

fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
	r = (i+1)*dr
	if (r < r_pulse):
		temp = 1e4
	else: temp = 0
	line = "%14.8e %10.4e %10.4e %10.4e 1.0\n" % (r,vel,rho,temp)
	fout.write(line)

fout.close()
