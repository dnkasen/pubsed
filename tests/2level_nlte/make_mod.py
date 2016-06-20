nx   = 100
temp = 5e4;
texp = 20.0;
rmin = 3e13
rmax = 1e14
rho  = 1e-20
vel  = 0

dr = (rmax - rmin)/(1.0*nx)
fout = open("vacuum.mod","w")

fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(rmin) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
	r = rmin + (i+1)*dr
	line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,vel,rho,temp)
	fout.write(line)

fout.close()

