nx    = 200
vel   = 1e9
texp  =  0
temp  = 1e4;


dpulse = 0.02
delta =  1
r0   = 1e5
r1   = r0 - delta
r2   = r0 + delta
dr   = 2*delta/(1.0*nx)
rho  = 1.0/dr/2.0


fout = open("constant.mod","w")

fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(r1) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
	r = r1+(i+1)*dr
	if (r < r0 + dpulse and r > r0 - dpulse):
		temp = 1e4
	else: temp = 0
	line = "%14.8e %10.4e %10.4e %10.4e 1.0\n" % (r,vel,rho,temp)
	fout.write(line)

fout.close()

