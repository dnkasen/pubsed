nx   = 100
vmax = 1e9
texp = 20*3600.0*24
temp = 1e4;
rmin = 5e14
vmin = rmin/texp
rho  = 1e-20
dv   = (vmax-vmin)/nx


fout = open("vacuum.mod","w")

fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(rmin) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
    v = vmin + (i+1.0)*dv
    r = v*texp
    line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,0,rho,temp)
    fout.write(line)

fout.close()

