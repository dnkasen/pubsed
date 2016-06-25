nx   = 100
vmax = 2e9
texp = 10*3600.0*24
temp = 6e3;
rho  = 1e-18
dv   = vmax/nx

fout = open("constant_homologous.mod","w")
fout.write("1D_sphere SNR\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")
for i  in range(nx):
    v = (i+1.0)*dv
    line = "%10.4e %10.4e %10.4e 1.0\n" % (v,rho,temp)
    fout.write(line)
fout.close()


fout = open("constant_standard.mod","w")
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")
for i  in range(nx):
    v = (i+1.0)*dv
    r = v*texp
    line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,v,rho,temp)
    fout.write(line)
fout.close()
