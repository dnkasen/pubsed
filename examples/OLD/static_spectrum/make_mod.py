# simple script to make a 1D ejecta profile
# This does homologous expansion

nx   = 100
vmax = 2e9
texp = 50*3600.0*24
temp = 6e3;
rho0 = 0.5e-15
dv   = vmax/nx
vphot = 1e9


fout = open("input_model.mod","w")
fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(0.0) + "\t" + str(texp) + " 1 \n")
fout.write("1.1\n")
for i  in range(nx):
    v = (i+1.0)*dv
    rho = rho0*(vphot/v)**2.0
    r = v*texp
    line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,v,rho,temp)
    fout.write(line)
fout.close()


