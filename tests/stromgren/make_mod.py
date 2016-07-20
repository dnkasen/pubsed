m_p = 1.67262158e-24   # mass of proton (g)
kpc    = 3.08e21       # kiloparsec in cm


nx   = 150
temp = 1e4;
rho  = 1e-3*m_p
rmax = 25*kpc

fout = open("constant.mod","w")

fout.write("1D_sphere standard\n")
fout.write(str(nx) + "\t" + str(0) + "\t" + str(0) + " 1 \n")
fout.write("1.1\n")

for i  in range(nx):
    r = rmax*(i+1)/(nx*1.0)
    v = 0
    line = "%10.4e %10.4e %10.4e %10.4e 1.0\n" % (r,v,rho,temp)
    fout.write(line)

fout.close()

