import pylab as py
import h5py


def compare(pdf):

	py.clf()

	fin = h5py.File('grid_00020.h5','r')
	
	r = py.array(fin['x'])
	Tgas = py.array(fin['T_gas'])
	Trad = py.array(fin['T_gas'])

	r = r - 0.5*(r[1] - r[0])

	L    = 2.45e43
	M    = 0.5*1.99e33
	c   = 2.99792458e10    # speed of light (cm/s)
	a   = 7.5657e-15       # radiation constant
	rin  = 1.0e14
	rout = 1.0e15
	kappa = 0.2

	tau = kappa*M/(4.0*3.1415*rin*rout)

	Erad = L*tau/(4.0*3.14159*c*rin**2)*((rin/r)**3 - (rin/rout)**3 + 4/tau*(rin/rout)**2)
	Tan  = (Erad/a)**0.25

	py.plot(r,Tan,color='red',linewidth=3)
	py.plot(r,Trad,'o')
	py.plot(r,Tgas,'--')
	py.legend(['analytic solution','Sedona'])

	py.title('sperhical diffusion')
	py.xlabel('radius (cm)',size=15)
	py.ylabel('temperature (K)',size=15)
	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

		

if __name__=='__main__': compare('')

