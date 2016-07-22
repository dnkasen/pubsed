import pylab as py
import h5py


def compare(pdf):

	py.clf()
	m_p = 1.67262158e-24   # mass of proton (g)
	k   = 1.380658e-16     # boltzmann constant (ergs/K)

	# read exact solution and mine
	z,x,rho,e,p,v,cs = py.loadtxt('exact_solution.dat',unpack=True,skiprows=2)
	r_me,rho_me,v_me,T_me = py.loadtxt('ray_00011.dat',unpack=True,usecols=[0,1,2,3])
	p_me = rho_me/m_p*k*T_me
	e_me = p_me/rho_me/(1.4 - 1)
	t = p*m_p/rho/k

	py.suptitle('sedov blastwave hydro test')

	py.subplot(2,2,1)
	py.plot(x,rho,color='red',linewidth=3)
	py.plot(r_me,rho_me,'o',markerfacecolor='none')
	py.ylabel('density')
	py.xlim(0,1.5)

	py.subplot(2,2,2)
	py.plot(x,v,color='red',linewidth=3)
	py.plot(r_me,v_me,'o',markerfacecolor='none')
	py.ylabel('velocity')

	py.subplot(2,2,3)
	py.plot(x,p,color='red',linewidth=3)
	py.plot(r_me,p_me,'o',markerfacecolor='none')
	py.ylabel('pressure')

	py.subplot(2,2,4)
	py.plot(x,e,color='red',linewidth=3)
	py.plot(r_me,e_me,'o',markerfacecolor='none')
	py.ylabel('energy per unit mass ')
	py.yscale('log')
	py.ylim(0,10000.0)


	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

		

if __name__=='__main__': compare('')

