import pylab as py
import h5py


def compare(pdf):

	py.clf()
	m_p = 1.67262158e-24   # mass of proton (g)
	k   = 1.380658e-16     # boltzmann constant (ergs/K)

	# read exact solution and mine
	z,x,rho,p,v = py.loadtxt('exact_solution.dat',unpack=True)
	r_me,rho_me,v_me,T_me = py.loadtxt('ray_00005.dat',unpack=True,usecols=[0,1,2,3])
	r_me = r_me - 1e3
	p_me = rho_me/m_p*k*T_me
	t = p*m_p/rho/k

	py.suptitle('sod shock tube hydro test')

	py.subplot(2,2,1)
	py.plot(x,rho,color='red',linewidth=3)
	py.plot(r_me,rho_me,'o',markerfacecolor='none')
	py.ylabel('density')

	py.subplot(2,2,2)
	py.plot(x,v,color='red',linewidth=3)
	py.plot(r_me,v_me,'o',markerfacecolor='none')
	py.ylabel('velocity')

	py.subplot(2,2,3)
	py.plot(x,p,color='red',linewidth=3)
	py.plot(r_me,p_me,'o',markerfacecolor='none')
	py.ylabel('pressure')

	py.subplot(2,2,4)
	py.plot(x,t,color='red',linewidth=3)
	py.plot(r_me,T_me,'o',markerfacecolor='none')
	py.ylabel('temperature')

	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

		

if __name__=='__main__': compare('')

