import pylab as py
import constants as pc
import h5py


def compare(pdf):

	py.clf()

	fin = h5py.File('levels_00001.h5','r')
	rp = 3.0e13

	r    = py.array(fin['x'])
	rmid = r - (r[1] - r[0])
	r   = r/rp
	rho = py.array(fin['rho'])
	nz = len(rho)

	n1 = py.zeros(nz)
	n2 = py.zeros(nz)
	for i in range(nz):
		name = 'zone_' + str(i)
		b = py.array(fin[name + '/Z_1/level_departure'])
		n = py.array(fin[name + '/Z_1/level_fraction'])
		n1[i] = n[0]
		n2[i] = n[1]
	py.plot(r,n1/n2,'o',color='black')


	# analytic solution
	T  = 5.0e4
	g1 = 2.0
	g2 = 8.0
	W = 0.5*(1 - (1 - (rp/rmid)**2)**0.5)
	eta = pc.ev_to_ergs*10.2/pc.k/T
	an = g1/g2*py.exp(eta)*(1 + py.exp(-1.0*eta)*(W - 1.0))/W
	py.plot(r,an,linewidth=4,color='red')

	py.ylim(0,110)
	py.ylabel(r'level population ratio: $n_1/n_2$',size=15)
	py.xlabel('radius/r_phot',size=15)
	py.legend(['Sedona result','analytic result'],loc=2)

	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

		

if __name__=='__main__': compare('')

