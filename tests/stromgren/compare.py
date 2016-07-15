import pylab as py
import constants as pc
import h5py


def compare(pdf):

	py.clf()

	fin = h5py.File('levels_00030.h5','r')

	r    = py.array(fin['x'])
	rmid = r - (r[1] - r[0])
	rho = py.array(fin['rho'])
	nz = len(rho)

	n1 = py.zeros(nz)
	n2 = py.zeros(nz)
	for i in range(nz):
		name = 'zone_' + str(i)
		b = py.array(fin[name + '/Z_1/level_departure'])
		n = py.array(fin[name + '/Z_1/level_fraction'])
		n1[i] = n[0]
		n2[i] = n[len(n)-1]
	py.plot(r/max(r),n1,color='black')
	py.plot(r/max(r),n2,'-',color='black')
	py.yscale('log')
	py.ylim(1e-5,2)
	
	#py.ylim(0,110)
	#py.ylabel(r'level population ratio: $n_1/n_2$',size=15)
	#py.xlabel('radius/r_phot',size=15)
	#py.legend(['Sedona result','analytic result'],loc=2)

	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

		

if __name__=='__main__': compare('')

