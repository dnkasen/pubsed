import pylab as py
import h5py

def compare(pdf):

	py.clf()
	rp = 3e13

	rr,s1,s2 = py.loadtxt('NLTE_analytic_solution.dat',unpack=True,skiprows=1)
	py.plot(rr/rp,s1,color='black',linewidth=4)
	py.plot(rr/rp,s2,color='red',linewidth=4)

	fin = h5py.File('levels_00001.h5','r')

	r   = py.array(fin['x'])
	r   = r/rp
	rho = py.array(fin['rho'])
	nz = len(rho)

	n1 = py.zeros(nz)
	n2 = py.zeros(nz)
	n3 = py.zeros(nz)
	for i in range(nz):
		name = 'zone_' + str(i)
		b = py.array(fin[name + '/Z_1/level_departure'])
		n = py.array(fin[name + '/Z_1/level_fraction'])
		n1[i] = n[0]
		n2[i] = n[1]
		n3[i] = n[2]
	py.plot(r,n1/n2,'o',color='black')
	py.plot(r,n1/n3,'o',color='red')

	py.legend(['n1/n2 analytic','n1/n3 analytic','n1/n2 sedona','n1/n3 sedona'],loc=2)
	py.xlabel('r/r_phot',size=15)
	py.ylabel('level population ratios',size=15)
	py.title('3 level NLTE test')

	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

		

if __name__=='__main__': compare('')