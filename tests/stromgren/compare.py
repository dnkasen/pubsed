import pylab as py
import h5py


def compare(pdf):

	py.clf()

	# compare problem boxsize = 6.6 kiloparsecs
	kpc = 3.08e21
	rbox = 6.6*kpc  
   
	# plot up final structure
	fin = h5py.File('levels_00010.h5','r')
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
	py.plot(r/rbox,n1,'o',markeredgewidth=2,markerfacecolor='none',color='black')
	py.plot(r/rbox,n2,'--',linewidth=3,color='black')

        # overplot analytic
        ran,xan = py.loadtxt('analytic_solution.dat',unpack=1)
        py.plot(ran/rbox,xan)
        py.plot(ran/rbox,1-xan)

	py.yscale('log')
	py.ylim(1e-5,2)
	py.xlim(0,1.5)
	py.ylabel(r'ion fraction',size=15)
	py.xlabel('r/rbox (rbox = 6.6 kpc)',size=15)
	py.legend(['HI (ground state) fraction','HII fraction'],loc=4)
        py.title('stromgren test: ionization structure')

	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

	# plot up convergence
	py.clf()

	i_arr = py.arange(1,11,1)
	r_arr = py.empty((0))
	for i in i_arr:
                fname = 'levels_0000' + str(i) + '.h5'
                if (i > 9): fname = 'levels_000' + str(i) + '.h5'
                fin = h5py.File(fname,'r')
                r    = py.array(fin['x'])
                for j in range(len(r)-1):
 		        name = 'zone_' + str(j)
 			n = py.array(fin[name + '/Z_1/level_fraction'])
 			if (n[len(n)-1] < 0.5): break
 		fin.close()
 	        r_arr = py.append(r_arr,r[j]/rbox)

        py.plot(i_arr,r_arr)
        py.plot(i_arr,r_arr,'o')
        py.ylabel('ionization front: r/rbox (rbox = 6.6 kpc)',size=15)
        py.xlabel('iteration #',size=15)
        py.title('stromgren test: convergence of ionization front')
        
        py.ylim(0.6,2.5)
 	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()


if __name__=='__main__': compare('')

