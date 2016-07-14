import pylab as py
import h5py



def compare(pdf):

	py.clf()

	# read initial energy
	r,v,d,tr,comp = py.loadtxt('constant.mod',unpack=True,skiprows=3)
	vel = v[0]
	rho = d[0]
	esum = 0.0
	n = len(r)
	for i in range(1,n-1):
		esum += tr[i]**4*(r[i]-r[i-1])
	print esum 

	r0 = 1e5
	print r0

	for i in py.arange(5,21,5):

		fname = 'grid_000' + str(i) + '.h5'
		if (i < 10): fname = 'grid_0000' + str(i) + '.h5'
		fin = h5py.File(fname,'r')		
		r = py.array(fin['x'])
		T = py.array(fin['T_rad'])

		#r = r - 0.5*(r[1] - r[0])
		py.plot(r-r0,T,'o')

		time = (py.array(fin['time']))[0]
		D = 3e10/rho*4.0/3.0
		T0 = (esum/(D*time*3.141519)**0.5)**0.25
		rc = r0 + time*vel
		Tan = T0*py.exp(-1.0*(r-rc)**2/(4.0*D*time))

		py.plot(r - r0,Tan)
		fin.close()

	py.title('Advected Pulse')
	py.legend(['Sedona','analytic diffusion'])

	py.xlabel('radius (cm)',size=15)
	py.ylabel('temperature (K)',size=15)



	if (pdf != ''): pdf.savefig()
	else:
		py.ion()
		py.show()
		j = raw_input()

		

if __name__=='__main__': compare('')

