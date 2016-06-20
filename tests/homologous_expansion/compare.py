import numpy as np
import pylab as py


def compare(pdf):
  
	py.clf()
	py.ion()

	cols = ['black','blue','red','green','orange']

	texp = 20.5
 	temp_0 = 1.0e4
	vmax = 0.5e9
	r0   = texp*vmax*3600.0*24.
	x,y = np.loadtxt('ray_00001.dat',unpack=1,usecols=[0,4])
	py.plot(x/r0,y,color=cols[0])
	x,y = np.loadtxt('ray_00010.dat',unpack=1,usecols=[0,4])
	py.plot(x/r0,y,color=cols[1])
	x,y = np.loadtxt('ray_00020.dat',unpack=1,usecols=[0,4])
	py.plot(x/r0,y,color=cols[2])
	x,y = np.loadtxt('ray_00030.dat',unpack=1,usecols=[0,4])
	py.plot(x/r0,y,color=cols[3])
	

	x = np.arange(0,2,0.01)*r0
	y = np.zeros(len(x))
	for i in range(len(x)):
		if (x[i] < r0): y[i] = temp_0

	for i in range(0,4):
		xn = x*(texp + (i)*10.0 + 0.5)/texp
		yn = y*texp/(texp + (i)*10.0 + 0.5)
		py.plot(xn/r0,yn,'--',color=cols[i])

	py.xlim(0,3)
	if (pdf != ''): pdf.savefig()
	else:
		py.show()
		j = raw_input('press>')


	py.clf()
	t = []
	r = []

	for i in range(1,62):
		if (i < 10): name = "ray_0000" + str(i)
		else: name = "ray_000" + str(i)
		name = name + '.dat'
		x,y = np.loadtxt(name,unpack=1,usecols=[0,4])
		tn = max(x)/1e9/3600.0/24.0
		t.append(tn)
		r.append(np.average(y[0:40]))


	py.plot(t,r,'o',markersize=8,markerfacecolor='none')
	x = np.arange(20,100,0.1)
	py.plot(x,temp_0*(texp/x),linewidth=2)
	py.legend(['sedona','analytic result'])

	py.xticks(size=15)
	py.yticks(size=15)
	py.title('homologous expansion test: mean radiation temp')
	py.ylabel('mean interior temperature',size=15)
	py.xlabel('time (days)',size=15)

	if (pdf != ''): pdf.savefig()
	else:
		py.show()
		j = raw_input('press>')

if __name__=='__main__': compare('')
