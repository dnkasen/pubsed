import numpy as np
import pylab as py
from matplotlib.backends.backend_pdf import PdfPages

outpdf = PdfPages('foo.pdf')


py.ion()

temp_0 = 1.0e4
texp_0 = 20.0   # days

t = []
r = []

for i in range(1,99):
    if (i < 10): name = "ray_0000" + str(i)
    else: name = "ray_000" + str(i)
    data = np.loadtxt(name)
    x = data[:,0]/3600.0/24.0
    y = data[:,4]
    texp = max(x)/1e9

    t.append(texp)
    r.append(np.average(y[0:40]))


py.plot(t,r,'o',markersize=8,markerfacecolor='none')
x = np.arange(20,70,1)
py.plot(x,temp_0*(texp_0/x),linewidth=2)
py.legend(['sedona','analytic result'])

py.xticks(size=15)
py.yticks(size=15)
py.title('homologous expansion test: mean radiation temp')
py.ylabel('mean interior temperature',size=15)
py.xlabel('time (days)',size=15)
outpdf.savefig()
py.show()
j = raw_input()

outpdf.close()

for i in range(1,99,10):
    if (i < 10): name = "ray_0000" + str(i)
    else: name = "ray_000" + str(i)
    data = np.loadtxt(name)
    x = data[:,0]
    trad = data[:,4]
    tgas = data[:,3]
    texp = max(x)/1e9/3600.0/24.0

    py.plot(x,tgas*(texp_0/texp),linewidth=2,color='black')
    py.plot(x,trad)

py.xlim(0,4e15)
py.title('homologous expansion test: T_rad profile')
py.ylabel('temperature')
py.xlabel('radius (cm)')
outpdf.savefig()
py.show()
j = raw_input()




