import pylab as py

def compare(pdf):
     
    ###########################################
    # compare the output
    ###########################################
    py.clf()
    x,f,c = py.loadtxt('spectrum_1.dat',unpack=1,skiprows=1)
    lam = 3e10/x*1e8
    ff = f/lam**2
    ff = ff/ff.mean()
    py.plot(lam,ff)

    py.xlim(1050,1350)
    py.ylim(0,2)

    x,f = py.loadtxt('analytic_line.dat',unpack=1)
    py.plot(x,f/f.mean(),linewidth=4,color='black')

    py.legend(['SEDONA','analytic Sobolev'])

    py.ion()
    if (pdf != ''): pdf.savefig()
    else:
	    py.show()
	    j = raw_input()
	


if __name__=='__main__': compare('')
