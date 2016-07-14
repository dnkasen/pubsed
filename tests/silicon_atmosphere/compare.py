import pylab as py


def compare(pdf):

    py.clf()

    x,f,c = py.loadtxt('spectrum_1.dat',unpack=1,skiprows=1)
    lam = 3e10/x*1e8
    ff = f/lam**2
    lam = lam[::-1]
    ff  = ff[::-1]
    b = py.interp(7500.0,lam,ff)
    ff = ff/b

    py.plot(lam,ff)

    ls,fs,es = py.loadtxt('synow_spectrum.dat',unpack=1)
    b = py.interp(7500.0,ls,fs)
    fs = fs/b

    py.plot(ls,fs)
    py.xlim(2000,10000)
    py.title('silicon atmosphere')
    py.xlabel('wavelength (angstroms)',size=15)
    py.ylabel('flux',size=15)
    py.legend(['sedona','Synow'])

    if (pdf != ''): pdf.savefig()
    else:
        py.ion()
        py.show()
        j = raw_input()

if __name__=='__main__': compare('')
