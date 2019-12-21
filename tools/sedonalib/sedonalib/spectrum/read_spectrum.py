from SpectrumFile import *

def read_spectrum(name):

    return SpectrumFile(name)


def get_spectrum(name,time=None,mu=None,phi=None):

    s = read_spectrum(name)
    return s.get_spectrum(time=time,mu=mu,phi=phi)

def plot_spectrum(name,time=None,mu=None,phi=None):

    import matplotlib.pyplot as plt
    x,f = get_spectrum(name,time=time,mu=mu,phi=phi)
    plt.plot(x,f)
