from .SpectrumFile import *

def read_spectrum(name):

    return SpectrumFile(name)


def get_spectrum(name,time=None,mu=None,phi=None):

    s = read_spectrum(name)
    return s.get_spectrum(time=time,mu=mu,phi=phi)

def plot_spectrum(name,time=None,mu=None,phi=None):

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    x,f = get_spectrum(name,time=time,mu=mu,phi=phi)

    if (len(f.shape) < 2):
        nmu = 1
        plt.plot(x,f,color='k',lw=3)

    else:
        nmu = f.shape[1]
        colors = plt.cm.jet(np.linspace(0,1,nmu))
        for i in range(nmu):
            plt.plot(x,f[:,nmu-1-i],color=colors[i],lw=3)

        dmu = 2.0/(1.0*nmu)
        c = np.arange(-1.0 + dmu/2.0, 1.0,dmu)

        cmap = mpl.cm.get_cmap('jet', nmu)
        norm = mpl.colors.BoundaryNorm(c,len(c)) #np.arange(-1.01,1.01,2.0/(1.0*len(c))),len(c)) # len(c)+1)+0.5,1) #len(c))
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # this line may be ommitted for matplotlib >= 3.1
        plt.colorbar(sm,ticks=c)

    plt.ion()
    plt.show()
    j = input()

def plot_lightcurve(name,band=None):

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    s = read_spectrum(name)
    if (band is None):
        t,f = s.get_bolometric_lc()
    else:
        t,f = s.get_band_lc(band)

    if (len(f.shape) < 2):
        nmu = 1
        plt.plot(t,f,color='k',lw=3)

    else:
        nmu = f.shape[1]
        colors = plt.cm.jet(np.linspace(0,1,nmu))
        for i in range(nmu):
            plt.plot(t,f[:,nmu-1-i],color=colors[i],lw=3)

        dmu = 2.0/(1.0*nmu)
        c = np.arange(-1.0 + dmu/2.0, 1.0,dmu)

        cmap = mpl.cm.get_cmap('jet', nmu)
        norm = mpl.colors.BoundaryNorm(c,len(c)) #np.arange(-1.01,1.01,2.0/(1.0*len(c))),len(c)) # len(c)+1)+0.5,1) #len(c))
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # this line may be ommitted for matplotlib >= 3.1
        plt.colorbar(sm,ticks=c)


    plt.ion()
    plt.show()
    j = input()