import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import h5py



def run_test(pdf="",runcommand=""):

    ###########################################
    # clean up old results and run the code
    ###########################################
    if (runcommand != ""):
        os.system("rm *_spectrum_* plt_* integrated_quantities.dat")
        os.system(runcommand)

    ###########################################
    # compare the output
    ###########################################

    plt.clf()
    failure = 0

    ### load and plot benchmark results from the runs that use a uniform grid

    # get data
    f = h5py.File('../2D_benchmark/optical_spectrum_final.h5','r')
    t_centers = np.array(f['time'])
    t_centers = t_centers/(3600.0*24.0)
    Lnu_mc = np.array(f['Lnu'])
    mu_centers  = np.array(f['mu'])

    n_times = np.size(t_centers)
    n_mu = np.size(mu_centers)

    t_min = 0.2
    t_max = 5
    index_t_min = np.searchsorted(t_centers,t_min)
    index_t_max = np.searchsorted(t_centers,t_max)

    # plot angle-integrated light curve
    L_tot_mc = np.zeros(n_times)
    for i in range(n_mu):
        L_tot_mc += Lnu_mc[:,0,i]
    L_tot_mc = L_tot_mc/(1.0*n_mu)
    plt.scatter(t_centers[index_t_min:index_t_max], L_tot_mc[index_t_min:index_t_max], color='green', s=1)
    plt.plot(t_centers[index_t_min:index_t_max], L_tot_mc[index_t_min:index_t_max], color='green', linewidth=1)

    # plot angle-dependent light curves
    for i in range(n_mu):
        plt.scatter(t_centers[index_t_min:index_t_max], Lnu_mc[index_t_min:index_t_max,0,i], color='blue', s=1)
        plt.plot(t_centers[index_t_min:index_t_max], Lnu_mc[index_t_min:index_t_max,0,i], color='blue', linewidth=1)

    ### load and plot results from the runs that use a nonuniform grid

    # get data
    f = h5py.File('optical_spectrum_final.h5','r')
    t_centers = np.array(f['time'])
    t_centers = t_centers/(3600.0*24.0)
    Lnu = np.array(f['Lnu'])
    mu_centers  = np.array(f['mu'])

    n_times = np.size(t_centers)
    n_mu = np.size(mu_centers)

    # plot angle-integrated light curve
    L_tot = np.zeros(n_times)
    for i in range(n_mu):
        L_tot += Lnu[:,0,i]
    L_tot = L_tot/(1.0*n_mu)
    plt.scatter(t_centers[index_t_min:index_t_max], L_tot[index_t_min:index_t_max], color='black', s=1)
    plt.plot(t_centers[index_t_min:index_t_max], L_tot[index_t_min:index_t_max], color='black', linewidth=1)

    # plot angle-dependent light curves
    for i in range(n_mu):
        plt.scatter(t_centers[index_t_min:index_t_max], Lnu[index_t_min:index_t_max,0,i], color='red', s=1)
        plt.plot(t_centers[index_t_min:index_t_max], Lnu[index_t_min:index_t_max,0,i], color='red', linewidth=1)
        # calculate error
        use = ((t_centers > t_min)*(t_centers < t_max))
        max_err,mean_err = get_error(Lnu[:,0,i],Lnu_mc[:,0,i],x=t_centers,x_comp=t_centers,use = use)
        if (max_err > 0.4): failure = 1
        if (mean_err > 0.1): failure = 1

    use = ((t_centers > t_min)*(t_centers < t_max))
    max_err,mean_err = get_error(L_tot_mc,L_tot,x=t_centers,x_comp=t_centers,use = use)
    if (max_err > 0.25): failure = 2
    if (mean_err > 0.1): failure = 2

    ## make plot
    plt.title('Ellipsoidal Outflow 2D - Nonuniform grid')
    lines = [Line2D([0,1], [0,1], color=col, linewidth=1) for col in ['green','blue','black','red']]
    plt.legend(lines, ['LC angle-integrated, uniform grid','LC angle-dependent, uniform grid','LC angle-integrated, nonuniform grid','LC angle-dependent, nonuniform grid'])
    plt.yscale('log')
    plt.xlim(0,5)
    plt.ylim(6e41,2e43)
    plt.xlabel('Time since explosion [days]',size=12)
    plt.ylabel('Luminosity [erg/s]',size=12)
    if (pdf != ''): pdf.savefig()
    else:
        plt.ion()
        plt.show()
        j = get_input('press any key> ')

    # this should return !=0 if failed
    return failure



#----------------------------------------------
# helper function for calculating the error
# between two numpy arrays
#----------------------------------------------

def get_error(a,b,x=[],x_comp=[],use=[]):

    #-------------------------------------------

    """ Function to calculate the error between two arrays

        Args:
        a: numpy array of result
        b: numpy array of comparison
        use: an array of 0's and 1's telling which element
             in the arrays to include
        x: optional array of x values to go along with a
        x_comp: optional array of x values to go along with b
        (if x and x_comp are set, will interpolate b values to x spacing)

        Returns:
            returns max_error, mean_error in percentages

        Example:
            say you have an array y that is a function of x
            you wnat to see how much it deviates from a reference array y_comp
            but only for values where x > 0.5. Use

            max_error, mean_error = get_error(y,y_comp,use=(x > 0.5))

    """
    #-------------------------------------------------

    # result array
    y = a
    # compare array
    y_comp = b

    # interpolate comparison if wanted
    if (len(x) != 0 and len(x_comp !=0)):
        y_comp = np.interp(x,x_comp,y_comp)

    # cut the array length if wanted
    if (len(use) > 0):
        y = y[use]
        y_comp = y_comp[use]
    err = abs(y - y_comp)

    max_err = max(err/y_comp)
    mean_err = np.mean(err)/np.mean(y_comp)

    return max_err,mean_err

#-----------------------------------------
# little function to just plot up and
# compare results. Assumes code has
# already been run and output files
# are present
#----------------------------------------
if __name__=='__main__':

    # Default to Python 3's input()
    get_input = input
    # If this is Python 2, use raw_input()
    if sys.version_info[:2] <= (2, 7):
        get_input = raw_input

    status = run_test('')
    if (status == 0):
        print ('SUCCESS')
    else:
        print ('FAILURE, code = ' + str(status))