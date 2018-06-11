#! /usr/bin/env python
# -*- coding:utf-8 -*-

"""Generate SN Ia toy models for Weizmann workshop code-comparison study
(Radiation Transfer and Explosive Thermonuclear Burning in Supernovae,
17-28 June 2018)

The model is defined by its total mass (--mtot) and asymptotic kinetic
energy (--ekin; alternatively it can be determined given the
composition based on Eq. 1 of W07). The density profile can either be
exponential (--densprof expon) or consist of a broken power law with
indices delta,n (--densprof power --densexp <delta>,<n>; see CS89,
K10).

The ejecta is divided into N zones with constant velocity width
(--dvel). The mass of each zone is computed given the zone volume
(radii determined from velocity assuming homologous expansion) and
density profile. Starting from the central zone, we keep adding mass
shells until the ejecta mass reaches 99.99% of the total mass.

The ejecta is supposed to consist of four distinct chemical zones: the
innermost zone consists of stable IGEs (mass set using --mige; 100% Fe
unless --xfracni is set to the relative fraction of stable Ni); then
comes the 56Ni zone (mass at t=0 set using --mni56); then the IME zone
(mass set using --mime; the IMEs to include are specified using --ime
and their relative fraction with --xfracime). Note that some trace
amount of Ti can be included in the 56Ni and IME zones with --xfracti
(we simply replace xfracti of the 56Ni and IME masses with
Ti). Finally, any remaining outermost layer is set to unburnt C/O (the
relative fraction of O is set using --xfraco). The ejecta must contain
some 56Ni and IMEs, but does not necessarily have to include stable
IGEs or unburnt C/O.

       |             ||       ||       ||             |
       | stable IGEs || 56Ni  || IMEs  || unburnt C/O |
       | (optional)  || (+Ti) || (+Ti) || (optional)  |
mass = 0.............................................mtot

The abundance profiles are connected using an analytical function
(--transprof) over a given mass range (--dmige for stable IGE -> 56Ni
connection; --dmni56 for 56Ni -> IME connection; --dmime for IME ->
unburnt C/O connection). Note that one can set dmige = dmni56 = dmime
using the --dmtrans option. The transition profile can either be a
linear function (--transprof linear), an inverse-exponential (aka
'logistic') function with an associated scale factor(--transprof
invexpon --transscl <scale factor>; see M18), or a cosine bell
(--transprof cosine).

The ejecta is evolved to a time (--tend) by solving the first law of
thermodynamics assuming a radiation-dominated gas, local energy
deposition from 56Ni decay, and no diffusion (i.e. the temperature in
each zone is solved independently from adjacent zones). Given these
assumptions, the final temperature can be determined analytically by
noting that the time-weighted internal energy (=t*E(t)) equals the
time-integrated time-weighted decay energy deposition rate
(=Int{t*Q(t) dt}), as noted by K13 (we ignore the time-weighted
internal energy shortly after explosion E(t0)*t0 << Int{Q(t) t dt}). A
minimum temperature can be set using --tempmin.

Last, an output file is generated (--fout) and the density/abundance
profiles are displayed (unless --noplot is set).

Parameters
----------
Typing:

   python mk_snia_toy_model.py -h

will print the usage and input parameters (with their default values))

Examples
--------
1) ejecta with default settings (see python mk_snia_toy_model.py -h):

   python mk_snia_toy_model.py

2) same as 1) but with broken power-law density profile

   python mk_snia_toy_model.py --densprof power --densexp 0,10

3) 1.4 Msun ejecta (default) with Ekin computed based on composition,
   consisting of 0.1 Msun stable IGEs (default), 0.6 Msun 56Ni
   (default), 0.6 Msun IMEs (Mg, Si, S, Ca, all with default relative
   mass fractions), and hence 0.1 Msun unburnt C/O in equal mass
   fractions (default), connected over a mass range 0.1 Msun
   (default) using a cosine bell:

   python mk_snia_toy_model.py --ekinw07 --transprof cosine

4) 1.0 Msun ejecta with Ekin=10^51 erg (default) consisting only of
   56Ni (0.5 Msun) and Si (0.5 Msun), connected over a mass range 0.1
   Msun (default):

   python mk_snia_toy_model.py --mtot 1.0 --mni56 0.5 --mime 0.5 --ime si

References
----------
CS89: Chevalier & Soker (1989), ApJ, 341, 867
J99: Jeffery (1999) arXiv:astro-ph/9907015
K10: Kasen (2010), ApJ, 708, 1025
K13: Katz et al. (2013), arXiv:1301.6766 [astro-ph]
M18: Magee et al. (2018), arXiv:1803.04436v1
W07: Woosley et al. (2007), ApJ, 662, 487

TODO
----
- define grid based on delta_mass as opposed to delta_vel

- adjust delta_vel (increase resolution) in composition transition zones

Revision history
----------------
27 Mar 2018 - first version of code (Stéphane Blondin, SB)

29 Mar 2018 - revised version (Boaz Katz, BK)
              o replaced temperature iteration with analytical calculation 
                (see Katz et al. 2013), and removed references to an initial 
                time t0 (ejecta evolved to final time T_END directly)
              o use a finer grid (in mass coordinates) for abundance profile
                calculations (change_mass_res() function)
              o correction to average density in transition region + special
                treatment of cell containing the break for broken power-law
                density profile
              o added values of various constants to output file
              o added new columns (X_IGE0 (at t=0), X_56Ni0, X_IME, X_CO) to 
                output file and rearranged columns to first display parameters
                that do not depend on the final time

03 Apr 2018 - revised version for testing by workshop participants (SB)
              o code clean-up and added references to radioactive data

05 Apr 2018 - revised version (SB, per Frank Timmes' suggestions)
              o added Python2/3 compatibility
              o removed unused variables for temperature iteration

15 May 2018 - revised version (SB)
              o added option to include some Ti in 56Ni & IME zones (--xfracti)
              o report actual abundances in output file header in addition to requested ones
              o version date stamp
              o rearrange IMEs order in output file by decreasing atomic mass

Author contact
--------------
Stéphane Blondin, stephane.blondin@lam.fr
"""

import sys
import os
import re
import numpy as np

### version number
VERSION = '2018-05-15'

### ensure Python2 (2.6 or 2.7) and Python3 compatibility
if sys.version_info.major == 2:
    input = raw_input # input() to mean raw_input() when running Python2
    
### constants
# (astro)physical constants
AMU = 1.660540e-24           # atomic mass unit (g)
ARAD = 7.5659125e-15         # radiation constant [erg/cm^3/K^4]
MSUN = 1.989e+33             # solar mass (g)
# 56Ni decay
EDECAY_56NI = 1.7206         # energy per 56Ni decay (MeV) - obtained by summing photon energies from http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=56NI&unc=nds
EDECAY_56CO = 3.6072         # energy per 56Co decay (MeV) - obtained by summing photon energies from http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=56CO&unc=nds
MASS_56NI = 55.94212855      # mass of 56Ni nucleus (AMU) - from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Ni&isotype=all
MASS_56CO = 55.93983880      # mass of 56Co nucleus (AMU) - from https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Co&isotype=all
THALF_56NI = 6.075           # 56Ni half-life (days) - from http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=56NI&unc=nds
THALF_56CO = 77.236          # 56Co half-life (days) - from http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=56CO&unc=nds
# conversion factors
DAY2SEC = 86400.0            # days -> sec conversion
MEV2ERG = 1.60217733e-6      # MeV -> erg conversion factor
# misc
EPSILON = 1e-5               # smallish number
MAXFRAC_TI = 1e-4            # maximum value for Ti fraction in 56Ni and IME zones

### defaults
MTOT_INIT = 1.40             # total mass (msun)
EKIN_INIT = 1.00             # asymptotic kinetic energy (1e51 erg)
DVEL_INIT = 100.0            # cell size (km/s)
DENSPROF_INIT = 'expon'      # "density profile: 'expon' (exponential) or 'power' (broken power-law)
DENSEXP_INIT = '0,10'        # exponents for broken power-law density profile: <delta>,<n> e.g. --densexp 0,10
MIGE_INIT = 0.1              # stable IGE mass (msun)
MNI56_INIT = 0.6             # 56Ni mass at t=0 (msun)
MIME_INIT = 0.6              # IME mass (msun)
DMIGE_INIT = 0.1             # mass interval over which stable IGE mass fraction transitions from 1 to 0 (msun)
DMNI56_INIT = 0.1            # mass interval over which 56Ni mass fraction transitions from 1 to 0 (msun)
DMIME_INIT = 0.1             # mass interval over which IME mass fraction transitions from 1 to 0 (msun)
DMFINE_INIT = 1e-4           # resolution of fine grid of masses used for transitions (msun)
TRANSPROF_INIT = 'linear'    # transition profile for mass fraction variation from 1 to 0: 'linear', 'invexpon' (inverse exponential) or 'cosine' (cosine bell)
TRANSSCL_INIT = 1.4e2        # scale factor for 'invexpon' (inverse exponential) transition profile; this default value of 140 ensures X>0.999 at the lower boundary
XIGEFRAC_NI = 0.1            # fraction of stable IGE mass as stable Ni; the rest gets set to stable Fe
XCOFRAC_O = 0.5              # fraction of unburnt C/O mass as O; the rest gets set to C
XFRACTI_INIT = 0.0           # fraction of mass in 56Ni and IME zones set to Ti
T_END = 1.0                  # final time for toy model (days)
TEMP_MIN = 1e3               # minimum allowed temperature (K)
FOUT_INIT = 'snia_toy.dat'   # output file name

### which IMEs to consider
#
# NOTE: can be modified but ensure Sum(XFRACIME_INIT)=1.0
#       (if only one IME is given then --xfracime is set to 1.0 automatically)
#
# in model DDC10 from Alexei Khokhlov:
#
# M(Ca+S+Si+Mg) = 0.466 Msun
# M(Ca) / M(Ca+S+Si+Mg) ~ 0.087
# M(S)  / M(Ca+S+Si+Mg) ~ 0.351
# M(Si) / M(Ca+S+Si+Mg) ~ 0.542
# M(Mg) / M(Ca+S+Si+Mg) ~ 0.020
#
IME_INIT      = 'ca,s,si,mg'              # comma-separated list of IMEs to include
XFRACIME_INIT = '0.087,0.351,0.542,0.020' # comma-separated list of relative IME fractions



###############################################################################

def change_mass_res(dm_oldres, x_oldres, dm_newres):
    
    """for mass grid with cell masses dm_oldres, and abundances
       x_oldres, find abundances at new resolution grid with cell masses 
       dm_newres
    """
    x_newres = dm_newres * 0.0
    l_new = 0
    l_old = 0
    Nnew = len(dm_newres) 
    Nold = len(dm_oldres)
    mold = dm_oldres[l_old]
    mnew = dm_newres[l_new]
    mxaccum = 0.0
    
    while (l_new < Nnew) and (l_old < Nold):
        if mnew <= mold:
            mxaccum += mnew * x_oldres[l_old]
            mold -= mnew
            x_newres[l_new] = mxaccum / dm_newres[l_new]
            mxaccum = 0.0
            l_new += 1
            if l_new < Nnew:
                mnew = dm_newres[l_new] 
        else:
            mxaccum += mold * x_oldres[l_old]
            mnew -= mold
            l_old += 1
            if l_old < Nold:
                mold = dm_oldres[l_old]

    if l_new < Nnew:
        x_newres[l_new] = mxaccum / dm_newres[l_new]
        
    return x_newres   

###############################################################################

if __name__ == '__main__':
    
    import argparse
    import matplotlib
    from matplotlib.pyplot import *

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    #
    # options
    #
    parser.add_argument('--mtot', default=MTOT_INIT, type=float, help='total mass (msun)')
    parser.add_argument('--ekin', default=EKIN_INIT, type=float, help='asymptotic Ekin (1e51 erg)')
    parser.add_argument('--ekinw07', action='store_true', help='compute Ekin based on W07, Eq. 1')
    parser.add_argument('--dvel', default=DVEL_INIT, type=float, help='cell size (km/s)')
    parser.add_argument('--densprof', default=DENSPROF_INIT, type=str, choices=['expon','power'], help="density profile: 'expon' (exponential) or 'power' (broken power-law)")
    parser.add_argument('--densexp', default=DENSEXP_INIT, type=str, help='exponents for broken power-law density profile: <delta>,<n> e.g. --densexp 0,10')    
    parser.add_argument('--tend', default=T_END, type=float, help='final time for toy model (d)')
    parser.add_argument('--tempmin', default=TEMP_MIN, type=float, help='minimum allowed temperature (K)')
    parser.add_argument('--mige', default=MIGE_INIT, type=float, help='stable IGE mass (msun)')
    parser.add_argument('--mni56', default=MNI56_INIT, type=float, help='56Ni mass at t=0 (msun)')
    parser.add_argument('--mime', default=MIME_INIT, type=float, help='IME mass (msun)')
    parser.add_argument('--dmige', default=DMIGE_INIT, type=float, help='mass interval over which stable IGE mass fraction transitions from 1 to 0 (msun)')
    parser.add_argument('--dmni56', default=DMNI56_INIT, type=float, help='mass interval over which 56Ni mass fraction transitions from 1 to 0 (msun)')
    parser.add_argument('--dmime', default=DMIME_INIT, type=float, help='mass interval over which IME mass fraction transitions from 1 to 0 (msun)')
    parser.add_argument('--dmtrans', default=None, type=float, help='to set dmige=dmni56=dmime=dmtrans in one go')
    parser.add_argument('--dmfine', default=DMFINE_INIT, type=float, help='resolution of fine grid of masses for transitions (msun)')
    parser.add_argument('--transprof', default=TRANSPROF_INIT, type=str, choices=['linear', 'invexpon','cosine'], help="transition profile for mass fraction variation from 1 to 0: 'linear', 'invexpon' (inverse exponential) or 'cosine' (cosine bell)")
    parser.add_argument('--transscl', default=TRANSSCL_INIT, type=float, help="scale factor for 'invexpon' (inverse exponential) transition profile")
    parser.add_argument('--xfracni', default=XIGEFRAC_NI, type=float, help='fraction of stable IGE mass as stable Ni; the rest gets set to stable Fe')
    parser.add_argument('--xfraco', default=XCOFRAC_O, type=float, help='fraction of unburnt C/O mass as O; the rest gets set to C')
    parser.add_argument('--xfracti', default=XFRACTI_INIT, type=float, help='fraction of mass in 56Ni and IME zones set to Ti')
    parser.add_argument('--ime', default=IME_INIT, type=str, help='comma-separated list of IMEs to include')    
    parser.add_argument('--xfracime', default=XFRACIME_INIT, type=str, help='comma-separated list of relative IME fractions')    
    parser.add_argument('--minxfrac', default=None, type=float, help='minimum mass fraction for output to file/plot')
    parser.add_argument('--fout', default=FOUT_INIT, type=str, help='output file name')
    parser.add_argument('--noplot', action='store_true', help='disable plotting of density/abundance profiles')
    parser.add_argument('--nowarn', action='store_true', help='disable warning messages')
    parser.add_argument('--debug', action='store_true', help='print various stuff for debugging')
    parser.add_argument('--test', action='store_true', help='for testing purposes')
    args = parser.parse_args()

    print('')
    print('#############################')
    print('       SN Ia toy model'       )
    print('#############################')

   
    #
    # check masses make sense
    #
    mtot = args.mtot
    mige = args.mige
    mni56 = args.mni56
    mime = args.mime
    if (1.0 - (mni56 + mime)/mtot) < EPSILON and mige > EPSILON:
        print('')
        print('WARNING - 56Ni mass + IME mass = total mass; setting IGE mass to 0')
        mige = 0.0
    mburnt = mige + mni56 + mime
    if mburnt > mtot:
        sys.exit("ERROR - burnt mass exceeds total mass! mtot, mburnt = {:.3f}, {:.3f} Msun".format(mtot, mburnt))
    elif mni56 < EPSILON:
        sys.exit("ERROR - 56Ni mass must be > 0! mni56 = {:.3f} Msun".format(mni56))
    elif mime < EPSILON:
        sys.exit("ERROR - IME mass must be > 0! mime = {:.3f} Msun".format(mime))
    else:
        munbco = mtot - mburnt # unburnt mass
        
    #
    # check IMEs
    #
    imes = args.ime.split(',')
    nime = len(imes)
    for ii, ime in enumerate(imes):
        if ime not in IME_INIT:
            sys.exit("ERROR - IME {:s} not in default IME_INIT: {:s}".format(ime, IME_INIT))

    if nime == 1:
        xfracimestr = ['1.0']
        xfracime = [1.0]
    else:
        xfracimestr = args.xfracime.split(',')[:nime]
        xfracime = [float(xx) for xx in xfracimestr]
    xfracimetot = sum(xfracime)
    if np.abs(1.0 - 1.0/xfracimetot) > EPSILON:
        sys.exit("ERROR - relative IME mass fractions don't sum up to 1! sum(xfracime) = {:.5f}".format(xfracimetot))

    #
    # check Ti fraction
    #
    xfracti = args.xfracti
    if (xfracti > MAXFRAC_TI):
        sys.exit("ERROR - xfracti ({:.4e}) cannot exceed MAXFRAC_TI ({:.4e})".format(xfracti, MAXFRAC_TI))
    else:
        mti_ni56 = xfracti * mni56 # Ti mass in 56Ni zone
        mti_ime  = xfracti * mime  # Ti mass in IME zone
        mti = mti_ni56 + mti_ime
        
    print('')
    print('INFO - user-defined ejecta mass and composition:')
    print('')
    print('            Mtot = {:.4e} Msun'.format(mtot))
    print('   M(stable IGE) = {:.4e} Msun of which {:.1f}% Fe and {:.1f}% Ni'.format(mige, (1.0-args.xfracni)*1e2, args.xfracni*1e2))
    print('         M(56Ni) = {:.4e} Msun'.format(mni56))
    sys.stdout.write('          M(IME) = {:.4e} Msun of which'.format(mime))
    for ii, ime in enumerate(imes):
        sys.stdout.write(' {:.1f}% {:s}'.format(xfracime[ii]*1e2, ime.capitalize()))
        if ii == nime-1:
            print('')
        else:
            if ii == nime-2:
                sys.stdout.write(' and')
            else:
                sys.stdout.write(',')
    print('  M(unburnt C/O) = {:.4e} Msun of which {:.1f}% C and {:.1f}% O'.format(munbco, (1.0-args.xfraco)*1e2, args.xfraco*1e2))

    if (xfracti > 0.0):
        print('')
        print('  NOTE: will replace {:.4e} Msun of 56Ni mass and {:.4e} Msun of IME mass with Ti'.format(mti_ni56, mti_ime))
    
    #
    # check mass intervals dmX
    #
    if args.dmtrans is not None:
        dmige = args.dmtrans
        dmni56 = args.dmtrans
        dmime = args.dmtrans
    else:
        dmige = args.dmige
        dmni56 = args.dmni56
        dmime = args.dmime

    # if there are no IGEs or unburnt C/O, set IGE or IME mass intervals to 0
    if mige < EPSILON:
        mige = 0.0
        dmige = 0.0
    if munbco < EPSILON:
        munbco = 0.0
        dmime = 0.0

    # requirements on IGE/56Ni/IME/CO mass given mass intervals
    if mige < 0.5*dmige:
        sys.exit("ERROR - Need to increase IGE mass or decrease dM(IGE) as M(IGE) < dM(IGE)/2! mime, dmige = {:.3f}, {:.3f} Msun".format(mige, dmige))
    if mni56 < 0.5*(dmige+dmni56):
        sys.exit("ERROR - Need to increase 56Ni mass or decrease dM(IGE)+dM(56Ni) as M(56Ni) < [dM(IGE)+dM(56Ni)]/2! mni56, dmige, dmni56 = {:.3f}, {:.3f}, {:.3f} Msun".format(mni56, dmige, dmni56))
    if mime < 0.5*(dmni56+dmime):
        sys.exit("ERROR - Need to increase 56Ni mass or decrease dM(56Ni)+dM(IME) as M(56Ni) < [dM(56Ni)+dM(IME)]/2! mime, dmni56, dmime = {:.3f}, {:.3f}, {:.3f} Msun".format(mime, dmni56, dmime))
    if munbco < 0.5*dmime:
        sys.exit("ERROR - Need to increase unburnt C/O mass or decrease dM(IME) as M(C/O) < dM(IME)/2! munbco, dmime = {:.3f}, {:.3f} Msun".format(munbco, dmime))

    # compute mass coordinate at which mass fraction starts decreasing from 1
    mcoord_ige = mige - 0.5*dmige # IGE mass fraction starts decreasing from 1 at this mass coordinate (unless M(IGE)=0!)
    mcoord_ni56 = mcoord_ige + mni56 + 0.5*(dmige-dmni56) # 56Ni mass fraction starts decreasing from 1 at this mass coordinate
    mcoord_ime = mcoord_ni56 + mime + 0.5*(dmni56-dmime) # IME mass fraction starts decreasing from 1 at this mass coordinate
    if args.debug:
        print('mcoord_ige, mcoord_ni56, mcoord_ime = {:.3f} {:.3f} {:.3f}'.format(mcoord_ige, mcoord_ni56, mcoord_ime))

    #
    # compute Ekin based on W07, Eq. 1 if --ekinw07 is set
    #
    # Ekin = 1.56 M(Ni) + 1.74 M(Fe) + 1.24 M(IME) - Eg + Eint
    #
    # (units=1e51 erg for Ekin, Eg, Eint; Msun for masses)
    #
    # NOTE: Eg and Eint correspond to MCh ejecta, so a warning is
    # issued if the requested total mass differs significantly from MCh
    if args.ekinw07:
        if np.abs(mtot-1.4) > 0.1:
            print('')
            print("WARNING - total mass differs significantly from MCh")
            zzz = input("          ===> apply Eq. 1 of W07 to determine Ekin anyway? [y/n] (default=n): ")
            if zzz == 'y':
                pass
            else:
                sys.exit("ERROR - exiting mk_snia_toy_model.py; adjust mtot or remove --ekinw07 option")

        ebind = 3.35 # gravitational binding energy for MCh WD from W07 (1e51 erg)
        eint = 2.89  # internal energy of MCh WD from W07 (1e51 erg)
        ekin = 1.56 * mni56 + 1.74 * mige + 1.24 * mime - ebind + eint
        print('')
        print('INFO - computed Ekin based on W07 = {:.4e} erg'.format(ekin*1e51))
    else:
        ekin = args.ekin
        print('')
        print('INFO - input Ekin = {:.4e} erg'.format(ekin*1e51))

    #
    # generate density profile at T_END
    #
    # NOTE: dens and vel are zone-centered
    #
    vel = []            # velocity coordinate in km/s
    vel_out_cgs = []    # velocity at zone outer edge in cm/s. This is what sedona expects in the .mod file
    rad = []            # radial coordinate in cm
    dens = []           # density in g/cm^3
    dmass = []          # shell mass in Msun

    # ejecta are evolved to final time T_END (days)
    tend = args.tend
    tend_sec = tend * DAY2SEC
    
    # set innermost shell properties
    dvel = args.dvel    # cell size in km/s
    v0 = 0.0         ; r0 = v0 * tend_sec * 1e5
    r_in = r0 # for sedona
    v1 = v0 + dvel   ; r1 = v1 * tend_sec * 1e5
    vcen = 0.5*(v0+v1)
    rcen = 0.5*(r0+r1)

    if args.densprof == 'expon':

        print('')
        print('INFO - using exponential density profile')

        # compute e-folding velocity for density profile (see J99, line after Eq. A6)
        # ve = sqrt(Ekin / 6Mtot) (units=cgs)
        ve_cgs = np.sqrt(ekin*1e51 / (6*mtot*MSUN))
        ve = ve_cgs * 1e-5 # cm/s -> km/s
        print('       computed e-folding velocity based on J99 = {:.0f} km/s'.format(ve))

        # compute central density at T_END (see J99, Eq. A7)
        # rho_c,0 = Mtot / (8 PI ve^3 t^3) (units=cgs)
        rhoc0 = mtot * MSUN / (8 * np.pi * ve_cgs**3 * tend_sec**3)
        print('       computed central density based on J99 = {:.2e} gcc at {:.0f} d'.format(rhoc0, tend))

        # compute rho @ zone center (rhocen) and mean density over [v0,v1] (rhoave = M/V = Int(rho dV) / V)
        z0 = v0/ve
        z1 = v1/ve
        zcen = 0.5*(z0+z1)
        rhocen = rhoc0 * np.exp(-zcen)
        rhoave = rhoc0 * 3.0 * (np.exp(-z0)*(z0**2+2.0*z0+2.0) - np.exp(-z1)*(z1**2+2.0*z1+2.0)) / (z1**3 - z0**3)
        
    elif args.densprof == 'power':

        densexp = args.densexp.split(',')
        exp_delta, exp_n = int(densexp[0]), int(densexp[1])
        print('')
        print('INFO - using broken power-law density profile with delta, n = {:d}, {:d}'.format(exp_delta, exp_n))
        if exp_delta >= 3 or exp_n <= 3:
            sys.exit("ERROR - we must have delta < 3 and n > 3 for broken power-law density profile! delta, n = {:d}, {:d}".format(exp_delta, exp_n))
        
        # compute transition velocity for broken power-law density profile
        fac3 = (1.0/(3.0-exp_delta) + 1.0/(exp_n-3.0))
        fac5 = (1.0/(5.0-exp_delta) + 1.0/(exp_n-5.0))
        fac = fac3 / fac5
        vt_cgs = np.sqrt(fac*2.0*ekin*1e51 / (mtot*MSUN))
        vt = vt_cgs * 1e-5 # cm/s -> km/s
        print('       computed transition velocity based on K10 = {:.0f} km/s'.format(vt))

        # compute central density at T_END
        rhoc0 = mtot*MSUN / (4 * np.pi * vt_cgs**3 * tend_sec**3) / fac3
        print('       computed central density based on K10 = {:.2e} gcc at {:.0f} d'.format(rhoc0, tend))

        # compute rho @ zone center (rhocen) and mean density over [v0,v1] (rhoave = M/V = Int(rho dV) / V)
        rhocen = rhoc0 * (vcen/vt)**(-exp_delta)
        rhoave = rhoc0 * 3.0 * (v1**(3.0-exp_delta) - v0**(3.0-exp_delta)) / (vt**(-exp_delta) * (3.0-exp_delta)) / (v1**3 - v0**3)

    else:
        sys.exit("ERROR - unknown density profile: {:s}!".format(args.densprof))
        
    if args.debug:
        rhodiff = 1.0 - rhocen/rhoave
        print('rhoave, rhocen, diff = {:.4e} {:.4e} {:.4e}'.format(rhoave, rhocen, rhodiff))

    dvol = 4./3.*np.pi*(r1**3 - r0**3)
    dm = dvol * rhoave / MSUN  # to be consistent with mean density

    vel.append(vcen)    # velocity at zone center
    vel_out_cgs.append(v1 * 1.e5) # velocity at zone outer edge in cm/s. This is what sedona expects in the .mod file
    rad.append(rcen)    # radius at zone center
    dens.append(rhoave) # mean density over [v0,v1]
    dmass.append(dm)    # mass in zone = Int(rho dV)
    
    while (1.0-sum(dmass)/mtot) > 1e-4:

        v0 += dvel       ; r0 = v0 * tend_sec * 1e5
        v1 = v0 + dvel   ; r1 = v1 * tend_sec * 1e5
        vcen = 0.5*(v0+v1)
        rcen = 0.5*(r0+r1)

        # compute rho @ zone center (rhocen) and mean density over [v0,v1] (rhoave = M/V = Int(rho dV) / V)
        if args.densprof == 'expon':
            z0 = v0/ve
            z1 = v1/ve
            zcen = 0.5*(z0+z1)       
            rhocen = rhoc0 * np.exp(-zcen)
            rhoave = rhoc0 * 3.0 * (np.exp(-z0)*(z0**2+2.0*z0+2.0) - np.exp(-z1)*(z1**2+2.0*z1+2.0)) / (z1**3 - z0**3)
        elif args.densprof == 'power':
            if v1 <= vt:
                rhocen = rhoc0 * (vcen/vt)**(-exp_delta)
                rhoave = rhoc0 * 3.0 * (v1**(3.0-exp_delta) - v0**(3.0-exp_delta)) / (vt**(-exp_delta) * (3.0-exp_delta)) / (v1**3 - v0**3)
            elif v0 >= vt:
                rhocen = rhoc0 * (vcen/vt)**(-exp_n)
                rhoave = rhoc0 * 3.0 * (v1**(3.0-exp_n) - v0**(3.0-exp_n)) / (vt**(-exp_n) * (3.0-exp_n)) / (v1**3 - v0**3)
            else:
                # special treatment for cell that contains the break
                if vcen <= vt:
                    rhocen = rhoc0 * (vcen/vt)**(-exp_delta)
                else:
                    rhocen = rhoc0 * (vcen/vt)**(-exp_n)
                numer0 = (vt**(3.0-exp_delta) - v0**(3.0-exp_delta)) / (vt**(-exp_delta) * (3.0-exp_delta))
                numer1 = (v1**(3.0-exp_n) - vt**(3.0-exp_n)) / (vt**(-exp_n) * (3.0-exp_n))
                rhoave = rhoc0 * 3.0 * (numer0 + numer1) / (v1**3 - v0**3)
                
        if args.debug:
            rhodiff = 1.0 - rhocen/rhoave
            print('rhoave, rhocen, diff = {:.4e} {:.4e} {:.4e}'.format(rhoave, rhocen, rhodiff))

        dvol = 4./3.*np.pi*(r1**3 - r0**3)
        dm = dvol * rhoave / MSUN  # to be consistent with mean density

        vel.append(vcen)    # velocity at zone center
        vel_out_cgs.append(v1 * 1.e5) # velocity at zone outer edge in cm/s. This is what sedona expects in the .mod file
        rad.append(rcen)    # radius at zone center
        dens.append(rhoave) # mean density over [v0,v1]
        dmass.append(dm)    # mass in zone = Int(rho dV)
        
    # convert lists to arrays
    vel = np.array(vel)
    rad = np.array(rad)
    dens = np.array(dens)
    dmass = np.array(dmass)
    nd = vel.size # number of zones
    if args.debug:
        print('nd = ',nd)

    # Lagrangian mass coordinate (corresponds to outer zone boundary)
    mass = np.cumsum(dmass)
            
    #
    # set abundances for stable IGEs, 56Ni, IMEs, unburnt C/O
    #
    if dmige+dmni56+dmime > EPSILON:
        print('')
        print('INFO - connecting abundance profiles with {:s} function'.format(args.transprof))
        print('')
        if mige > EPSILON and dmige > EPSILON:
            print('       stable IGE -> 56Ni zone over mass interval [{:.4e},{:.4e}] Msun'.format(mcoord_ige, mcoord_ige+dmige))
        if dmni56 > EPSILON:
            print('       56Ni -> IME zone over mass interval [{:.4e},{:.4e}] Msun'.format(mcoord_ni56, mcoord_ni56+dmni56))
        if munbco > EPSILON and dmime > EPSILON:
            print('       IME -> unburnt C/O zone over mass interval [{:.4e},{:.4e}] Msun'.format(mcoord_ime, mcoord_ime+dmime))

    # first calculate the abundance profiles on a high resolution grid of masses
    dmfine = args.dmfine
    mass_fine = np.arange(dmfine, mass[-1]+dmfine, dmfine)
    N_fine = len(mass_fine)
    dm_fine = np.ones(N_fine)*dmfine         
    xige_fine = np.zeros(N_fine)
    xni56_fine = np.zeros(N_fine)
    xime_fine = np.zeros(N_fine)
    xunbco_fine = np.zeros(N_fine)

    for i in range(N_fine):
        if mass_fine[i] <= mcoord_ige:
            xige_fine[i] = 1.0
        elif mass_fine[i] <= mcoord_ige + dmige:
            if args.transprof == 'linear':
                xige_fine[i] = (mcoord_ige - mass_fine[i]) / dmige + 1.0
            elif args.transprof == 'invexpon':
                xige_fine[i] = 1.0 / (np.exp(args.transscl * (mass_fine[i] - (mcoord_ige + dmige/2.0))) + 1.0)
            elif args.transprof == 'cosine':
                xige_fine[i] = 1.0 - (1.0 - np.cos(np.pi*(mass_fine[i] - mcoord_ige) / dmige)) / 2.0
            xni56_fine[i] = 1.0 - xige_fine[i]
        elif mass_fine[i] < mcoord_ni56:
            xni56_fine[i] = 1.0
        elif mass_fine[i] <= mcoord_ni56 + dmni56:
            if args.transprof == 'linear':
                xni56_fine[i] = (mcoord_ni56 - mass_fine[i]) / dmni56 + 1.0
            elif args.transprof == 'invexpon':
                xni56_fine[i] = 1.0 / (np.exp(args.transscl * (mass_fine[i] - (mcoord_ni56 + dmni56/2.0))) + 1.0)
            elif args.transprof == 'cosine':
                xni56_fine[i] = 1.0 - (1.0 - np.cos(np.pi*(mass_fine[i] - mcoord_ni56) / dmni56)) / 2.0
            xime_fine[i] = 1.0 - xni56_fine[i]
        elif mass_fine[i] <= mcoord_ime:
            xime_fine[i] = 1.0
        elif mass_fine[i] <= mcoord_ime + dmime:
            if args.transprof == 'linear':
                xime_fine[i] = (mcoord_ime - mass_fine[i]) / dmime + 1.0
            elif args.transprof == 'invexpon':
                xime_fine[i] = 1.0 / (np.exp(args.transscl * (mass_fine[i] - (mcoord_ime + dmime/2.0))) + 1.0)
            elif args.transprof == 'cosine':
                xime_fine[i] = 1.0 - (1.0 - np.cos(np.pi*(mass_fine[i] - mcoord_ime) / dmime)) / 2.0
            xunbco_fine[i] = 1.0 - xime_fine[i]
        else:
            xunbco_fine[i] = 1.0
        if args.debug:
            print(mass_fine[i], xige_fine[i], xni56_fine[i], xime_fine[i], xunbco_fine[i])
    
    # Now map the high resolution grid to the actual grid        
    xige = change_mass_res(dm_fine, xige_fine, dmass)
    xni56 = change_mass_res(dm_fine, xni56_fine, dmass)
    xime = change_mass_res(dm_fine, xime_fine, dmass)
    xunbco = change_mass_res(dm_fine, xunbco_fine, dmass)

    # replace part of 56Ni and IME mass with Ti
    xti = (xni56 + xime) * xfracti
    xni56 = xni56 * (1.0 - xfracti)
    xime = xime * (1.0 - xfracti)
    
    print('')
    print('INFO - final ejecta has {:d} zones with Vmax = {:.4e} km/s and'.format(nd, vel.max()))
    print('')
    print('            Mtot = {:.4e} Msun'.format(np.sum(dmass)))
    print('            Ekin = {:.4e} erg'.format(5e9 * np.sum(dmass*MSUN * vel**2))) # 5e9 = 0.5 * 1e10 i.e. 1/2 factor * (km/s->cm/s)^2
    print('   M(stable IGE) = {:.4e} Msun of which {:.1f}% Fe and {:.1f}% Ni'.format(np.sum(dmass*xige), (1.0-args.xfracni)*1e2, args.xfracni*1e2))
    print('     M(56Ni,t=0) = {:.4e} Msun'.format(np.sum(dmass*xni56)))
    sys.stdout.write('          M(IME) = {:.4e} Msun of which'.format(np.sum(dmass*xime)))
    for ii, ime in enumerate(imes):
        sys.stdout.write(' {:.1f}% {:s}'.format(xfracime[ii]*1e2, ime.capitalize()))
        if ii == nime-1:
            print('')
        else:
            if ii == nime-2:
                sys.stdout.write(' and')
            else:
                sys.stdout.write(',')
    print('  M(unburnt C/O) = {:.4e} Msun of which {:.1f}% C and {:.1f}% O'.format(np.sum(dmass*xunbco), (1.0-args.xfraco)*1e2, args.xfraco*1e2))

    if (xfracti > 0.0):
        print('')
        print('     NOTE: M(Ti) = {:.4e} Msun in 56Ni and IME zones'.format(np.sum(dmass*xti)))
            
    #
    # account for 56Ni decay between t~0 and T_END
    #
    decay_const_ni56 = np.log(2) / THALF_56NI / DAY2SEC
    decay_const_co56 = np.log(2) / THALF_56CO / DAY2SEC
    
    t1 = np.exp(-decay_const_ni56 * tend_sec)
    t2 = np.exp(-decay_const_co56 * tend_sec)
    t3 = decay_const_ni56 * (t2-t1) / (decay_const_ni56 - decay_const_co56)

    xni56_old = xni56.copy()
    xni56 = xni56_old * t1
    xco56 = xni56_old * t3          # assumes X(56Co)=0 at t=0
    xfe56 = xni56_old * (1.0-t1-t3) # assumes X(56Co)=X(56Fe from 56Ni decay)=0 at t=0

    print('')
    print('INFO - accounted for 56Ni decay at t = {:.2f} d:'.format(tend))
    print('')
    print('         M(56Ni) = {:.4e} Msun'.format(np.sum(dmass*xni56)))
    print('         M(56Co) = {:.4e} Msun'.format(np.sum(dmass*xco56)))
    print('         M(56Fe) = {:.4e} Msun'.format(np.sum(dmass*xfe56)))

    #
    # set individual IGE abundances
    #
    xni_stable = xige * args.xfracni
    xfe_stable = xige * (1.0 - args.xfracni)
    xni = xni_stable + xni56
    xco = xco56.copy()
    xfe = xfe_stable + xfe56 # xfe56 stands for 56Fe from 56Co decay
    
    #
    # set individual IME abundances (Mg, Si, S, Ca)
    #

    # initialize individual IME mass fractions
    ximeindiv = {} # dictionary containing IME name and associated mass fraction array, e.g. ximeindiv['si']
    for ime in IME_INIT:
        ximeindiv[ime] = np.zeros(nd)

    # set individual IME mass fractions
    for ii, ime in enumerate(imes):
        ximeindiv[ime] = xfracime[ii] * xime

    #
    # set unburnt C/O abundances
    #
    xo = xunbco * args.xfraco
    xc = xunbco * (1.0 - args.xfraco)
            
    #
    # check mass fraction normalization
    #
    xtot = xni + xco + xfe + xo + xc
    for ime in imes:
        xtot += ximeindiv[ime]
        
    for i in range(nd):
        t1 = 1.0 - 1.0/xtot[i]
        if np.abs(t1) > 1e-3:
            if not args.nowarn:
                print('WARNING - Mass fraction not normalized at depth '+str(i)+' : (1 - 1/Xtot) is '+str(t1))
        
    #
    # compute temperate at T_END (days), assuming radiation-dominated gas, no diffusion, local deposition 
    #
    print('')
    print('INFO - computing final temperature at t = {:.2f} d'.format(tend))

    tauni = 1.0 / decay_const_ni56
    tauco = 1.0 / decay_const_co56
    fco = decay_const_ni56 / (decay_const_ni56 - decay_const_co56)
    
    expni = np.exp(-decay_const_ni56 * tend_sec)
    expco = np.exp(-decay_const_co56 * tend_sec)
    
    # integration of exp(-t/tauni)*t from 0 to tend: 
    inttexpni = tauni * (tauni-expni * (tauni+tend_sec))
    
    # integration of exp(-t/tauco)*t from 0 to tend:
    inttexpco = tauco * (tauco-expco * (tauco+tend_sec))
    
    # time-weighted integral of deposition from Ni (per Ni nucleous at t=0)
    qtdtni = inttexpni / tauni * EDECAY_56NI * MEV2ERG
    
    # time-weighted integral of deposition from Co (per Ni nucleous at t=0)
    qtdtco = fco / tauco * (inttexpco-inttexpni) * EDECAY_56CO * MEV2ERG
    
    # total
    qtdt = qtdtni + qtdtco
    
    # time-weighted integral of deposition (per Ni mass at t=0)
    qtdt_per_mass = qtdt / (MASS_56NI*AMU)

    print('')
    print('  computed time-weighted integral of decay energy:')
    print('')
    print('  qtdt_per_nucleous = {:.4e} erg s'.format(qtdt))
    print('  qtdt_per_mass     = {:.4e} erg s/g'.format(qtdt_per_mass))
    
    # calculate temperature analytically from t*E(t) = Int{t*Q(t) dt} (see K13)
    # NOTE: we ignore initial internal energy, t0*E(t0) << Int{Q(t) t dt}
    temp = ( qtdt_per_mass * dens * xni56_old / tend_sec / ARAD )**(0.25)
    
    # set minimal temperature
    temp = (temp >= args.tempmin)*temp + (args.tempmin > temp)*args.tempmin
    
    print('')
    sys.stdout.write('  ===> maximum temperature is {:.4e} K'.format(temp.max()))
    if temp.argmax().size == 1:
        print(' at {:.4e} km/s'.format(vel[temp.argmax()]))
    else:
        print('')

    # display final abundances
    if args.minxfrac is not None:
        print('')
        print('INFO - will set mass fractions > {:.4e}'.format(args.minxfrac))
        ### IGEs
        xfe[np.where(xfe < args.xmin)] = args.xmin
        xni[np.where(xni < args.xmin)] = args.xmin
        xco[np.where(xco < args.xmin)] = args.xmin
        ### IMEs
        for ime in imes:
            ximetmp = ximeindiv[ime]
            ximetmp[np.where(ximetmp < args.xmin)] = args.xmin
            ximeindiv[ime] = ximetmp.copy()
        ### unbunrt C/O
        xo[np.where(xo < args.xmin)] = args.xmin
        xc[np.where(xc < args.xmin)] = args.xmin

    print('')
    print('INFO - final elemental abundances at t = {:.2f} d'.format(tend))
    print('')
    print('           M(Ni) = {:.4e} Msun'.format(np.sum(dmass*xni)))
    print('           M(Co) = {:.4e} Msun'.format(np.sum(dmass*xco)))
    print('           M(Fe) = {:.4e} Msun'.format(np.sum(dmass*xfe)))
    if (xfracti > 0.0):
        print('           M(Ti) = {:.4e} Msun in 56Ni and IME zones'.format(np.sum(dmass*xti)))
    for ime in imes:
        print('           M({:2s}) = {:.4e} Msun'.format(ime.capitalize(), np.sum(dmass*ximeindiv[ime])))
    print('           M(C ) = {:.4e} Msun'.format(np.sum(dmass*xc)))
    print('           M(O ) = {:.4e} Msun'.format(np.sum(dmass*xo)))


    # Prepare for output to sedona
    
    n_elems = len(imes) + 7 # IMEs plus Ni56, Ni58, Co, Fe, Ti, C, O
    ime_elem_ids = {'ca':'20.40', 's':'16.32', 'si':'14.28', 'mg':'12.24'}

    with open(args.fout, 'w') as f:
    
    #
    # output model to file
    #
        
        ### header
        f.write('1D_sphere SNR\n')
        f.write('{:d} {:.4e} {:.4e} {:d}\n'.format(nd, r_in, tend_sec, n_elems))
        ### element list
        f.write('28.56 28.58 27.56 26.56 22.44 ') #IGEs plus titanium. I *think* Ti should be 44, but not sure
        for ime in imes:
            f.write(ime_elem_ids[ime] + ' ')
        f.write('8.16 6.12')
        f.write('\n')
        for i in range(nd):
            f.write('{:.4e} {:.4e} {:.4e}'.format(vel_out_cgs[i], dens[i], temp[i])) 
            f.write(' {:.4e} {:.4e} {:.4e} {:.4e} {:.4e}'.format(xni56[i], xni[i], xco[i], xfe[i], xti[i]))
            for ime in imes:
                f.write(' {:.4e}'.format(ximeindiv[ime][i]))
            f.write(' {:.4e} {:.4e}'.format(xo[i], xc[i]))
            f.write('\n')


    if args.debug:
        os.system('cat '+args.fout)

    #
    # plot
    #
    if not args.noplot:
        print('')
        print('INFO - on to plot')
        #        plt.ion()
        fig = figure(figsize=(10, 8))
        fig.subplots_adjust(left=.1, bottom=.1, right=.95, top=.95, hspace=.3)
        
        # density profile
        ax = fig.add_subplot(321)
        ax.set_xlabel('Velocity [10$^3$ km s$^{-1}$]')
        ax.set_ylabel('Log$_{10}$ Density [g cm$^{-3}$]')
        ax.set_xlim(0., vel.max()/1e3)
        ax.plot(vel/1e3, np.log10(dens), marker='.', label='$t = {:.2f}$ day'.format(tend))
        ax.legend()
        
        ax = fig.add_subplot(322)
        ax.set_xlabel('Mass [M$_\odot$]')
        ax.set_ylabel('Log$_{10}$ Density [g cm$^{-3}$]')
        ax.set_xlim(0., mtot)
        ax.plot(mass, np.log10(dens), marker='.', label='$t = {:.2f}$ day'.format(tend))

        # temperature profile
        ax = fig.add_subplot(323)
        ax.set_xlabel('Velocity [10$^3$ km s$^{-1}$]')
        ax.set_ylabel('Log$_{10}$ Temperature [K]')
        ax.set_xlim(0., vel.max()/1e3)
        ax.plot(vel/1e3, np.log10(temp), marker='.', label='$t = {:.2f}$ day'.format(tend))
        ax.legend()

        ax = fig.add_subplot(324)
        ax.set_xlabel('Mass [M$_\odot$]')
        ax.set_ylabel('Log$_{10}$ Temperature [K]')
        ax.set_xlim(0., mtot)
        ax.plot(mass, np.log10(temp), marker='.', label='$t = {:.2f}$ day'.format(tend))
            
        # abundance profiles - show grid points to check resolution
        ax = fig.add_subplot(325)
        ax.set_xlabel('Velocity [10$^3$ km s$^{-1}$]')
        ax.set_ylabel('Mass Fraction $X_i$')
        ax.set_xlim(0., vel.max()/1e3)
        ax.set_ylim(-.05,1.05)
        ax.plot(vel/1e3, xige, marker='.', label='stable IGE')
        ax.plot(vel/1e3, xni56, marker='.', label='$^{56}$Ni')
        ax.plot(vel/1e3, xime, marker='.', label='IME')
        ax.plot(vel/1e3, xunbco, marker='.', label='unburnt C/O')
        ax.legend(fontsize='small', loc='center right')
        
        ax = fig.add_subplot(326)
        ax.set_xlabel('Mass [M$_\odot$]')
        ax.set_ylabel('Mass Fraction $X_i$')
        ax.set_xlim(0., mtot)
        ax.set_ylim(-.05,1.05)
        ax.plot(mass, xige, marker='.', label='stable IGE')
        ax.plot(mass, xni56, marker='.', label='$^{56}$Ni')
        ax.plot(mass, xime, marker='.', label='IME')
        ax.plot(mass, xunbco, marker='.', label='unburnt C/O')
        
        show()
        #print('')
        #zzz = input("===> Hit <return> to quit <===")
        #ax.clear()
        close()

    print('')
    print('############################' + '#' * len(args.fout))
    print(' THE END - see output file ' + args.fout)
    print('############################' + '#' * len(args.fout))
    print('')
