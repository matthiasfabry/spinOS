"""
presenting spinOS: the SPectroscopic and INterferometric Orbital Solution finder.

Goal:
    spinOS computes the best fit orbital solution given:
     1) radial and systemic velocity data for either or both components of a spectroscopic binary and/or
     2) astrometric data containing the separtions and error ellipses, and
     3) an initial guess of the parameters:
        - e        the eccentricity
        - i        the inclination (deg)
        - omega    the longitude of the periastron with respect to the ascending node of the secondary (deg)
        - Omega    the longitude of the ascending node of the seconday measured east of north (deg)
        - t0       the time of periastron passage (day)
        - k1       the semiamplitude of the radial velocity curve of the primary (km/s)
        - k2       the semiamplitude of the radial velocity curve of the secondary (km/s)
        - p        the period of the binary (day)
        - gamma1   the peculiar velocity of the primary (km/s)
        - gamma2   the peculiar velocity of the secondary (km/s)
        - d        the distance (pc)
        For each parameter, a tag True/False should be supplied to decide for the minimizer to vary this parameter. (So
        False means it will keep it fixed at the supplied value)

This program minimizes a binary orbit model to your supplied data and afterwards plots the data and the minimized model.
The program then gives a best fit value for the parameters itemized above, as well as the component masses.
Errors can be calculated using a Markov Chain Monte Carlo (MCMC) method, the reported errors are half of the difference
between the 15.87 and 84.13 percentiles found in the MCMC sampling.

Usage:
To use spinOS, simply run:
 $ python3 spinOS.py pointerfile.txt [plotonly] [seppa] [domcmc] [steps]
where pointerfile.txt is a plain text file with a list to your datafiles and guessfile and:
plotonly (default = false): boolean to indicate only to plot the data with the model created with the guesses
seppa (default = True): boolean to indicate whether your astrometric data is in separation/position angle or not
            (if you have east/north data, put this to false)
domcmc (default = False): boolean to calculate an MCMC error estimation
steps (default = 1000): integer counting the number of MCMC steps

A pointerfile looks like this:
    RV1file WRstarvels.txt
    #RV2file Ostarvels.txt
    ASfile relative_astrometry.txt
    guessfile guesses.txt
the keys on the first column should be 'RV1file', 'RV2file', 'ASfile' and 'guessfile', each with the relative paths to
the respective data. guessfile is mandatory for spinOS to output something. You can exclude data (if unwanted or
unavailable) by commenting (#) lines or omitting them from the list.

spinOS expects the data to be in the following format: All data files should be plain text files, with:
for RV data:
 JD(days) RV(km/s) error_on_RV(km/s)
eg:
 45000 25.1 2.1
 45860 -4.2 1.1
 etc...

for AS data:
 JD(days) E_separation(mas) N_separation(mas) semimajor_ax_errorellipse(mas)
                                                            semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 -2.5 2.4 0.1 0.8 60
 48050 2.1 8.4 0.4 0.5 90
 etc...

Dependencies:
    python 3.7
    numpy 1.17.2
    scipy 1.3.1
    lmfit 0.9.14
    matplotlib 3.1.1
    emcee 3.0.0 (if MCMC error calculation is performed)

Author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
    21 Nov 2019

Version:
    1.5

Acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.

"""
import sys

import matplotlib.pyplot as plt

import binarySystem
import spinOSloader as spl
import spinOSminimizer as spm
import spinOSplotter as spp

# os.system('clear')
print('Hello, this is spinOS, your personal orbital solution finder. I will start promptly!\n')
# read in files


try:
    plotonly = sys.argv[2] == 'True'
except IndexError:
    plotonly = False

try:
    seppa = not (sys.argv[4] == 'False')
except IndexError:
    seppa = True

try:
    domcmc = sys.argv[3] == 'True'
except IndexError:
    domcmc = False

try:
    steps = int(sys.argv[4])
except (IndexError, ValueError):
    steps = 1000

wd, guessdict, datadict = spl.spinOSparser(sys.argv[1], seppa)
# compute best elements
if plotonly:
    bestpars = guessdict['guesses']
    redchisq = 0.
    dof = 0
    minimizationresult = None
else:
    minimizationresult = spm.LMminimizer(guessdict, datadict, domcmc, steps=steps)
    bestpars = minimizationresult.params.valuesdict()
    redchisq = minimizationresult.redchi
    dof = minimizationresult.nfree

# compute model of these elements
system = binarySystem.System(bestpars)

# plot resulting RV curve and resulting apparent orbit
fig1, fig2, rvax, asax = spp.make_plots()
spp.plot_relative_orbit(asax, system)
spp.plot_rv_curves(rvax, system)
spp.plot_rv_data(rvax, datadict, system)
spp.plot_as_data(asax, datadict)

# calculate the resulting masses
primary_mass = system.primary_mass()
secondary_mass = system.secondary_mass()
print('I have come to an optimal solution! These are:')
print('P = {} days'.format(system.p))
print('e = {}'.format(system.e))
print('i = {} (deg)'.format(system.i))
print('omega = {} (deg)'.format(system.secondary.omega))
print('Omega = {} (deg)'.format(system.Omega))
print('K1 = {} (km/s)'.format(system.primary.k))
print('K2 = {} (km/s)'.format(system.secondary.k))
print('t0 = {} (jd)'.format(system.t0))
print('gamma1 = {} (km/s)'.format(system.primary.gamma))
print('gamma2 = {} (km/s)'.format(system.secondary.gamma))
print('d = {} (pc)'.format(system.d))
print('M1 = {} (Msun)'.format(primary_mass))
print('M2 = {} (Msun)\n'.format(secondary_mass))
if not plotonly:
    print('The minimization algorithm stopped on a reduced chi squared of {}'.format(redchisq))
    print('with {} degrees of freedom'.format(dof))
if domcmc:
    figc = spp.plot_corner_diagram(minimizationresult)
    figc.tight_layout()
fig1.tight_layout()
fig2.tight_layout()
plt.show()
print('\nThis was spinOS, thanks for letting me help you!\n')
