"""
presenting spinOS: the SPectroscopic and INterferometric Orbital Solution finder.

author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

date:
    23 Oct 2019

version:
    1.0

goal:
    spinOS computes the best fit orbital solution given:
     1) radial and systemic velocity data for either or both components of a spectroscopic binary and/or
     2) astrometric data containing the separtions and position angles, and
     3) an initial guess of the parameters:
        - e        the eccentricity
        - i        the inclination
        - omega    the longitude of the periastron with respect to the ascending node of the secondary
        - Omega    the longitude of the ascending node of the seconday measured east of north
        - T        the time of periastron passage
        - K1       the semiamplitude of the radial velocity curve of the primary
        - K2       the semiamplitude
        - P        the period of the binary
        - gamma1   the (apparent) systemic velocity of the primary
        - gamma2   the (apparent) systemic velocity of the secondary
        - d        the distance
        For each parameter, a tag True/False should be supplied to decide for the minimizer to keep this parameter fixed

acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.

"""
import sys
import os
import numpy as np
import spinOSloader
import spinOSminimizer
import spinOSplotter
import orbit
import constants as c
import matplotlib.pyplot as plt

# os.system('clear')
print('Hello, this is spinOS, your personal orbital solution finder. I will start promptly!\n')
# read in files

data_dict, tag = spinOSloader.spinOSparser(sys.argv[1])

# compute best elements
if tag == spinOSloader.SpinOStag.PLOTONLY:
    bestpars = data_dict['guesses']
else:
    bestpars = spinOSminimizer.LMminimizer(data_dict, tag)

# compute model of these elements
system = orbit.System(bestpars)

# plot resulting RV curve and resulting apparent orbit
fig, rvax, relax = spinOSplotter.make_plots_horizontal(system)
spinOSplotter.plot_relative_orbit(relax, system)
spinOSplotter.plot_rv_curves(rvax, system)

# calculate the resulting masses
primary_mass = np.power(1 - system.e ** 2, 1.5) * (
        system.primary.k + system.secondary.k) ** 2 * system.secondary.k * system.p * 86400 / (
                       2 * np.pi * c.G * system.sini ** 3)
secondary_mass = np.power(1 - system.e ** 2, 1.5) * (
        system.primary.k + system.secondary.k) ** 2 * system.primary.k * system.p * 86400 / (
                         2 * np.pi * c.G * system.sini ** 3)
print('I have come to an optimal solution! These are:')
print('P = {} days'.format(system.p))
print('e = {}'.format(system.e))
print('i = {} (deg)'.format(system.i*180/np.pi))
print('omega = {} (deg)'.format(system.secondary.omega*180/np.pi))
print('Omega = {} (deg)'.format(system.Omega*180/np.pi))
print('K1 = {} (km/s)'.format(system.primary.k))
print('K2 = {} (km/s)'.format(system.secondary.k))
print('t0 = {} (hjd mod P)'.format(system.t0))
print('gamma1 = {} (km/s)'.format(system.primary.gamma))
print('gamma2 = {} (km/s)'.format(system.secondary.gamma))
print('d = {} (pc)'.format(system.d))
print('M1 = {} (Msun)'.format(primary_mass / c.m_sun))
print('M2 = {} (Msun)'.format(secondary_mass / c.m_sun))
plt.show()
print('This was spinOS, thank you for letting me help you')
