"""
presenting spinOS: the SPectroscopic and INterferometric Orbital Solution finder.

author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

date:
    23 Oct 2019

version:
    alpha

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
relative_orbit, primary_orbit, secondary_orbit = orbit.generate_orbits(bestpars)

# plot resulting RV curve and resulting apparent orbit
spinOSplotter.make_plots(relative_orbit, primary_orbit, secondary_orbit, data_dict)

# calculate the resulting masses
primary_mass = np.power(1 - primary_orbit.e ** 2, 1.5) * (
        primary_orbit.k + secondary_orbit.k) ** 2 * secondary_orbit.k * primary_orbit.P * 84600 / (
                       2 * np.pi * c.G * secondary_orbit.sini ** 3)
secondary_mass = np.power(1 - secondary_orbit.e ** 2, 1.5) * (
        primary_orbit.k + secondary_orbit.k) ** 2 * primary_orbit.k * secondary_orbit.P * 84600 / (
                         2 * np.pi * c.G * secondary_orbit.sini ** 3)
print('I have come to an optimal solution! These are:')
print('P = {} days'.format(primary_orbit.P))
print('e = {}'.format(primary_orbit.e))
print('i = {} (deg)'.format(primary_orbit.i_deg))
print('omega = {} (deg)'.format(primary_orbit.omega_deg))
print('Omega = {} (deg)'.format(primary_orbit.Omega_deg))
print('K1 = {} (km/s)'.format(primary_orbit.k))
print('K2 = {} (km/s)'.format(secondary_orbit.k))
print('t0 = {} (hjd mod P)'.format(primary_orbit.t0))
print('gamma1 = {} (km/s)'.format(primary_orbit.gamma))
print('gamma2 = {} (km/s)'.format(secondary_orbit.gamma))
print('d = {} (pc)')
print('M1 = {}'.format(primary_mass / c.m_sun))
print('M2 = {}'.format(secondary_mass / c.m_sun))
plt.show()
print('This was spinOS, thank you for letting me help you')
