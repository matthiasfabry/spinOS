"""
Module that performs a non-linear least squares minimization of the spectrescopic and the astrometric data
using the Levenberg=Marquardt algorithm
"""

import scipy.optimize as spopt
import numpy as np
from spinOSloader import SpinOStag
import orbit
import time

n1, n2, nr = 0, 0, 0


def model(hjds, e, i, omega, Omega, t0, k1, k2, p, gamma1, gamma2, d):
    parameters = {'e': e, 'i': i, 'omega': omega, 'Omega': Omega, 't0': t0, 'k1': k1, 'k2': k2, 'p': p,
                  'gamma1': gamma1, 'gamma2': gamma2, 'd': d}
    system = orbit.System(parameters)
    primary_rvs = system.primary.radial_velocity_of_hjd(hjds[:n1])
    secondary_rvs = system.secondary.radial_velocity_of_hjd(hjds[n1:n1 + n2])
    easts = system.relative.east_of_hjd(hjds[n1 + n2:n1 + n2 + nr])
    norths = system.relative.north_of_hjd(hjds[n1 + n2 + nr:n1 + n2 + 2 * nr])
    return np.concatenate((primary_rvs, secondary_rvs, easts, norths))


def convert_error_ellipses(majors, minors, angles):
    num = 1000
    east_errors = np.zeros(nr)
    north_errors = np.zeros(nr)
    for i in range(nr):
        cosa = np.cos(angles[i])
        sina = np.sin(angles[i])
        temp_majors = np.random.randn(num) * majors[i]
        temp_minors = np.random.randn(num) * minors[i]
        rotated_temp = np.matmul(np.array([[cosa, sina], [-sina, cosa]]), [temp_majors, temp_minors])
        east_errors[i] = np.std(rotated_temp[0])
        north_errors[i] = np.std(rotated_temp[1])
    return east_errors, north_errors


def LMminimizer(datadict: dict, tag: SpinOStag):
    guess = [datadict['guesses']['e'], datadict['guesses']['i'], datadict['guesses']['omega'],
             datadict['guesses']['Omega'], datadict['guesses']['t0'], datadict['guesses']['k1'],
             datadict['guesses']['k2'], datadict['guesses']['p'], datadict['guesses']['gamma1'],
             datadict['guesses']['gamma2'], datadict['guesses']['d']]
    print(guess)
    bounds = ([0., 20*np.pi/180, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [1., np.pi, 2 * np.pi, 2 * np.pi, 10000., 2000., 2000., 10000., 500., 500., 10000.])
    global n1, n2, nr  # we need these globals to keep track of the number of measurements
    if tag == SpinOStag.SB2AS:
        n1 = len(datadict['RV1'][:, 0])
        n2 = len(datadict['RV2'][:, 0])
        nr = len(datadict['AS'][:, 0])
        hjds = np.concatenate(
            (datadict['RV1'][:, 0], datadict['RV2'][:, 0], datadict['AS'][:, 0], datadict['AS'][:, 0]))
        measurements = np.concatenate(
            (datadict['RV1'][:, 1], datadict['RV2'][:, 1], datadict['AS'][:, 1], datadict['AS'][:, 2]))
        east_errors, north_errors = convert_error_ellipses(datadict['AS'][:, 3], datadict['AS'][:, 4],
                                                           datadict['AS'][:, 5])
        errors = np.concatenate((datadict['RV1'][:, 2], datadict['RV2'][:, 2], east_errors, north_errors))
    elif tag == SpinOStag.SB2:
        n1 = len(datadict['RV1'][:, 0])
        n2 = len(datadict['RV2'][:, 0])
        hjds = np.concatenate((datadict['RV1'][:, 0], datadict['RV2'][:, 0]))
        measurements = np.concatenate((datadict['RV1'][:, 1], datadict['RV2'][:, 1]))
        errors = np.concatenate((datadict['RV1'][:, 2], datadict['RV2'][:, 2]))
    elif tag == SpinOStag.SB1AS:
        n1 = len(datadict['RV1'][:, 0])
        nr = len(datadict['AS'][:, 0])
        hjds = np.concatenate((datadict['RV1'][:, 0], datadict['AS'][:, 0], datadict['AS'][:, 0]))
        measurements = np.concatenate((datadict['RV1'][:, 1], datadict['AS'][:, 1], datadict['AS'][:, 2]))
        east_errors, north_errors = convert_error_ellipses(datadict['AS'][:, 3], datadict['AS'][:, 4],
                                                           datadict['AS'][:, 5])
        errors = np.concatenate((datadict['RV1'][:, 2], datadict['AS'][:, 2], east_errors, north_errors))
    elif tag == SpinOStag.SB1:
        n1 = len(datadict['RV1'][:, 0])
        hjds = datadict['RV1'][:, 0]
        measurements = datadict['RV1'][:, 1]
        errors = datadict['RV1'][:, 2]
    else:
        nr = len(datadict['AS'][:, 0])
        hjds = np.concatenate((datadict['AS'][:, 0], datadict['AS'][:, 0]))
        measurements = np.concatenate((datadict['AS'][:, 1], datadict['AS'][:, 2]))
        east_errors, north_errors = convert_error_ellipses(datadict['AS'][:, 3], datadict['AS'][:, 4],
                                                           datadict['AS'][:, 5])
        errors = np.concatenate((east_errors, north_errors))

    print('Starting Minimization...')
    tic = time.time()
    bestpars, _ = spopt.curve_fit(model, hjds, measurements, p0=guess, sigma=errors, bounds=bounds, verbose=2)
    print('Minimization Complete in {}!'.format(time.time() - tic))
    np.savetxt('bestpars.txt', bestpars)
    bestpars = {'e': bestpars[0], 'i': bestpars[1], 'omega': bestpars[2], 'Omega': bestpars[3], 't0': bestpars[4],
                'k1': bestpars[5], 'k2': bestpars[6], 'p': bestpars[7], 'gamma1': bestpars[8], 'gamma2': bestpars[9],
                'd': bestpars[10]}
    return bestpars
