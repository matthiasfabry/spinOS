"""
Module that performs a non-linear least squares minimization of the spectrescopic and the astrometric data
using the Levenberg=Marquardt algorithm
"""

import scipy.optimize as spopt
import numpy as np
from spinOSloader import SpinOStag
import orbit

n1, n2, nr = 0, 0, 0


def model(hjds, e, i, omega, Omega, t0, k1, k2, p, gamma1, gamma2, d):
    parameters = [e, i, omega, Omega, t0, k1, k2, p, gamma1, gamma2, d]
    relative_orbit, primary_orbit, secondary_orbit = orbit.generate_orbits(parameters)
    primary_rvs = primary_orbit.radial_velocity_of_hjd(hjds[:n1])
    secondary_rvs = secondary_orbit.radial_velocity_of_hjd(hjds[n1:n1 + n2])
    easts = relative_orbit.east_of_hjd(hjds[n1 + n2:n1 + n2 + nr])
    norths = relative_orbit.north_of_hjd(hjds[n1 + n2 + nr:n1 + n2 + 2 * nr])
    return np.concatenate((primary_rvs, secondary_rvs, easts, norths))


def LMminimizer(datadict: dict, tag: SpinOStag):
    bounds = ([0.] * 11, [1., 180., 360., 360., 10000., 2000., 2000., 10000., 2000., 2000., 10000.])
    global n1, n2, nr
    for i in range(11):
        if datadict['flags'][i]:
            bounds[0][i] = datadict['guesses'][i]
            bounds[1][i] = datadict['guesses'][i] + 1e-8
    if tag == SpinOStag.SB2AS:
        n1 = len(datadict['RV1'][:, 0])
        n2 = len(datadict['RV2'][:, 0])
        nr = len(datadict['AS'][:, 0])
        hjds = np.concatenate(
            (datadict['RV1'][:, 0], datadict['RV2'][:, 0], datadict['AS'][:, 0], datadict['AS'][:, 0]))
        measurements = np.concatenate(
            (datadict['RV1'][:, 1], datadict['RV2'][:, 1], datadict['AS'][:, 1], datadict['AS'][:, 2]))
        errors = np.concatenate(
            (datadict['RV1'][:, 2], datadict['RV2'][:, 2], datadict['AS'][:, 3],
             datadict['AS'][:, 4]))  # TODO: modify errors to represent the true north/east errors (not the ellipse)
    elif tag == SpinOStag.SB2:
        n1 = len(datadict['RV1'][:, 0])
        n2 = len(datadict['RV2'][:, 0])
        hjds = np.concatenate((datadict['RV1'][:, 0], datadict['RV2'][:, 0]))
        measurements = np.concatenate((datadict['RV1'][:, 1], datadict['RV2'][:, 1]))
        errors = np.concatenate((datadict['RV1'][:, 2], datadict['RV2'][:, 2]))
    elif tag == SpinOStag.SB1AS:
        n1 = len(datadict['RV1'][:, 0])
        nr = len(datadict['AS'][:, 0])
        hjds = np.concatenate((datadict['RV1'][:, 0], datadict['AS'][:, 0]))
        measurements = np.concatenate((datadict['RV1'][:, 1], datadict['AS'][:, 1]))
        errors = np.concatenate((datadict['RV1'][:, 2], datadict['AS'][:, 2]))
    elif tag == SpinOStag.SB1:
        n1 = len(datadict['RV1'][:, 0])
        hjds = datadict['RV1'][:, 0]
        measurements = datadict['RV1'][:, 1]
        errors = datadict['RV1'][:, 2]
    else:
        nr = len(datadict['AS'][:, 0])
        hjds = datadict['AS'][:, 0]
        measurements = datadict['AS'][:, 1]
        errors = datadict['AS'][:, 2]
    print('Starting Minimization...')
    print(datadict['guesses'])
    bestpars, _ = spopt.curve_fit(model, hjds, measurements, p0=datadict['guesses'], sigma=errors, bounds=bounds)
    print('Minimization Complete!')
    return bestpars
