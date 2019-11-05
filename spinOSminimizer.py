"""
Module that performs a non-linear least squares minimization of the spectrescopic and the astrometric data
using the Levenberg=Marquardt algorithm
"""

import lmfit as lm
import scipy.optimize as spopt
import numpy as np
import orbit
import time


def fcn2min(params, hjds, data, errors):
    parvals = params.valuesdict()
    parameters = {'e': parvals['e'], 'i': parvals['i'], 'omega': parvals['omega'], 'Omega': parvals['Omega'],
                  't0': parvals['t0'], 'k1': parvals['k1'], 'k2': parvals['k2'], 'p': parvals['p'],
                  'gamma1': parvals['gamma1'], 'gamma2': parvals['gamma2'], 'd': parvals['d']}
    system = orbit.System(parameters)
    try:
        chisq_rv1 = ((system.primary.radial_velocity_of_hjds(hjds['RV1']) - data['RV1']) / errors['RV1'])
    except KeyError:
        chisq_rv1 = np.asarray(list())
    try:
        chisq_rv2 = ((system.secondary.radial_velocity_of_hjds(hjds['RV2']) - data['RV2']) / errors['RV2'])
    except KeyError:
        chisq_rv2 = np.asarray(list())
    try:
        chisq_east = ((system.relative.east_of_hjds(hjds['AS']) - data['easts']) / errors['easts'])
    except KeyError:
        chisq_east = np.asarray(list())
    try:
        chisq_north = ((system.relative.north_of_hjds(hjds['AS']) - data['norths']) / errors['norths'])
    except KeyError:
        chisq_north = np.asarray(list())
    return np.concatenate((chisq_rv1, chisq_rv2, chisq_east, chisq_north))


def convert_error_ellipses(majors, minors, angles):
    num = 1000
    east_errors = np.zeros(len(majors))
    north_errors = np.zeros(len(majors))
    for i in range(len(east_errors)):
        cosa = np.cos(angles[i])
        sina = np.sin(angles[i])
        temp_majors = np.random.randn(num) * majors[i]
        temp_minors = np.random.randn(num) * minors[i]
        rotated_temp = np.matmul(np.array([[cosa, sina], [-sina, cosa]]), [temp_majors, temp_minors])
        east_errors[i] = np.std(rotated_temp[0])
        north_errors[i] = np.std(rotated_temp[1])
    return east_errors, north_errors


def LMminimizer(guessdict: dict, search, datadict: dict):
    guesses = guessdict['guesses']
    varying = guessdict['varying']
    # setup Parameter objects for the solver
    params = lm.Parameters()
    params.add_many(
        ('e', guesses['e'], varying['e'], (1 - search) * guesses['e'], min(1, (1 + search) * guesses['e'])),
        ('i', guesses['i'], varying['i'], (1 - search) * guesses['i'], min(np.pi, (1 + search) * guesses['i'])),
        ('omega', guesses['omega'], varying['omega'], (1 - search) * guesses['omega'],
         min(2 * np.pi, (1 + search) * guesses['omega'])),
        ('Omega', guesses['Omega'], varying['Omega'], (1 - search) * guesses['Omega'],
         min(2 * np.pi, (1 + search) * guesses['Omega'])),
        ('t0', guesses['t0'], varying['t0'], (1 - search) * guesses['t0'], (1 + search) * guesses['t0']),
        ('k1', guesses['k1'], varying['k1'], (1 - search) * guesses['k1'], (1 + search) * guesses['k1']),
        ('k2', guesses['k2'], varying['k2'], (1 - search) * guesses['k2'], (1 + search) * guesses['k2']),
        ('p', guesses['p'], varying['p'], (1 - search) * guesses['p'], (1 + search) * guesses['p']),
        ('gamma1', guesses['gamma1'], varying['gamma1'], (1 - search) * guesses['gamma1'],
         (1 + search) * guesses['gamma1']),
        ('gamma2', guesses['gamma2'], varying['gamma2'], (1 - search) * guesses['gamma2'],
         (1 + search) * guesses['gamma2']),
        ('d', guesses['d'], varying['d'], (1 - search) * guesses['d'], (1 + search) * guesses['d']))
    # setup data for the solver
    hjds = dict()
    data = dict()
    errors = dict()
    try:
        hjds['RV1'] = datadict['RV1'][:, 0]
        data['RV1'] = datadict['RV1'][:, 1]
        errors['RV1'] = datadict['RV1'][:, 2]
    except KeyError:
        pass
    try:
        hjds['RV2'] = datadict['RV2'][:, 0]
        data['RV2'] = datadict['RV2'][:, 1]
        errors['RV2'] = datadict['RV2'][:, 2]
    except KeyError:
        pass
    try:
        hjds['AS'] = datadict['AS'][:, 0]
        data['east'] = datadict['AS'][:, 1]
        data['north'] = datadict['AS'][:, 2]
        errors['east'], errors['north'] = convert_error_ellipses(datadict['AS'][:, 3], datadict['AS'][:, 4],
                                                                 datadict['AS'][:, 5])
    except KeyError:
        pass
    # build a minimizer object
    minimizer = lm.Minimizer(fcn2min, params, fcn_args=(hjds, data, errors))
    print('Starting Minimization...')
    tic = time.time()
    result = minimizer.minimize()
    toc = time.time()
    print('Minimization Complete in {} s!\n'.format(toc - tic))
    parvals = result.params.valuesdict()
    bestpars = {'e': parvals['e'], 'i': parvals['i'] * 180 / np.pi, 'omega': parvals['omega'] * 180 / np.pi,
                'Omega': parvals['Omega'] * 180 / np.pi,
                't0': parvals['t0'], 'k1': parvals['k1'], 'k2': parvals['k2'], 'p': parvals['p'],
                'gamma1': parvals['gamma1'], 'gamma2': parvals['gamma2'], 'd': parvals['d']}
    return bestpars
