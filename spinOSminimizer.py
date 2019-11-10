"""
Module that performs a non-linear least squares minimization of the spectrescopic and the astrometric data
using the lmfit package.
"""

import lmfit as lm
import numpy as np
import orbit
import time

RV1, RV2, AS = False, False, False


def LMminimizer(guessdict: dict, datadict: dict, search: float):
    """
    Minimizes the provided data to a binary star model, with initial provided guesses and a search radius
    :param guessdict: dictionary containing guesses and 'to-vary' flags for the 11 parameters
    :param datadict: dictionary containing observational data of RV and/or separations
    :param search: float between 0 and 1 that indicated how far a guess will be varied, so the minimizer will maximally
                   vary your guess between (1-search)*guess and (1+search)*guess
    :return: result from the lmfit minimization routine. It is a MinimizerResult object.
    """
    # get guesses and vary flags
    guesses = guessdict['guesses']
    varying = guessdict['varying']

    # setup Parameters object for the solver
    params = lm.Parameters()
    # populate with parameter data
    params.add_many(
        ('e', guesses['e'], varying['e'], (1 - search) * guesses['e'], min(1, (1 + search) * guesses['e'])),
        ('i', guesses['i'], varying['i'], (1 - search) * guesses['i'], min(np.pi, (1 + search) * guesses['i'])),
        ('omega', guesses['omega'], varying['omega'], (1 - search) * guesses['omega'],
         min(2 * np.pi, (1 + search) * guesses['omega'])),
        ('Omega', guesses['Omega'], varying['Omega'], (1 - search) * guesses['Omega'],
         min(2 * np.pi, (1 + search) * guesses['Omega'])),
        ('t0', guesses['t0'], varying['t0']),
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
    # we need to store this on module level so the function to minimize knows quickly which data is included or not
    global RV1, RV2, AS
    try:
        hjds['RV1'] = datadict['RV1'][:, 0]
        data['RV1'] = datadict['RV1'][:, 1]
        errors['RV1'] = datadict['RV1'][:, 2]
        RV1 = True
    except KeyError:
        pass
    try:
        hjds['RV2'] = datadict['RV2'][:, 0]
        data['RV2'] = datadict['RV2'][:, 1]
        errors['RV2'] = datadict['RV2'][:, 2]
        RV2 = True
    except KeyError:
        pass
    try:
        hjds['AS'] = datadict['AS'][:, 0]
        data['east'] = datadict['AS'][:, 1]
        data['north'] = datadict['AS'][:, 2]
        errors['east'], errors['north'] = convert_error_ellipse(datadict['AS'][:, 3], datadict['AS'][:, 4],
                                                                datadict['AS'][:, 5])
        AS = True
    except KeyError:
        pass
    # build a minimizer object
    minimizer = lm.Minimizer(fcn2min, params, fcn_args=(hjds, data, errors), iter_cb=callback)
    print('Starting Minimization...')
    tic = time.time()
    result = minimizer.minimize()
    toc = time.time()
    print('Minimization Complete in {} s!\n'.format(toc - tic))
    return result


def fcn2min(params, hjds, data, errors):
    """
    Define the function to be minimized by the minimizer. It is simply to array of weighted distances from the model to
    the data, schematically:
        fun = array((data[hjd]-model[hjd])/error_on_data(hjd))
    The function will find out which data is omitted.
    :param params: Parameters object from the package lmfit, containing the 11 parameters to fit.
    :param hjds: dictionary of the days of the observations
    :param data: dictionary of the measurements, be it RV or AS data
    :param errors: dictionary of the errors on the measurements
    :return: array with the weighter errors of the data to the model defined by the parameters
    """
    # create the system belonging to the parameters
    system = orbit.System(params.valuesdict())

    if RV1:
        # Get weighted distance for RV1 data
        chisq_rv1 = ((system.primary.radial_velocity_of_hjds(hjds['RV1']) - data['RV1']) / errors['RV1'])
    else:
        # is RV1 not there, make empty list for this part of the data
        chisq_rv1 = np.asarray(list())
    if RV2:
        # Same for RV2
        chisq_rv2 = ((system.secondary.radial_velocity_of_hjds(hjds['RV2']) - data['RV2']) / errors['RV2'])
    else:
        chisq_rv2 = np.asarray(list())
    if AS:
        # same for AS
        chisq_east = ((system.relative.east_of_hjds(hjds['AS']) - data['east']) / errors['east'])
        chisq_north = ((system.relative.north_of_hjds(hjds['AS']) - data['north']) / errors['north'])
    else:
        chisq_east = np.asarray(list())
        chisq_north = np.asarray(list())

    # concatentate the four parts (RV1, RV2, ASeast, ASnorth)
    res = np.concatenate((chisq_rv1, chisq_rv2, chisq_east, chisq_north))
    return res


def callback(params, it, resid, *args):
    pass


def convert_error_ellipse(major, minor, angle):
    """
    Converts error ellipses to actual east and north errors by a sampling the error ellipse in a monte-carlo way and
    then taking the variance in the east and north directions.
    :param major: length of the major axis of the error ellipse
    :param minor: length of the minor axis of the error ellipse
    :param angle: position angle east of north of the major axis
    :return: east and north error
    """
    num = 1000
    east_error = np.zeros(len(major))
    north_error = np.zeros(len(major))
    for i in range(len(east_error)):
        cosa = np.cos(angle[i])
        sina = np.sin(angle[i])
        temp_major = np.random.randn(num) * major[i]
        temp_minor = np.random.randn(num) * minor[i]
        rotated_temp = np.matmul(np.array([[cosa, sina], [-sina, cosa]]), [temp_major, temp_minor])
        east_error[i] = np.std(rotated_temp[0])
        north_error[i] = np.std(rotated_temp[1])
    return east_error, north_error
