"""
Module that performs a non-linear least squares minimization of the spectrescopic and the astrometric data
using the lmfit package.
This module is developed with lmfit 0.9.14 and numpy 1.17.2, and requires emcee 3.0.0.

Author:
Matthias Fabry, Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
13 Nov 2019
"""

import lmfit as lm
import numpy as np
import binarySystem as bsys
import time

RV1, RV2, AS = False, False, False


def LMminimizer(guessdict: dict, datadict: dict, domcmc: bool):
    """
    Minimizes the provided data to a binary star model, with initial provided guesses and a search radius
    :param doseppaconversion: boolean to indicate whether to convert sep/pa data to east/north
    :param domcmc: boolean to indicate whether to do an MCMC posterior error estimation
    :param guessdict: dictionary containing guesses and 'to-vary' flags for the 11 parameters
    :param datadict: dictionary containing observational data of RV and/or separations
    :return: result from the lmfit minimization routine. It is a MinimizerResult object.
    """
    # get guesses and vary flags
    guesses = guessdict['guesses']
    varying = guessdict['varying']

    # setup Parameters object for the solver
    params = lm.Parameters()
    # populate with parameter data
    params.add_many(
        ('e', guesses['e'], varying['e'], 0, 1),
        ('i', guesses['i'], varying['i']),
        ('omega', guesses['omega'], varying['omega']),
        ('Omega', guesses['Omega'], varying['Omega']),
        ('t0', guesses['t0'], varying['t0']),
        ('k1', guesses['k1'], varying['k1'], 0),
        ('k2', guesses['k2'], varying['k2'], 0),
        ('p', guesses['p'], varying['p'], 0),
        ('gamma1', guesses['gamma1'], varying['gamma1']),
        ('gamma2', guesses['gamma2'], varying['gamma2']),
        ('d', guesses['d'], varying['d'], 0))

    # setup data for the solver
    hjds = dict()
    data = dict()
    errors = dict()
    # we need to store this on module level so the function to minimize knows quickly which data is included or not
    global RV1, RV2, AS
    try:
        hjds['RV1'] = datadict['RV1']['hjds']
        data['RV1'] = datadict['RV1']['RVs']
        errors['RV1'] = datadict['RV1']['errors']
        RV1 = True
    except KeyError:
        pass
    try:
        hjds['RV2'] = datadict['RV2']['hjds']
        data['RV2'] = datadict['RV2']['RVs']
        errors['RV2'] = datadict['RV2']['errors']
        RV2 = True
    except KeyError:
        pass
    try:
        hjds['AS'] = datadict['AS']['hjds']
        data['east'] = datadict['AS']['easts']
        data['north'] = datadict['AS']['norths']
        errors['east'] = datadict['AS']['easterrors']
        errors['north'] = datadict['AS']['northerrors']
        AS = True
    except KeyError:
        pass
    # build a minimizer object
    minimizer = lm.Minimizer(fcn2min, params, fcn_args=(hjds, data, errors))
    print('Starting Minimization...')
    tic = time.time()
    result = minimizer.minimize()
    toc = time.time()
    print('Minimization Complete in {} s!\n'.format(toc - tic))
    if domcmc:
        mcminimizer = lm.Minimizer(fcn2min, params=result.params, fcn_args=(hjds, data, errors))
        print('Starting MCMC sampling using the minimized paramters...')
        tic = time.time()
        newresults = mcminimizer.minimize(method='emcee')
        toc = time.time()
        print('MCMC complete in {} s!\n'.format(toc - tic))
        newresults.params.pretty_print()
        return newresults
    result.params.pretty_print()
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
    system = bsys.System(params.valuesdict())

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


