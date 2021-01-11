"""
Module that performs a non-linear least squares minimization of the spectroscopic and/or astrometric data
using the lmfit package.
This module is developed with lmfit 0.9.14 and numpy 1.17.2, and requires emcee 3.0.0.

Author:
Matthias Fabry, Instituut voor Sterrekunde, KU Leuven, Belgium

"""
import time

import lmfit as lm
import numpy as np

from modules.binary_system import System

RV1, RV2, AS = False, False, False
LAS, LRV = 0, 0


def LMminimizer(guess_dict: dict, datadict: dict, domcmc: bool, steps: int = 1000, as_weight: float = None,
                lock_g: bool = None, lock_q: bool = None):
    """
    Minimizes the provided data to a binary star model, with initial provided guesses and a search radius
    :param as_weight: weight to give to the astrometric data, optional.
    :param domcmc: boolean to indicate whether to do an MCMC posterior error estimation
    :param guess_dict: dictionary containing guesses and 'to-vary' flags for the 11 parameters
    :param datadict: dictionary containing observational data of RV and/or separations
    :param steps: integer giving the number of steps the MCMC should perform
    :param lock_g: boolean to indicate whether to lock gamma1 to gamma2
    :param lock_q: boolean to indicate whether to lock k2 to k1/q, and that q is supplied rather than k2 in that field
    :return: result from the lmfit minimization routine. It is a MinimizerResult object.
    """

    # setup data for the solver
    hjds = dict()
    data = dict()
    errors = dict()
    # we need to store this on module level so the function to minimize knows quickly which data is included or not
    global RV1, RV2, AS
    RV1 = RV2 = AS = False
    global LAS, LRV
    LAS = LRV = 0
    try:
        hjds['RV1'] = datadict['RV1']['hjds']
        data['RV1'] = datadict['RV1']['RVs']
        errors['RV1'] = datadict['RV1']['errors']
        RV1 = True
        LRV += len(data['RV1'])
    except KeyError:
        pass
    try:
        hjds['RV2'] = datadict['RV2']['hjds']
        data['RV2'] = datadict['RV2']['RVs']
        errors['RV2'] = datadict['RV2']['errors']
        RV2 = True
        LRV += len(data['RV2'])
    except KeyError:
        pass
    try:
        hjds['AS'] = datadict['AS']['hjds']
        data['east'] = datadict['AS']['easts']
        data['north'] = datadict['AS']['norths']
        errors['east'] = datadict['AS']['easterrors']
        errors['north'] = datadict['AS']['northerrors']
        AS = True
        LAS += 2*len(data['east'])
    except KeyError:
        pass

    # setup Parameters object for the solver
    params = lm.Parameters()
    # populate with parameter data
    params.add_many(
        ('e', guess_dict['e'][0], guess_dict['e'][1], 0, 1 - 1e-5),
        ('i', guess_dict['i'][0], guess_dict['i'][1]),
        ('omega', guess_dict['omega'][0], guess_dict['omega'][1]),
        ('Omega', guess_dict['Omega'][0], guess_dict['Omega'][1]),
        ('t0', guess_dict['t0'][0], guess_dict['t0'][1]),
        ('p', guess_dict['p'][0], guess_dict['p'][1], 0),
        ('mt', guess_dict['mt'][0], guess_dict['mt'][1], 0),
        ('d', guess_dict['d'][0], guess_dict['d'][1], 0),
        ('k1', guess_dict['k1'][0], guess_dict['k1'][1], 0),
        ('gamma1', guess_dict['gamma1'][0], guess_dict['gamma1'][1]),
        ('k2', guess_dict['k2'][0], guess_dict['k2'][1], 0),
        ('gamma2', guess_dict['gamma2'][0], guess_dict['gamma2'][1])
    )

    if lock_g:
        params['gamma2'].set(expr='gamma1')
    if lock_q:
        params.add('q', value=params['k1']/params['k2'], vary=False)
        params['k2'].set(expr='k1/q')

    # put e to a non zero value to avoid conditioning problems in MCMC
    if params['e'].value < 1e-8:
        print('Warning: eccentricity is put to 1e-8 to avoid conditioning issues!')
        params['e'].set(value=1e-8)

    if RV1 and RV2:
        if not AS:
            for key in 'd', 'i', 'Omega', 'mt':
                params[key].set(vary=False)
    elif RV1:
        for key in 'k2', 'gamma2', 'd':
            params[key].set(vary=False)
        if not AS:
            for key in 'i', 'Omega', 'mt':
                params[key].set(vary=False)
    elif AS:
        for key in 'k1', 'gamma1', 'k2', 'gamma2':
            params[key].set(vary=False)
    else:
        raise ValueError('No data supplied! Cannot minimize.')

    # build a minimizer object
    minimizer = lm.Minimizer(fcn2min, params, fcn_args=(hjds, data, errors, as_weight))
    print('Starting Minimization with {}{}{}...'.format('primary RV data, ' if RV1 else '',
                                                        'secondary RV data, ' if RV2 else '',
                                                        'astrometric data' if AS else ''))
    tic = time.time()
    result = minimizer.minimize()
    toc = time.time()
    print('Minimization Complete in {} s!\n'.format(toc - tic))
    lm.report_fit(result.params)
    print('\n')
    rms_rv1, rms_rv2, rms_as = 0, 0, 0
    system = System(result.params.valuesdict())
    if RV1:
        # weigh with number of points for RV1 data
        rms_rv1 = np.sqrt(
            np.sum((system.primary.radial_velocity_of_hjds(hjds['RV1']) - data['RV1']) ** 2) / len(data['RV1']))
    if RV2:
        # Same for RV2
        rms_rv2 = np.sqrt(
            np.sum((system.secondary.radial_velocity_of_hjds(hjds['RV2']) - data['RV2']) ** 2) / len(data['RV2']))
    if AS:
        # same for AS
        omc2E = np.sum((system.relative.east_of_hjds(hjds['AS']) - data['east']) ** 2)
        omc2N = np.sum((system.relative.north_of_hjds(hjds['AS']) - data['north']) ** 2)
        rms_as = np.sqrt(omc2E + omc2N / (len(data['east']) + len(data['north'])))
    if domcmc:
        mcminimizer = lm.Minimizer(fcn2min, params=result.params, fcn_args=(hjds, data, errors))
        print('Starting MCMC sampling using the minimized parameters...')
        tic = time.time()
        newresults = mcminimizer.emcee(steps=steps)
        toc = time.time()
        print('MCMC complete in {} s!\n'.format(toc - tic))
        lm.report_fit(newresults.params)
        print('\n')
        return newresults, rms_rv1, rms_rv2, rms_as
    return result, rms_rv1, rms_rv2, rms_as


def fcn2min(params, hjds, data, errors, weight=None):
    """
    Define the function to be minimized by the minimizer. It is simply to array of weighted distances from the model to
    the data, schematically:
        fun = array((data[hjd]-model[hjd])/error_on_data(hjd))
    The function will find out which data is omitted.
    :param weight: multiplicative weight to give to the astrometric points, optional. If None, no additional weight is
    applied
    :param params: Parameters object from the package lmfit, containing the 11 parameters to fit.
    :param hjds: dictionary of the days of the observations
    :param data: dictionary of the measurements, be it RV or AS data
    :param errors: dictionary of the errors on the measurements
    :return: array with the weighted errors of the data to the model defined by the parameters
    """
    # create the system belonging to the parameters
    system = System(params.valuesdict())

    if RV1:
        # Get weighted distance for RV1 data
        chisq_rv1 = ((system.primary.radial_velocity_of_hjds(hjds['RV1']) - data['RV1']) / errors['RV1'])
        if weight:
            chisq_rv1 *= (1-weight) * (LAS + LRV) / LRV
    else:
        # is RV1 not there, make empty list for this part of the data
        chisq_rv1 = np.asarray(list())
    if RV2:
        # Same for RV2
        chisq_rv2 = ((system.secondary.radial_velocity_of_hjds(hjds['RV2']) - data['RV2']) / errors['RV2'])
        if weight:
            chisq_rv2 *= (1 - weight) * (LAS + LRV) / LRV
    else:
        chisq_rv2 = np.asarray(list())
    if AS:
        # same for AS
        chisq_east = ((system.relative.east_of_hjds(hjds['AS']) - data['east']) / errors['east'])
        chisq_north = ((system.relative.north_of_hjds(hjds['AS']) - data['north']) / errors['north'])
        if weight:
            chisq_east *= weight * (LAS + LRV) / LAS
            chisq_north *= weight * (LAS + LRV) / LAS
    else:
        chisq_east = np.asarray(list())
        chisq_north = np.asarray(list())

    # concatentate the four parts (RV1, RV2, ASeast, ASnorth)
    res = np.concatenate((chisq_rv1, chisq_rv2, chisq_east, chisq_north))
    return res
