"""
Copyright 2020-2024 Matthias Fabry
This file is part of spinOS.

spinOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

spinOS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with spinOS.  If not, see <https://www.gnu.org/licenses/>.


Module that performs a non-linear least squares minimization of the
spectroscopic and/or astrometric
data using the lmfit package.
"""
import time

import lmfit as lm
import numpy as np
import emcee as mc
import multiprocessing as mp

from modules.binary_system import BinarySystem

RV1 = RV2 = AS = False
LAS = LRV = 0


from scipy.stats import gaussian_kde, norm


class GaussianPrior:
    def __init__(self, mu, sigma):
        self.mu, self.sigma = mu, sigma
    def logpdf(self, x):
        return norm.logpdf(x, self.mu, self.sigma)


class KDEPrior:
    def __init__(self, samples):
        self._kde = gaussian_kde(samples)
    def logpdf(self, x):
        return float(self._kde.logpdf(x)[0])


class PriorSet(dict):
    """name -> prior object with a .logpdf(x) method. Missing keys contribute 0
    (i.e. the parameter keeps lmfit's implicit flat/bounds-only prior)."""
    def logpdf(self, params):
        return sum(prior.logpdf(params[name].value) for name, prior in self.items())


def priors_from_chain(flatchain, names, kind='kde'):
    priors = PriorSet()
    for name in names:
        samples = flatchain[name].values
        if kind == 'kde':
            priors[name] = KDEPrior(samples)
        else:
            priors[name] = GaussianPrior(*norm.fit(samples))
    return priors



def determine_datasets(data_dict):
    global RV1, RV2, AS
    RV1 = RV2 = AS = False
    global LAS, LRV
    LAS = LRV = 0

    rv1s = None
    rv2s = None
    aas = None
    if 'RV1' in data_dict and data_dict['RV1'] is not None:
        rv1s = data_dict['RV1']
        RV1 = True
        LRV = len(data_dict['RV1'])
    if 'RV2' in data_dict and data_dict['RV2'] is not None:
        rv2s = data_dict['RV2']
        RV2 = True
        LRV += len(data_dict['RV2'])
    if 'AS' in data_dict and data_dict['AS'] is not None:
        aas = data_dict['AS']
        AS = True
        LAS = 2 * len(data_dict['AS'])

    return rv1s, rv2s, aas


def build_master_param_set(guess_dict, lock_g=False, lock_q=False):
    params = lm.Parameters()
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
        ('gamma2', guess_dict['gamma2'][0], guess_dict['gamma2'][1]))
    if lock_g:
        params['gamma2'].set(expr='gamma1')
    if lock_q:
        params.add('q', value=params['k1'] / params['k2'], vary=False)
        params['k2'].set(expr='k1/q')
    if params['e'].value < 1e-8:
        params['e'].set(value=1e-8)
    return params


def restrict_to_RV(params, sb2: bool):
    """sb2 is auto-detected from data_dict via determine_datasets()/global RV2."""
    for key in ('i', 'Omega', 'mt', 'd'):
        params[key].set(vary=False)
    if not sb2:
        for key in ('k2', 'gamma2'):
            params[key].set(vary=False)
    return params


def restrict_to_AS(params):
    for key in ('k1', 'gamma1', 'k2', 'gamma2'):
        params[key].set(vary=False)
    return params


def sample_sigmas(p0, best_pars, sigmas):
    pert = np.zeros_like(p0)
    for i, sigma in enumerate(sigmas):
        pert[:, i] = np.random.normal(best_pars[i], sigma, size=p0.shape[0])
    return pert


SHARED_PARAMS = ('e', 'omega', 't0', 'p')

STAGE_CONFIG = {
    'RV': dict(restrict=restrict_to_RV, lnprob=lnprob_RV,
               data_key=lambda rv1s, rv2s, aas: (rv1s, rv2s)),
    'AS': dict(restrict=restrict_to_AS, lnprob=lnprob_AS,
               data_key=lambda rv1s, rv2s, aas: (aas,)),
}


def sequential_MCMC(guess_dict, error_dict, data_dict, direction='RV_AS',
                     steps=1000, walkers=100, burn=100, thin=1,
                     prior_kind='kde', lock_g=False, lock_q=False):
    """
    guess_dict/error_dict are assumed to already sit at a local minimum
    (e.g. from a prior leastsq run) — no local optimization is performed here.

    direction: 'RV_AS' fits RV first and feeds its posterior as priors into
    the astrometric fit; 'AS_RV' does the reverse.
    """
    rv1s, rv2s, aas = determine_datasets(data_dict)  # sets RV1, RV2, AS, LAS, LRV
    sb2 = RV2

    stage1_name, stage2_name = ('RV', 'AS') if direction == 'RV_AS' else ('AS', 'RV')

    def make_params(stage_name):
        params = build_master_param_set(guess_dict, lock_g, lock_q)
        if stage_name == 'RV':
            return restrict_to_RV(params, sb2)
        return restrict_to_AS(params)

    def init_walkers(params, nwalkers):
        """Ball of walkers around the (already-converged) guess, using the
        per-parameter errors already supplied in error_dict."""
        varying = [name for name in params if params[name].vary]
        best = np.array([guess_dict[name][0] for name in varying])
        sigma = np.array([error_dict[name] for name in varying])
        p0 = np.tile(best, nwalkers).reshape(nwalkers, len(varying))
        return p0 + sample_sigmas(p0, best, sigma)

    # ---- stage 1: flat-prior MCMC directly around the local minimum ----
    params1 = make_params(stage1_name)
    args1 = STAGE_CONFIG[stage1_name]['data_key'](rv1s, rv2s, aas)
    lnprob1 = STAGE_CONFIG[stage1_name]['lnprob']

    pos1 = init_walkers(params1, walkers)
    mcmc1 = lm.Minimizer(lnprob1, params1, fcn_args=(*args1, PriorSet()))
    result1 = mcmc1.emcee(steps=steps, nwalkers=walkers, burn=burn, thin=thin, pos=pos1)

    # ---- build priors for the parameters shared with stage 2 ----
    priors = priors_from_chain(result1.flatchain, SHARED_PARAMS, kind=prior_kind)

    # ---- stage 2: same local minimum, shared entries now carry priors ----
    params2 = make_params(stage2_name)
    pos2 = init_walkers(params2, walkers)
    mcmc2 = lm.Minimizer(STAGE_CONFIG[stage2_name]['lnprob'], params2,
                          fcn_args=(*STAGE_CONFIG[stage2_name]['data_key'](rv1s, rv2s, aas), priors))
    result2 = mcmc2.emcee(steps=steps, nwalkers=walkers, burn=burn, thin=thin, pos=pos2)

    return {stage1_name: result1, stage2_name: result2, 'priors': priors}


def _gaussian_lnlike(resid):
    return -0.5 * np.sum(resid ** 2)


def lnprob_RV(params, rv1s, rv2s, priors=PriorSet()):
    resid = residuals_RV(params, rv1s, rv2s)
    return _gaussian_lnlike(resid) + priors.logpdf(params)


def lnprob_AS(params, aas, priors=PriorSet()):
    resid = residuals_AS(params, aas)
    return _gaussian_lnlike(resid) + priors.logpdf(params)


def LMminimizer(guess_dict: dict, data_dict: dict, method: str = 'leastsq', hops: int = 10,
                steps: int = 1000, walkers: int = 100, burn: int = 100, thin: int = 1,
                as_weight: float = None, lock_g: bool = None, lock_q: bool = None):
    """
    Minimizes the provided data to a binary star model, with initial
    provided guesses and a search
    radius
    :param as_weight: weight to give to the astrometric data, optional.
    :param hops: int designating the number of hops if basinhopping is selected
    :param method: string to indicate what method to be used, 'leastsq' or
    'bqsinhopping' or 'emcee'
    :param guess_dict: dictionary containing guesses and 'to-vary' flags for
    the 11 parameters
    :param data_dict: dictionary containing observational data of RV and/or
    separations
    :param steps: integer giving the number of steps each walker in the MCMC
    should perform
    :param walkers: integer giving the number of independent walkers to be
    running
    :param burn: integer giving the number of samples to be discarded (
    "burned") at the start
    :param thin: integer indicating to accept only 1 every thin samples
    :param lock_g: boolean to indicate whether to lock gamma1 to gamma2
    :param lock_q: boolean to indicate whether to lock k2 to k1/q, and that q is supplied rather
    than k2 in that field.
    :return: result from the lmfit minimization routine. It is a
    MinimizerResult object.
    """

    # protect users
    if method == 'emcee' and burn >= steps:
        print('You are burning all steps of the MCMC chain! please put burn < '
              'steps')
        return

    # setup data for the solver
    rv1s, rv2s, aas = determine_datasets(data_dict)

    # setup Parameters object for the solver
    params = lm.Parameters()
    # populate with parameter data
    params.add_many(('e', guess_dict['e'][0], guess_dict['e'][1], 0, 1 - 1e-5),
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
                    ('gamma2', guess_dict['gamma2'][0], guess_dict['gamma2'][1]))

    if lock_g:
        params['gamma2'].set(expr='gamma1')
    if lock_q:
        params.add('q', value=params['k1'] / params['k2'], vary=False)
        params['k2'].set(expr='k1/q')

    # put e to a non-zero value to avoid conditioning problems in MCMC
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
        elif AS and ('q' in params.valuesdict().keys()) and params.valuesdict()['q'] != 0:
            params['i'].set(expr='180-180/pi*asin(sqrt(1-e**2)*k1*(q+1)/q*'
                                 '(p*86400/(2*pi*6.67430e-20*mt*1.9885e30))**(1/3))')

    elif AS:
        for key in 'k1', 'gamma1', 'k2', 'gamma2':
            params[key].set(vary=False)
    else:
        raise ValueError('No data supplied! Cannot minimize.\n')

    # build a minimizer object
    minimizer = lm.Minimizer(fcn2min, params, fcn_args=(rv1s, rv2s, aas, as_weight))
    print('Starting Minimization with {}{}{}...'.format('primary RV data, ' if RV1 else '',
                                                        'secondary RV data, ' if RV2 else '',
                                                        'astrometric data' if AS else ''))
    tic = time.time()
    if method == 'leastsq':
        result = minimizer.minimize()
    elif method == 'basinhopping':
        result = minimizer.minimize(method=method, disp=True, niter=hops, T=5,
                                    minimizer_kwargs={'method': 'Nelder-Mead'})
    elif method == 'emcee':
        localresult = minimizer.minimize()
        mcminimizer = lm.Minimizer(fcn2min, params=localresult.params,
                                   fcn_args=(rv1s, rv2s, aas, as_weight))
        print('Starting MCMC sampling using the minimized parameters...')
        #TODO: allow for non-uniform priors! lm.emcee only does uniform priors
        result = mcminimizer.emcee(steps=steps, nwalkers=walkers, burn=burn, thin=thin)
    else:
        print('this minimization method not implemented')
        return
    toc = time.time()
    print('Minimization Complete in {} s!\n'.format(np.round(toc - tic, 3)))
    lm.report_fit(result.params)
    rms_rv1, rms_rv2, rms_as = 0, 0, 0
    system = BinarySystem(result.params.valuesdict())
    if RV1:
        # weigh with number of points for RV1 data
        rms_rv1 = np.sqrt(
            np.sum((system.primary.radial_velocity_of_hjd(rv1s[:, 0]) - rv1s[:, 1]) ** 2) / len(
                rv1s[:, 1]))
    if RV2:
        # Same for RV2
        rms_rv2 = np.sqrt(
            np.sum((system.secondary.radial_velocity_of_hjd(rv2s[:, 0]) - rv2s[:, 1]) ** 2) / len(
                rv2s[:, 1]))
    if AS:
        # same for AS
        omc2E = np.sum((system.relative.east_of_hjd(aas[:, 0]) - aas[:, 1]) ** 2)
        omc2N = np.sum((system.relative.north_of_hjd(aas[:, 0]) - aas[:, 2]) ** 2)
        rms_as = np.sqrt((omc2E + omc2N) / LAS)
    print('Minimization complete, check parameters tab for resulting orbit!\n')
    return result, rms_rv1, rms_rv2, rms_as


def residuals_RV(params, rv1s, rv2s, weight=None):
    system = BinarySystem(params.valuesdict())
    if RV1:
        chisq_rv1 = (system.primary.radial_velocity_of_hjd(rv1s[:, 0]) - rv1s[:, 1]) / rv1s[:, 2]
        if weight:
            chisq_rv1 *= (1 - weight) * (LAS + LRV) / LRV
    else:
        chisq_rv1 = np.asarray([])
    if RV2:
        chisq_rv2 = (system.secondary.radial_velocity_of_hjd(rv2s[:, 0]) - rv2s[:, 1]) / rv2s[:, 2]
        if weight:
            chisq_rv2 *= (1 - weight) * (LAS + LRV) / LRV
    else:
        chisq_rv2 = np.asarray([])
    return np.concatenate((chisq_rv1, chisq_rv2))


def residuals_AS(params, aas, weight=None):
    system = BinarySystem(params.valuesdict())
    if not AS:
        return np.asarray([])
    chisq_east = (system.relative.east_of_hjd(aas[:, 0]) - aas[:, 1]) / aas[:, 3]
    chisq_north = (system.relative.north_of_hjd(aas[:, 0]) - aas[:, 2]) / aas[:, 4]
    if weight:
        chisq_east *= weight * (LAS + LRV) / LAS
        chisq_north *= weight * (LAS + LRV) / LAS
    return np.concatenate((chisq_east, chisq_north))


def fcn2min(params, rv1s, rv2s, aas, weight=None):
    """Unchanged public signature — now just delegates."""
    return np.concatenate((residuals_RV(params, rv1s, rv2s, weight),
                            residuals_AS(params, aas, weight)))
