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
spectroscopic and/or astrometric data using the lmfit package.
"""
import multiprocessing as mp
import time
from dataclasses import dataclass

import lmfit as lm
import numpy as np

from modules.binary_system import BinarySystem

from scipy.stats import gaussian_kde, norm


class GaussianPrior:
    def __init__(self, mu, sigma):
        self.mu, self.sigma = mu, sigma

    def logpdf(self, x):
        return norm.logpdf(x, self.mu, self.sigma)


class PriorSet(dict):
    """name -> prior object with a .logpdf(x) method. Missing keys contribute 0
    (i.e. the parameter keeps lmfit's implicit flat/bounds-only prior)."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.joint_priors = []

    def add_joint_prior(self, prior):
        self.joint_priors.append(prior)

    def logpdf(self, params):
        scalar_terms = sum(prior.logpdf(params[name].value) for name, prior in self.items())
        joint_terms = sum(prior.logpdf(params) for prior in self.joint_priors)
        return scalar_terms + joint_terms


class JointKDEPrior:
    def __init__(self, names, kde):
        self.names = tuple(names)
        self._kde = kde

    def logpdf(self, params):
        point = np.asarray([params[name].value for name in self.names], dtype=float)
        return float(self._kde.logpdf(point))


@dataclass(frozen=True)
class DatasetContext:
    has_rv1: bool
    has_rv2: bool
    has_as: bool
    n_as: int
    n_rv: int


def gaussian_priors_from_specs(prior_specs):
    """
    Build a PriorSet from a dict of name -> (mu, sigma) tuples. Missing or None values are ignored.
    :param prior_specs: dict of name -> (mu, sigma) tuples, where mu and sigma can be None to indicate no prior.
    :return:
    """
    priors = PriorSet()
    if prior_specs is None:
        return priors
    for name, spec in prior_specs.items():
        if spec is None:
            continue
        mu, sigma = spec
        if mu is None or sigma is None:
            continue
        sigma = float(sigma)
        if sigma <= 0:
            raise ValueError(f'Prior sigma for {name} must be positive.')
        priors[name] = GaussianPrior(float(mu), sigma)
    return priors


def joint_kde_prior_from_chain(flatchain, names):
    shared_names = [name for name in names if name in flatchain]
    if not shared_names:
        return None, ()

    matrix = np.vstack([np.asarray(flatchain[name].values, dtype=float) for name in shared_names])
    finite_mask = np.all(np.isfinite(matrix), axis=0)
    matrix = matrix[:, finite_mask]
    n_dim, n_samples = matrix.shape
    if n_samples <= n_dim:
        return None, ()
    try:
        kde = gaussian_kde(matrix)
    except np.linalg.LinAlgError:
        return None, ()
    return JointKDEPrior(shared_names, kde), tuple(shared_names)


def merge_priors(primary_priors=None, fallback_priors=None):
    priors = PriorSet()
    if fallback_priors is not None:
        priors.update(fallback_priors)
        if isinstance(fallback_priors, PriorSet):
            for joint_prior in fallback_priors.joint_priors:
                priors.add_joint_prior(joint_prior)
    if primary_priors is not None:
        priors.update(primary_priors)
        if isinstance(primary_priors, PriorSet):
            for joint_prior in primary_priors.joint_priors:
                priors.add_joint_prior(joint_prior)
    return priors


def mcmc_worker_count(nwalkers):
    cpu_total = 4
    return max(1, min(cpu_total, nwalkers))


def determine_datasets(data_dict):
    rv1s = None
    rv2s = None
    aas = None
    has_rv1 = has_rv2 = has_as = False
    n_as = n_rv = 0
    if 'RV1' in data_dict and data_dict['RV1'] is not None:
        rv1s = data_dict['RV1']
        has_rv1 = True
        n_rv = len(data_dict['RV1'])
    if 'RV2' in data_dict and data_dict['RV2'] is not None:
        rv2s = data_dict['RV2']
        has_rv2 = True
        n_rv += len(data_dict['RV2'])
    if 'AS' in data_dict and data_dict['AS'] is not None:
        aas = data_dict['AS']
        has_as = True
        n_as = 2 * len(data_dict['AS'])

    return rv1s, rv2s, aas, DatasetContext(has_rv1=has_rv1, has_rv2=has_rv2, has_as=has_as, n_as=n_as, n_rv=n_rv)


def build_master_param_set(guess_dict, lock_g=False, lock_q=False):
    params = lm.Parameters()
    params.add_many(
        ('e', guess_dict['e'][0], guess_dict['e'][1], 0, 1 - 1e-5),
        ('i', guess_dict['i'][0], guess_dict['i'][1], 0, 180),
        ('omega', guess_dict['omega'][0], guess_dict['omega'][1], 0, 360),
        ('Omega', guess_dict['Omega'][0], guess_dict['Omega'][1], 0, 360),
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


def constrain_params(params, context: DatasetContext):
    if context.has_rv1 and context.has_rv2:
        if not context.has_as:
            for key in 'd', 'i', 'Omega', 'mt':
                params[key].set(vary=False)
        else:
            if params['d'].vary:
                params['mt'].set(vary=False)
            elif params['mt'].vary:
                params['d'].set(vary=False)
    elif context.has_rv1:
        for key in 'k2', 'gamma2', 'd':
            params[key].set(vary=False)
        if not context.has_as:
            for key in 'i', 'Omega', 'mt':
                params[key].set(vary=False)
        elif context.has_as and ('q' in params.valuesdict().keys()) and params.valuesdict()['q'] != 0:
            params['i'].set(expr='180-180/pi*asin(sqrt(1-e**2)*k1*(q+1)/q*'
                                 '(p*86400/(2*pi*6.67430e-20*mt*1.9885e30))**(1/3))')
    elif context.has_as:
        for key in 'k1', 'gamma1', 'k2', 'gamma2':
            params[key].set(vary=False)
    else:
        raise ValueError('No data supplied! Cannot minimize or do MCMC.\n')


def restrict_to_RV(params, sb2: bool):
    """sb2 indicates whether RV2 observations are available."""
    for key in ('i', 'Omega', 'mt', 'd'):
        params[key].set(vary=False)
    if not sb2:
        for key in ('k2', 'gamma2'):
            params[key].set(vary=False)

def restrict_to_AS(params):
    for key in ('k1', 'gamma1', 'k2', 'gamma2'):
        params[key].set(vary=False)


def sample_sigmas(p0, sigmas):
    pert = np.zeros_like(p0)
    for i, sigma in enumerate(sigmas):
        pert[:, i] = np.random.normal(0, sigma, size=p0.shape[0])
    return pert


def _gaussian_lnlike(resid):
    return -0.5 * np.sum(resid ** 2)


def lnprob_RV(params, rv1s, rv2s, context, priors=None):
    resid = residuals_RV(params, rv1s, rv2s, context)
    return _gaussian_lnlike(resid) + log_prior(params, priors)


def lnprob_AS(params, aas, context, priors=None):
    resid = residuals_AS(params, aas, context)
    return _gaussian_lnlike(resid) + log_prior(params, priors)


STAGE_CONFIG = {
    'RV': dict(restrict=restrict_to_RV, lnprob=lnprob_RV,
               data_key=lambda rv1s, rv2s, aas: (rv1s, rv2s)),
    'AS': dict(restrict=restrict_to_AS, lnprob=lnprob_AS,
               data_key=lambda rv1s, rv2s, aas: (aas,)),
}


def varying_param_names(params):
    """
    Names of parameters emcee will actually sample — same rule lmfit's own
    Minimizer.prepare_fit() uses: vary=True and no expr constraint.
    Does NOT filter or copy `params` itself; `params` always keeps every
    parameter, fixed or varying.
    """
    return [name for name, par in params.items() if par.vary and par.expr is None]


def prepare_stage_params(guess_dict, data_dict, stage,
                         lock_g=False, lock_q=False):
    """
    Select the appropriate parameters for a specific stage of the fitting process.
    :param guess_dict: incoming local minimum
    :param data_dict: all datasets to compute residuals against
    :param stage: 'AS' or 'RV'
    :param lock_g: whether to lock the gammas
    :param lock_q: whether to lock mass ratio
    :return:
    """
    rv1s, rv2s, aas, context = determine_datasets(data_dict)
    sb2 = context.has_rv2
    if stage == 'RV' and context.n_rv == 0:
        raise ValueError('RV MCMC requested, but no RV data are supplied.')
    if stage == 'AS' and context.n_as == 0:
        raise ValueError('AS MCMC requested, but no astrometric data are supplied.')
    params = build_master_param_set(guess_dict, lock_g, lock_q)
    if stage == 'RV':
        restrict_to_RV(params, sb2)
    else:
        restrict_to_AS(params)
    constrain_params(params, context)
    args = STAGE_CONFIG[stage]['data_key'](rv1s, rv2s, aas) + (context,)
    return params, args


def stage_ndata(dataset, args):
    context = args[-1]
    return context.n_rv if dataset == 'RV' else context.n_as


def init_walkers(params, nwalkers, guess_dict, error_dict):
    """
    Ball of walkers around the supplied local minimum, sized per varying
    parameter from error_dict. Returns pos, shape (nwalkers, nvarying),
    covering ONLY the varying dimensions -- this is what emcee/lmfit expect
    for the `pos` argument. It has no bearing on `params`, which must
    separately already hold correct values for every parameter (fixed
    ones included) before being passed to Minimizer.
    """
    varying = varying_param_names(params)
    best = np.array([guess_dict[name][0] for name in varying])
    sigma = np.array([error_dict[name] for name in varying])
    p0 = np.tile(best, nwalkers).reshape(nwalkers, len(varying))
    delta = sample_sigmas(p0, sigma)
    return p0 + delta


def single_MCMC(guess_dict, error_dict, data_dict, dataset,
                steps=1000, walkers=100, burn=100, thin=1,
                priors=None, lock_g=False, lock_q=False, num_cores=None):
    """dataset: 'RV' or 'AS'. priors: optional PriorSet, empty (flat) by default."""
    priors = priors if priors is not None else PriorSet()
    params, args = prepare_stage_params(guess_dict, data_dict, dataset, lock_g=lock_g, lock_q=lock_q)
    lnprob = STAGE_CONFIG[dataset]['lnprob']

    pos = init_walkers(params, walkers, guess_dict, error_dict)
    print("Running MCMC sampling for {} dataset with {} walkers, {} steps, {} burn-in, {} thinning..."
          .format(dataset, walkers, steps, burn, thin))
    print('no of free parameters: {}'.format(len(varying_param_names(params))))
    worker_count = mcmc_worker_count(walkers) if num_cores is None else min(num_cores, mcmc_worker_count(walkers))

    if worker_count == 1:
        result = lm.Minimizer(lnprob, params, fcn_args=(*args, priors)).emcee(
            steps=steps, nwalkers=walkers, burn=burn, thin=thin, pos=pos)
    else:
        print('using {} worker processes for MCMC'.format(worker_count))
        with mp.get_context('spawn').Pool(processes=worker_count) as pool:
            result = lm.Minimizer(lnprob, params, fcn_args=(*args, priors)).emcee(
                steps=steps, nwalkers=walkers, burn=burn, thin=thin, pos=pos, workers=pool)

    result.ndata = stage_ndata(dataset, args)
    return result

def sequential_MCMC(guess_dict, error_dict, data_dict,
                    direction='RV_AS', priors=None,
                    steps1=1000, walkers1=100, burn1=100, thin1=1,
                    steps2=1000, walkers2=100, burn2=100, thin2=1,
                    lock_g=False, lock_q=False, num_cores=None):

    stage1_name, stage2_name = ('RV', 'AS') if direction == 'RV_AS' else ('AS', 'RV')

    common = dict(lock_g=lock_g, lock_q=lock_q)
    result1 = single_MCMC(guess_dict, error_dict, data_dict, dataset=stage1_name,
                          priors=priors, steps=steps1, walkers=walkers1, burn=burn1, thin=thin1, num_cores=num_cores,
                          **common)

    # print acceptance fractions
    print(lm.fit_report(result1))

    stage2_params, _ = prepare_stage_params(guess_dict, data_dict, stage2_name,
                                            lock_g=lock_g, lock_q=lock_q)
    posterior_joint_prior, joint_names = joint_kde_prior_from_chain(
        result1.flatchain, varying_param_names(stage2_params))
    stage2_priors = merge_priors(fallback_priors=priors)
    if posterior_joint_prior is not None:
        for name in joint_names:
            stage2_priors.pop(name, None)
        stage2_priors.add_joint_prior(posterior_joint_prior)
    result2 = single_MCMC(guess_dict, error_dict, data_dict, dataset=stage2_name,
                          steps=steps2, walkers=walkers2, burn=burn2, thin=thin2, priors=stage2_priors,
                          num_cores=num_cores, **common)

    print(lm.fit_report(result2))

    print("MCMC complete! errors are placed in the parameters tab")

    return {'stage1': result1, 'stage2': result2, 'priors': stage2_priors,
            'initial_priors': priors, 'posterior_joint_prior': posterior_joint_prior,
            'posterior_joint_names': joint_names}


def log_prior(params, priors=None):
    prior_set = priors if priors is not None else PriorSet()
    return prior_set.logpdf(params)


def LMminimizer(guess_dict: dict, data_dict: dict, method: str = 'leastsq', hops: int = 10,
                as_weight: float = None, lock_g: bool = None, lock_q: bool = None):
    """
    Minimizes the provided data to a binary star model, with initial
    provided guesses
    :param as_weight: weight to give to the astrometric data, optional.
    :param hops: int designating the number of hops if basinhopping is selected
    :param method: string to indicate what method to be used, 'leastsq' or 'bqsinhopping'
    :param guess_dict: dictionary containing guesses and 'to-vary' flags for the 11 parameters
    :param data_dict: dictionary containing observational data of RV and/or separations
    :param lock_g: boolean to indicate whether to lock gamma1 to gamma2
    :param lock_q: boolean to indicate whether to lock k2 to k1/q, and that q is supplied rather
    than k2 in that field.
    :return: result from the lmfit minimization routine. It is a
    MinimizerResult object.
    """

    # setup data for the solver
    rv1s, rv2s, aas, context = determine_datasets(data_dict)

    # setup Parameters object for the solver
    params = build_master_param_set(guess_dict, lock_g=lock_g, lock_q=lock_q)

    # build a minimizer object
    minimizer = lm.Minimizer(fcn2min, params, fcn_args=(rv1s, rv2s, aas, context, as_weight))
    print('Starting Minimization with {}{}{}...'.format('primary RV data, ' if context.has_rv1 else '',
                                                        'secondary RV data, ' if context.has_rv2 else '',
                                                        'astrometric data' if context.has_as else ''))
    tic = time.time()
    if method == 'leastsq':
        result = minimizer.minimize()
    elif method == 'basinhopping':
        result = minimizer.minimize(method=method, disp=True, niter=hops, T=5,
                                    minimizer_kwargs={'method': 'Nelder-Mead'})
    else:
        print('this minimization method not implemented')
        return
    toc = time.time()
    print('Minimization Complete in {} s!\n'.format(np.round(toc - tic, 3)))
    lm.report_fit(result.params)
    rms_rv1, rms_rv2, rms_as = 0, 0, 0
    system = BinarySystem(result.params.valuesdict())
    if context.has_rv1:
        # weigh with number of points for RV1 data
        rms_rv1 = np.sqrt(
            np.sum((system.primary.radial_velocity_of_hjd(rv1s[:, 0]) - rv1s[:, 1]) ** 2) / len(
                rv1s[:, 1]))
    if context.has_rv2:
        # Same for RV2
        rms_rv2 = np.sqrt(
            np.sum((system.secondary.radial_velocity_of_hjd(rv2s[:, 0]) - rv2s[:, 1]) ** 2) / len(
                rv2s[:, 1]))
    if context.has_as:
        # same for AS
        omc2E = np.sum((system.relative.east_of_hjd(aas[:, 0]) - aas[:, 1]) ** 2)
        omc2N = np.sum((system.relative.north_of_hjd(aas[:, 0]) - aas[:, 2]) ** 2)
        rms_as = np.sqrt((omc2E + omc2N) / context.n_as)
    print('Minimization complete, check parameters tab for resulting orbit!\n')
    print(lm.fit_report(result))
    return result, rms_rv1, rms_rv2, rms_as


def residuals_RV(params, rv1s, rv2s, context: DatasetContext, weight=None):
    system = BinarySystem(params.valuesdict())
    if context.has_rv1:
        chisq_rv1 = (system.primary.radial_velocity_of_hjd(rv1s[:, 0]) - rv1s[:, 1]) / rv1s[:, 2]
        if weight and context.n_rv > 0:
            chisq_rv1 *= (1 - weight) * (context.n_as + context.n_rv) / context.n_rv
    else:
        chisq_rv1 = np.asarray([])
    if context.has_rv2:
        chisq_rv2 = (system.secondary.radial_velocity_of_hjd(rv2s[:, 0]) - rv2s[:, 1]) / rv2s[:, 2]
        if weight and context.n_rv > 0:
            chisq_rv2 *= (1 - weight) * (context.n_as + context.n_rv) / context.n_rv
    else:
        chisq_rv2 = np.asarray([])
    return np.concatenate((chisq_rv1, chisq_rv2))


def residuals_AS(params, aas, context: DatasetContext, weight=None):
    system = BinarySystem(params.valuesdict())
    if not context.has_as:
        return np.asarray([])
    chisq_east = (system.relative.east_of_hjd(aas[:, 0]) - aas[:, 1]) / aas[:, 3]
    chisq_north = (system.relative.north_of_hjd(aas[:, 0]) - aas[:, 2]) / aas[:, 4]
    if weight and context.n_as > 0:
        chisq_east *= weight * (context.n_as + context.n_rv) / context.n_as
        chisq_north *= weight * (context.n_as + context.n_rv) / context.n_as
    return np.concatenate((chisq_east, chisq_north))


def fcn2min(params, rv1s, rv2s, aas, context: DatasetContext, weight=None):
    return np.concatenate((residuals_RV(params, rv1s, rv2s, context, weight),
                            residuals_AS(params, aas, context, weight)))
