import lmfit
import scipy
from multiprocessing.pool import Pool
import numpy as np
import os


def cost_fun(params, **kwargs):
    return scipy.optimize.rosen_der([params['a'], params['b']])


class FastPool:
    # vectorised pool, it reduces multiprocessing overhead
    def __init__(self, pool):
        self.pool = pool

    def map(self, f, arg):
        npool = len(self.pool._pool)

        arg_list = np.array_split(list(arg), npool)
        vf = np.vectorize(f, signature='(n)->()')

        res_list = self.pool.map(vf, arg_list)

        return np.hstack(res_list)


if __name__ == '__main__':
    params = lmfit.Parameters()
    params.add('a', 1, min=-5, max=5, vary=True)
    params.add('b', 1, min=-5, max=5, vary=True)

    fitter = lmfit.Minimizer(cost_fun, params)
    with Pool(processes=4) as pool:
        MC_results = fitter.emcee(workers=FastPool(pool))

# not worth for this rosenbrock function...
