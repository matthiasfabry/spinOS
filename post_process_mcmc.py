import numpy as np
import h5py as h5
import corner
import matplotlib.pyplot as plt

import modules.utils as utils
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False
plt.switch_backend('agg')
plt.style.use('rsc/spinOS.mplstyle')

#%%

system = ''
results = h5.File(f'{system}/mcmc_chain.hdf5')

#%% make a custom corner plot
samples = results['flatchain'][:]
auto_corr = results['autocorr'][:]
maxtau = np.nanmax(auto_corr)
print(n_tau := (len(samples)/maxtau))
burn = int(3*maxtau)
trim = int(len(samples) // 10000)
print((len(samples) - 3*maxtau)/trim)
labels = [utils.convert_to_latex_string(name) for name in results['var_names'][:]]

medians = np.median(samples[burn::trim], axis=0)
sigma_plus = np.percentile(samples[burn::trim], 84, axis=0) - medians
sigma_min = medians - np.percentile(samples[burn::trim], 16, axis=0)

fig, axs = plt.subplots(nrows=len(labels), ncols=len(labels), figsize=(10, 10),
                        gridspec_kw={'hspace': 0.05, 'wspace': 0.05})

corner.corner(samples[burn::trim], labels=labels, labelpad=0.1,
              truths=medians, fig=fig, quantiles=[0.16, 0.5, 0.84], range=[.99 for _ in labels])
plt.savefig(f'{system}/mcmc_chain.png', dpi=300)
plt.close()


