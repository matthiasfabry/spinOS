# spinOS: the SPectroscopic and INterferometric Orbital Solution finder.

[![GitHub license](https://img.shields.io/github/license/matthiasfabry/spinOS)](https://github.com/matthiasfabry/spinOS/blob/master/COPYING.txt)
[![Other](https://img.shields.io/badge/ASCL-2102.001-blue)](https://ascl.net/2102.001)
![Github version](https://img.shields.io/badge/version-v3.0.0-red)

## Goal:

spinOS calculates binary orbital elements. Given a set of radial velocity measurements of a
spectroscopic binary and/or relative position measurement of an astrometric binary, spinOS fits an
orbital model by minimizing a chi squared metric. You need to supply an initial guess for the
parameters that define the binary orbit:

- p:       the period of the binary
- e:       the eccentricity of the orbit
- i:       the inclination of the orbit (with respect to the plane of the sky)
- omega:   the argument of periastron of the secondary, with respect to its ascending node
- Omega:   the longitude of the ascending node of the secondary, measured East from North
- t0:      the time of periastron passage
- d:       the distance to the system
- k1:      the semiamplitude of the RV curve of the primary
- k2:      the semiamplitude of the RV curve of the secondary
- gamma1:  the peculiar velocity of the primary
- gamma2:  the peculiar velocity of the secondary
- mt:      the total mass of the system

This package provides a GUI to easily visualize your data. The GUI allows for easy plotting of data
and models, as well as minimization of the model to your supplied data. The program then gives a
best fit value for the parameters itemized above.

## Usage:

To use the GUI, simply run:

    python spinOS.py [dir]

## Data
In the GUI data tab, put your working directory, where all your data files are relative to where you
launched spinOS from. You can select files to be loaded in, or add Datasets and data entries
manually in the respective tabs of the Data Manager. When loading in data, see formatting below.

## Parameters
In the System/parameters tab, you can play with the systemic parameters, and see the changes of the
orbit on the plots by pressing 'refresh'. With the checkbuttons, indicate which parameters should be
minimized. Below, some inferred parameters of the model are presented. Use the load guesses button
to load all guesses from your guessfile indicated in the data tab. The save buttons save either the
guesses to files you specify in the output names tab (warning: may overwrite!)

## Fitting
In the minimize tab, you can apply a custom weighting of the astrometric data to the chi-squared
value (typically, you would want to increase this if you trust astrometry better than spectroscopy,
and have only few astrometric measurements). You can minimize the model to the selected data with
the minimize button, selecting a method first. Levenberg-Marquardt does local non-linear least
squares minimization, while basinhopping tries to find a global minimum by iterating local
Nelder-Mead simplex minimizations. Note that basinhopping can take considerably longer with lots of
free parameters. When minimizing, errors are estimated as the diagonal elements of the
correlation matrix.

## MCMC
To better estimate uncertainties (taking into account correlations), you can go to the MCMC tab.
Starting from a local minimum with estimated errors, MCMC recreated the parameters' posterior
distributions using an ensemble sampler. When MCMC sampling, the first _burn_ (default = 100)
samples are discarded, and then only 1 every _thin_ (default = 1) samples are retained in the final
results. The philosophy behind this is that the underlying sampler does not draw independent samples
from the posterior distribution, it first needs to 'settle' to the maximum likelihood region (hence
the burn-in), and then a random walk will only yield independent results twice every time the
characteristic autocorrelation time has passed (hence the thinning). These parameters are difficult
to estimate beforehand. While a burn-in is necessary, thinning is not needed, since the reported
errors are calculated from the quantiles of the chain's samples. Thinning would only be needed
for visualization purposes, to reduce the number of points in the corner plot.
The MCMC tab also allows for the use of gaussian priors on the parameters (e.g., for implementing constraints from other
missions/datasets not modeled in spinOS).
Finally, in V3.0.0 MCMC is reimplemented in a separate tab to satisfy Bayesian statistics, and removes the need to
introduce statistical weights to the chi-squared metric, which is a frequentist approach (as was implemented in spinOS
v2.X).


## Plotting
In the plot controls tab, various checkbuttons are provided to plot certain elements on the plot
windows on the right. The phase slider allows for overplotting a dot at the phase indicated (for
illustrative purposes, eg, for visually connecting the apparent orbit with the RV plot).

## Data formatting:

This application expects the data to be in the following format: All data files should be plain text
files, formatted as:
for RV data:

    JD(days) RV(km/s) error_on_RV(km/s)

_e.g._:

    45000 25.1 2.1
    45860 -4.2 1.1
    etc...

for AS data:
either (set button E/N):

    JD(days) E_separation(mas) N_separation(mas) semimajor_ax_errorellipse(mas) semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)

_e.g._:

    48000 -2.5 2.4 0.1 0.8 60
    48050 2.1 8.4 0.4 0.5 90
    etc...

or (set button Sep/PA):

    JD(days) separation(mas) PA(deg) semimajor_ax_errorellipse(mas) semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)

_e.g._:

    48000 3.5 316 0.1 0.8 60
    48050 8.7 76 0.4 0.5 90
    etc...

for the guess file, format should be _e.g._:

    e 0.648 True
    i 86.53 True
    omega 211.0 True
    Omega 67.3 True
    t0 56547.1 True
    k1 31.0 False
    k2 52.0 True
    p 3252.0 True
    gamma1 15.8 False
    gamma2 5.6 False
    mt 30.0 True

All eleven parameters should be guessed if you load in guesses. (for their meaning see above)

## Dependencies:

A LaTeX distribution that allows matplotlib.rc(usetex=True)
the supplied yaml file can be used to configure your python environment
e.g. run 
`conda create -n spinOS -f environment.yaml`

## Author:

Matthias FABRY

Department of Astrophysics and Planetary Science, Villanova University, PA, USA

## Date:

7 October 2026

## Licence:

Copyright Matthias Fabry. This software is released under the GNU GPL-3.0-or-later License, which is attached to this
repository

## Acknowledgements:

We thank the authors of scipy, lmfit, corner for the development of their respective packages.