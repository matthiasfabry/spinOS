"""
Copyright 2020, 2021 Matthias Fabry
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

presenting spinOS: the SPectroscopic and INterferometric Orbital Solution finder.

Goal:
    spinOS computes the best fit orbital solution given:
     1) radial and systemic velocity data for either or both components of a spectroscopic binary and/or
     2) astrometric data containing the separtions and error ellipses, and
     3) an initial guess of the parameters:
        - e        the eccentricity
        - i        the inclination (deg)
        - omega    the longitude of the periastron with respect to the ascending node of the secondary (deg)
        - Omega    the longitude of the ascending node of the seconday measured east of north (deg)
        - t0       the time of periastron passage (day)
        - k1       the semiamplitude of the radial velocity curve of the primary (km/s)
        - k2       the semiamplitude of the radial velocity curve of the secondary (km/s)
        - p        the period of the binary (day)
        - gamma1   the peculiar velocity of the primary (km/s)
        - gamma2   the peculiar velocity of the secondary (km/s)
        - d        the distance (pc)

This program minimizes a binary orbit model to your supplied data and afterwards plots the data and the minimized model.
The program then gives a best fit value for the parameters itemized above, as well as the component masses.
Errors can be calculated using a Markov Chain Monte Carlo (MCMC) method, the reported errors are half of the difference
between the 15.87 and 84.13 percentiles found in the MCMC sampling.

Usage:
To use spinOS, simply run:

 $ python3 spinOS.py -i <pointerfile> [-p] [-s] [-m] [-t <steps>]]

where:
 <pointerfile> is a plain text file with a list to your datafiles and guessfile,

and include the switches:
 -p to indicate only to plot the data with the model created with the guesses
 -s boolean to indicate whether your astrometric data is in separation/position angle or not
            (if you have east/north data, omit this)
 -m to calculate an MCMC error estimation
 -t <steps> (default = 1000): integer denoting the number of MCMC samples te be taken


A pointerfile looks like this:
    RV1file secondary_vels.txt
    # RV2file primary_vels.txt
    ASfile relative_astrometry.txt
    guessfile guesses.txt
the keys on the first column should be 'RV1file', 'RV2file', 'ASfile' and 'guessfile', each with the relative paths to
the respective data. guessfile is mandatory for spinOS to output something. You can exclude data (if unwanted or
unavailable) by commenting (#) lines or omitting them from the list.

spinOS expects the data to be in the following format: All data files should be plain text files, with:
for RV data:
 JD(days) RV(km/s) error_on_RV(km/s)
eg:
 45000 25.1 2.1
 45860 -4.2 1.1
 etc...

for AS data:
 JD(days) E_separation(mas) N_separation(mas) semimajor_ax_errorellipse(mas)
                                                            semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 -2.5 2.4 0.1 0.8 60
 48050 2.1 8.4 0.4 0.5 90
 etc...

or:
 JD(days) separation(mas) PA(deg) semimajor_ax_errorellipse(mas)
                                                            semiminor_ax_errorellipse(mas) angle_E_of_N_of_major_ax(deg)
eg:
 48000 3.5 316 0.1 0.8 60
 48050 8.7 76 0.4 0.5 90
 etc...

if relative astrometry is formatted in separation, position angle.

The guessfile must be formatted as:
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
 d 1250.0 False
All ten parameters should be guessed. For each parameter, a tag True/False should be supplied to decide for the
minimizer to vary this parameter. (So False means it will keep it fixed at the supplied value)


Dependencies:
    python 3.8.2
    numpy 1.18.1
    scipy 1.4.1
    lmfit 1.0.0
    matplotlib 3.1.3
    emcee 3.0.2 (if MCMC error calculation is performed)
    corner 2.0.1 (if MCMC error calculation is performed)

Author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
    3 Apr 2020

Version:
    2.3

Acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.

"""
import sys

from modules import spinOSio as spl, spinOSminimizer as spm, spinOSplotter as spp, binary_system


def run(opts):
    # os.system('clear')
    print('Hello, this is spinOS, your personal orbital solution finder. I will start promptly!\n')

    # parse command line options
    pointer = None
    plotonly = False
    do_seppa_conv = False
    domcmc = False
    steps = 1000
    for opt, arg in opts:
        if opt == '-h':
            print('spinOS.py -i <pointer> [-p] [-s] [-m [-t <steps>]]')
            sys.exit()
        elif opt == "-i":
            pointer = arg
        elif opt == "-p":
            plotonly = True
        elif opt == "-s":
            do_seppa_conv = True
        elif opt == "-m":
            domcmc = True
        elif opt == "-t":
            steps = arg

    # read in files

    wd, guessdict, datadict = spl.spinOSparser(pointer, do_seppa_conv)

    # compute best elements
    if plotonly:
        bestpars = guessdict
        redchisq = 0.
        dof = 0
        minimizationresult = None
    else:
        minimizationresult, _, _, _ = spm.LMminimizer(guessdict, datadict, domcmc, steps=steps)
        bestpars = minimizationresult.params.valuesdict()
        redchisq = minimizationresult.redchi
        dof = minimizationresult.nfree

    # compute model of these elements
    system = binary_system.System(bestpars)

    # plot resulting RV curve and resulting apparent orbit
    fig1, fig2, rvax, asax = spp.make_plots()
    spp.plot_relative_orbit(asax, system)
    spp.plot_rv_curves(rvax, system)
    spp.plot_rv_data(rvax, datadict, system)
    spp.plot_as_data(asax, datadict)
    rvax.legend()
    asax.legend()

    # calculate the resulting masses
    primary_mass = system.primary_mass()
    secondary_mass = system.secondary_mass()
    print('\nI have come to an optimal solution! These are:')
    print('P = {} days'.format(system.p))
    print('e = {}'.format(system.e))
    print('i = {} (deg)'.format(system.i))
    print('omega = {} (deg)'.format(system.secondary.omega))
    print('Omega = {} (deg)'.format(system.Omega))
    print('K1 = {} (km/s)'.format(system.primary.k))
    print('K2 = {} (km/s)'.format(system.secondary.k))
    print('t0 = {} (jd)'.format(system.t0))
    print('gamma1 = {} (km/s)'.format(system.primary.gamma))
    print('gamma2 = {} (km/s)'.format(system.secondary.gamma))
    print('d = {} (pc)'.format(system.d))
    print('M1 = {} (Msun)'.format(primary_mass))
    print('M2 = {} (Msun)\n'.format(secondary_mass))
    if not plotonly:
        print('The minimization algorithm stopped on a reduced chi squared of {}'.format(redchisq))
        print('with {} degrees of freedom'.format(dof))
    if domcmc:
        figc = spp.plot_corner_diagram(minimizationresult)
        figc.savefig(wd+'mcmc_corner.png')
    fig1.tight_layout()
    fig1.savefig(wd+'RVplot.png')
    fig2.tight_layout()
    fig2.savefig(wd+'ASplot.png')
    print('\nThis was spinOS, thanks for letting me help you!\n')
