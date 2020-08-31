"""
This is spinOS, the SPectroscopic and INterferometric Orbiral Solution finder, put in a graphical user interface

Goal:
This is a graphical user interface implementation of the command line version spinOS. It uses (spectroscopic or
otherwise) Radial Velocity (RV) measurements of either of the two components of the binary, as well as (
interferometric or otherwise) relative astrometry (AS) measurements of the binary. You then need to supply a guess
for the parameters that define the binary orbit:
    -) e:       the eccentricity of the orbit
    -) i:       the inclination of the orbit (with respect to the plane of the sky)
    -) omega:   the argument of periastron of the secondary, with respect to its ascending node
    -) Omega:   the longitude of the ascending node of the secondary, measured East from North
    -) t0:      the time of periastron passage (this number should be between 0 and the supposed period)
    -) p:       the period of the binary
    -) d:       the distance to the system
and:
    -) k1:      the semiamplitude of the RV curve of the primary
    -) k2:      the semiamplitude of the RV curve of the secondary
    -) gamma1:  the peculiar velocity of the primary
    -) gamma2:  the peculiar velocity of the secondary
if you have an SB2 with or without astrometric data, or:
    -) k1:      the semiamplitude of the RV curve of the primary
    -) gamma1:  the peculiar velocity of the primary
    -) mt:      the total dynamical mass of the system (which sets the apparent size of the orbit)
if you an SB1 with or without astrometric data, or:
    -) mt:      the total dynamical mass of the system (which sets the apparent size of the orbit)
if you have an astrometric orbit only.

This application allows for easy plotting of data and models, as well as minimization of the model to your supplied
data. The program then gives a best fit value for the parameters itemized above, as well as the component masses.
Error are estimated as the diagonal elements of the correlation matrix, or as half of the difference
between the 15.87 and 84.13 percentiles found in an Markov Chain Monte Carlo sampling.

Usage:
To start spinOS, simply run:
 python3 spinOS.py [<dir>]

In the data tab, put your working directory, where all your data files are
    (typically <objectname/> don't forget the slash!)

The application expects the data to be in the following format: All data files should be plain text files, formatted as:
for RV data:
 JD(days) RV(km/s) error_on_RV(km/s)
eg:
 45000 25.1 2.1
 45860 -4.2 1.1
 etc...

for AS data:
either:
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

for the guess file (which is optional), format should be eg:
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
All eleven parameters should be guessed. (for their meaning see above)


In the System/parameters tab, you can play with the systemic parameters, and see immediately the changes of the orbit
on the plots. With the checkbuttons, indicate which parameters should be minimized. Below, some inferred parameters of
the model are presented. Use the load guesses button to load all guesses from your guessfile indicated in the data tab.
The save buttons save either the guesses to guesses.txt (warning: this overwrites!) or minimized parameters to
params_run<i>.txt (i is number of minimization runs performed in this session).

In the minimize tab, you can apply a custum weighting of the astrometric data to the chi-squared value (typically, you
would want to increase this).
You can minimize the model to the selected data with the minimize button, with or without an mcmc error estimation.
If the last minimization run contained an MCMC analysis, you can create a corner plot with the button provided. It will
be saved at corner<i>.png (i is run number).

In the plot controls tab, various checkbuttons are provided to plot certain elements on the plot windows on the right.
The phase slider allows for overplotting a dot at the phase indicated (for illustrative purposes, eg, connecting the
apparent orbit with the RV plot).


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
    5 Jun 2020

Version:
    2.4

Acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.
"""
import sys

import modules.spinOSApp as app
import modules.spinOScommandline as cmline
import getopt


def hhelp():
    """
    help function
    """
    print('spinOS.py [<dir>]                                             for the GUI')
    print('or')
    print('spinOS.py [<dir>] -i <pointer> [-p] [-s] [-m [-t <steps>]]    for the commandline utility.')


if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], "w:hi:psmt:")
    except getopt.GetoptError:
        hhelp()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            hhelp()
            sys.exit()

    if not opts:
        try:
            wd = args[0]
        except IndexError:
            wd = None
        app.run(wd)
    else:
        cmline.run(opts)
