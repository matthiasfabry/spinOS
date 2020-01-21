"""
This module provides functions to plot radial velocity curves and apparent orbits on the sky.
This module is developed with matplotlib 3.1.1. and numpy 1.18.1.

Author:
Matthias Fabry, Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
21 Jan 2020
"""
import corner
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection

plt.rc('text', usetex=True)
plt.rc('font', size=16)


def make_plots():
    """
    Makes two plot objects, and formats one to display RV curves, and the other for relative astrometry
    :return: RV figure, AS figure, RV axis, AS axis
    """
    fig1 = plt.figure(figsize=(10, 5))
    fig2 = plt.figure(figsize=(10, 5))
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111, aspect=1)
    setup_rvax(ax1)
    setup_asax(ax2)
    fig1.tight_layout()
    fig2.tight_layout()
    return fig1, fig2, ax1, ax2


def setup_asax(asax):
    """
    sets up a given axis for the plotting of the relative orbit
    :param asax: axis to format
    """

    asax.set_xlim((-10, 10))
    asax.invert_xaxis()
    asax.set_ylim((-10, 10))
    asax.set_xlabel(r'$East (mas)$')
    asax.set_ylabel(r'$North (mas)$')
    asax.axhline(linestyle=':', color='black')
    asax.axvline(linestyle=':', color='black')
    asax.grid()


def setup_rvax(rvax):
    """
    sets up a given axis for the plotting of radial velocity curves
    :param rvax: axis to format
    """
    rvax.set_xlabel(r'$orbital$ $phase$')
    rvax.set_ylabel(r'$RV (km s^{-1})$')
    rvax.set_xlim((-0.18, 1.18))
    rvax.set_ylim((-50, 50))
    rvax.axhline(linestyle=':', color='black')
    rvax.grid()


def plot_rv_curves(ax, system, rv1line=None, rv2line=None):
    """
    Plots RV curves for both components of a given system on a given axis
    :param rv1line: line object containing the rv1 data
    :param rv2line: line object containgin the rv2 data
    :param ax: axis to plot on
    :param system: system to get orbital parameters from
    """
    num = 200
    phases = np.linspace(-0.15, 1.15, num=num)
    vrads1 = system.primary.radial_velocity_of_phases(phases)
    vrads2 = system.secondary.radial_velocity_of_phases(phases)
    if rv1line is None:
        rv1line, = ax.plot(phases, vrads1, label='primary', color='b', ls='--')
    else:
        rv1line.set_ydata(vrads1)
    if rv2line is None:
        rv2line, = ax.plot(phases, vrads2, label='secondary', color='r', ls='--')
    else:
        rv2line.set_ydata(vrads2)
    ax.axis('auto')
    return rv1line, rv2line


def plot_relative_orbit(ax, system, asline=None, nodeline=None, peridot=None):
    """
    Plots the relative orbit of a given system on a given axis
    :param asline: line object that contains the as data
    :param nodeline: line object that contains the nodeline data
    :param peridot: line object that contains the periastron data
    :param ax: axis to plot on
    :param system: system to get orbital parameters from
    """
    num = 200
    ecc_anoms = np.linspace(0, 2 * np.pi, num)
    norths = system.relative.north_of_ecc(ecc_anoms)
    easts = system.relative.east_of_ecc(ecc_anoms)
    if asline is None:
        peridot, = ax.plot([system.relative.east_of_ecc(0)], [system.relative.north_of_ecc(0)], color='b', marker='s',
                           fillstyle='full', label='periastron', markersize=8)
        asline, = ax.plot(easts, norths, label='relative orbit', color='k')
        nodeline, = ax.plot([system.relative.east_of_true(-system.relative.omega),
                             system.relative.east_of_true(-system.relative.omega + np.pi)],
                            [system.relative.north_of_true(-system.relative.omega),
                             system.relative.north_of_true(-system.relative.omega + np.pi)], color='0.5', ls='--',
                            label='line of nodes')
    else:
        peridot.set_xdata(system.relative.east_of_ecc(0))
        peridot.set_ydata(system.relative.north_of_ecc(0))
        asline.set_xdata(easts)
        asline.set_ydata(norths)
        nodeline.set_xdata([system.relative.east_of_true(-system.relative.omega),
                            system.relative.east_of_true(-system.relative.omega + np.pi)])
        nodeline.set_ydata([system.relative.north_of_true(-system.relative.omega),
                            system.relative.north_of_true(-system.relative.omega + np.pi)])
    ax.axis('image')
    return asline, nodeline, peridot


def plot_dots(rvax, asax, phase, system, rv1dot=None, rv2dot=None, asdot=None):
    rv1 = system.primary.radial_velocity_of_phase(phase)
    if rv1dot is None:
        rv1dot, = rvax.plot(phase, rv1, 'rx')
    else:
        rv1dot.set_xdata(phase)
        rv1dot.set_ydata(rv1)
    rv2 = system.secondary.radial_velocity_of_phase(phase)
    if rv2dot is None:
        rv2dot, = rvax.plot(phase, rv2, 'bx')
    else:
        rv2dot.set_xdata(phase)
        rv2dot.set_ydata(rv2)
    N = system.relative.north_of_ph(phase)
    E = system.relative.east_of_ph(phase)
    if asdot is None:
        asdot, = asax.plot(E, N, 'ro')
    else:
        asdot.set_xdata(E)
        asdot.set_ydata(N)
    return rv1dot, rv2dot, asdot


def plot_rv_data(rvax, datadict, system, rv1dataline=None, rv2dataline=None):
    """
    Plots the given rv data for a given system on the given axes
    :param rv1dataline: line object that contains the rv1data
    :param rv2dataline: line object that contains the rv2data
    :param rvax: RV axis to plot RV data on
    :param datadict: dictionary with observational data
    :param system: system to get orbital parameters from
    """
    for key, data in datadict.items():
        if key == 'RV1' or key == 'RV2':
            phases, rv, err = system.create_phase_extended_RV(datadict[key], 0.15)
            if key == 'RV1':
                color = 'blue'
                if rv1dataline is None:
                    rv1dataline = rvax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o', ms=5,
                                                color=color)

                else:
                    rv1dataline.set_ydata(rv)
            else:
                color = 'red'
                if rv2dataline is None:
                    rv2dataline = rvax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o', ms=5,
                                                color=color)
                else:
                    rv2dataline.set_ydata(rv)
            rvax.axis('auto')
    return rv1dataline, rv2dataline


def plot_as_data(asax, datadict, asdataline=None, asellipses=None):
    """
    Plots the given as data for a given system on the given axes
    :param asdataline: line object containing the as data
    :param asellipses: collection containing the error ellipses
    :param asax: AS axis to plot astrometric data on
    :param datadict: dictionary with observational data
    """
    for key, data in datadict.items():
        if key == 'AS':
            if asdataline is None:
                asdataline, = asax.plot(data['easts'], data['norths'], 'r.', markersize=1)
            else:
                asdataline.set_xdata(data['easts'])
                asdataline.set_ydata(data['norths'])
                asax.remove_collection(asellipses)
            asellipses = EllipseCollection(2 * data['majors'], 2 * data['minors'], data['pas'] - 90,
                                           offsets=np.column_stack((data['easts'], data['norths'])),
                                           transOffset=asax.transData,
                                           units='x', edgecolors='r', facecolors=(0, 0, 0, 0))
            asax.add_collection(asellipses)
            plotmin = min(min(data['easts']), min(data['norths']))
            plotmax = max(max(data['easts']), max(data['norths']))
            asax.set_xlim([plotmax + 5, plotmin - 5])
            asax.set_ylim([plotmin - 5, plotmax + 5])
            asax.axis('image')
    return asdataline, asellipses


def plot_corner_diagram(mcmcresult):
    labels = []
    thruths = []
    for key in mcmcresult.var_names:
        if key == 'e':
            labels.append(r'$e$')
        elif key == 'i':
            labels.append(r'$i$ (deg)')
        elif key == 'omega':
            labels.append(r'$\omega$ (deg)')
        elif key == 'Omega':
            labels.append(r'$\Omega$ (deg)')
        elif key == 't0':
            labels.append(r'$T_0$ (MJD)')
        elif key == 'k1':
            labels.append(r'$K_1$ (km/s)')
        elif key == 'k2':
            labels.append(r'$K_2$ (km/s)')
        elif key == 'p':
            labels.append(r'$P$ (d)')
        elif key == 'gamma1':
            labels.append(r'$\gamma_1$ (km/s)')
        elif key == 'gamma2':
            labels.append(r'$\gamma_2$ (km/s)')
        elif key == 'd':
            labels.append(r'$d$ (pc)')

        if mcmcresult.params[key].vary:
            thruths.append(mcmcresult.params.valuesdict()[key])
    return corner.corner(mcmcresult.flatchain, labels=labels, truths=thruths)
