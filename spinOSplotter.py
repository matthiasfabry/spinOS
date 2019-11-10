"""
This module provides functions to plot radial velocity curves and apparent orbits on the sky.
This module is developed with matplotlib 3.1.1

Author:
Matthias Fabry, Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
29 Oct 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection

plt.rc('text', usetex=True)
plt.rc('font', size=16)


def make_plots():
    """
    Makes two plot objects, and formats one to display RV curves, and the other for relative astrometry
    :return: RV figure, AS figure, RV axis, AS axis
    """
    fig1 = plt.figure()
    fig2 = plt.figure()
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
    asax.autoscale()


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
    rvax.autoscale(axis='y')


def plot_rv_curves(ax, system):
    """
    Plots RV curves for both components of a given system on a given axis
    :param ax: axis to plot on
    :param system: system to get orbital parameters from
    """
    num = 200
    leftboundary_E = system.eccentric_anom_of_phase(-0.15)
    rightboundary_E = system.eccentric_anom_of_phase(1.15)
    ecc_anoms = np.linspace(leftboundary_E, rightboundary_E, num)
    vrads1, vrads2 = np.zeros(num), np.zeros(num)
    phases1, phases2 = np.zeros(num), np.zeros(num)
    for i in range(num):
        vrads1[i] = system.primary.radial_velocity_of_ecc_anom(ecc_anoms[i])
        vrads2[i] = system.secondary.radial_velocity_of_ecc_anom(ecc_anoms[i])
        phases1[i] = system.phase_of_ecc_anom(ecc_anoms[i])
        phases2[i] = system.phase_of_ecc_anom(ecc_anoms[i])
    ax.plot(phases1, vrads1, label='primary', color='r')
    ax.plot(phases2, vrads2, label='secondary', color='b')
    ax.legend()
    ax.relim()
    ax.autoscale_view()


def plot_relative_orbit(ax, system):
    """
    Plots the relative orbit of a given system on a given axis
    :param ax: axis to plot on
    :param system: system to get orbital parameters from
    """
    num = 200
    ecc_anoms = np.linspace(0, 2 * np.pi, num)
    norths = system.relative.north_of_ecc(ecc_anoms)
    easts = system.relative.east_of_ecc(ecc_anoms)
    ax.plot(easts, norths, label='relative orbit')
    ax.plot([system.relative.east_of_true(-system.relative.omega),
             system.relative.east_of_true(-system.relative.omega + np.pi)],
            [system.relative.north_of_true(-system.relative.omega),
             system.relative.north_of_true(-system.relative.omega + np.pi)], 'k--',
            label='line of nodes')
    ax.plot([system.relative.east_of_ecc(0), system.relative.east_of_ecc(np.pi)],
            [system.relative.north_of_ecc(0), system.relative.north_of_ecc(np.pi)], 'k-',
            label='major axis')
    ax.plot([system.relative.east_of_ecc(0)], [system.relative.north_of_ecc(0)], marker='s', fillstyle='none',
            label='periastron')
    ax.relim()
    ax.autoscale_view()


def plot_data(rvax, asax, datadict, system):
    """
    Plots the given data for a given system on the given axes
    :param rvax: RV axis to plot RV data on
    :param asax: AS axis to plot astrometric data on
    :param datadict: dictionary with observational data
    :param system: system to get orbital parameters from
    """
    for key, data in datadict.items():
        if key == 'RV1' or key == 'RV2':
            phases, rv, err = system.create_phase_extended_RV(datadict[key], 0.15)
            if key == 'RV1':
                color = 'orange'
            else:
                color = 'green'
            rvax.errorbar(phases, rv, yerr=err, ls='', capsize=0.1, marker='o', ms=2, color=color)
        elif key == 'AS':
            ellipses = EllipseCollection(data[:, 3], data[:, 4], data[:, 5] - 90,
                                         offsets=np.column_stack((data[:, 1], data[:, 2])), transOffset=asax.transData,
                                         units='x', edgecolors='r', facecolors='w')
            asax.add_collection(ellipses)
    asax.autoscale_view()
