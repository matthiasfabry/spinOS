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


def setup_relax(relax):
    relax.set_xlim((-10, 10))
    relax.invert_xaxis()
    relax.set_ylim((-10, 10))
    relax.set_xlabel(r'$East (mas)$')
    relax.set_ylabel(r'$North (mas)$')
    relax.axhline(linestyle=':', color='black')
    relax.axvline(linestyle=':', color='black')
    relax.grid()
    relax.autoscale()


def setup_rvax(rvax):
    rvax.set_xlabel(r'$orbital$ $phase$')
    rvax.set_ylabel(r'$RV (km s^{-1})$')
    rvax.set_xlim((-0.03, 1.03))
    rvax.set_ylim((-50, 50))
    rvax.axhline(linestyle=':', color='black')
    rvax.grid()
    rvax.autoscale(axis='y')


def make_plots():
    fig1 = plt.figure()
    fig2 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111, aspect=1)
    setup_rvax(ax1)
    setup_relax(ax2)
    fig1.tight_layout()
    fig2.tight_layout()
    return fig1, fig2, ax1, ax2


def plot_rv_curves(ax, system):
    num = 200
    ecc_anoms = np.linspace(0, 2 * np.pi, num)
    vrads1, vrads2 = np.zeros(num), np.zeros(num)
    phases1, phases2 = np.zeros(num), np.zeros(num)
    for i in range(num):
        vrads1[i] = system.primary.radial_velocity_of_ecc_anom(ecc_anoms[i])
        vrads2[i] = system.secondary.radial_velocity_of_ecc_anom(ecc_anoms[i])
        phases1[i] = system.phase_of_ecc_anom(ecc_anoms[i])
        phases2[i] = system.phase_of_ecc_anom(ecc_anoms[i])
    ax.plot(phases1, vrads1, label='primary')
    ax.plot(phases2, vrads2, label='secondary')
    ax.legend()
    ax.relim()
    ax.autoscale_view()


def plot_relative_orbit(ax, system):
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


def plot_data(rvax, relax, datadict, system):
    for key, data in datadict.items():
        if key == 'RV1':
            phases = system.phase_of_hjds(data[:, 0])
            rvax.errorbar(phases, data[:, 1], yerr=data[:, 2], ls='', capsize=0.1, marker='o', ms=2)
        elif key == 'RV2':
            phases = system.phase_of_hjds(data[:, 0])
            rvax.errorbar(phases, data[:, 1], yerr=data[:, 2], ls='', marker='o', ms=2)
        elif key == 'AS':
            ellipses = EllipseCollection(data[:, 3], data[:, 4], data[:, 5] - 90,
                                         offsets=np.column_stack((data[:, 1], data[:, 2])), transOffset=relax.transData,
                                         units='x', edgecolors='r', facecolors='w')
            relax.add_collection(ellipses)
    relax.autoscale_view()
