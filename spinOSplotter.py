"""
This module provides functions to plot radial velocity curves and apparent orbits on the sky.
It requires Orbit objects to function.

This module is developed with matplotlib 3.1.1


Author:
Matthias Fabry, Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
29 Oct 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from orbit import Orbit, System
from matplotlib.widgets import Slider
from matplotlib.patches import Ellipse

plt.rc('text', usetex=True)
plt.rc('font', size=16)


def setup_relax(relax):
    relax.set_xlim((10, -10))
    relax.set_ylim((-10, 10))
    relax.set_xlabel(r'East (mas)')
    relax.set_ylabel(r'North (mas)')
    relax.axhline(linestyle=':', color='black')
    relax.axvline(linestyle=':', color='black')
    relax.grid()


def setup_rvax(rvax):
    rvax.set_xlabel(r'$phase$')
    rvax.set_ylabel(r'$RV (km s^{-1})$')
    rvax.set_xlim((-0.05, 1.05))
    rvax.set_ylim((-50, 50))
    rvax.axhline(linestyle=':', color='black')
    rvax.grid()


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
    ax.clear()
    setup_rvax(ax)
    num = 200
    ecc_anoms = np.linspace(0, 2 * np.pi, num)
    vrads1, vrads2 = np.zeros(num), np.zeros(num)
    phases1, phases2 = np.zeros(num), np.zeros(num)
    for i in range(num):
        vrads1[i] = system.primary.radial_velocity_of_ecc_anom(ecc_anoms[i])
        vrads2[i] = system.secondary.radial_velocity_of_ecc_anom(ecc_anoms[i])
        phases1[i] = system.primary.phase_of_ecc_anom(ecc_anoms[i])
        phases2[i] = system.secondary.phase_of_ecc_anom(ecc_anoms[i])

    ax.plot(phases1, vrads1, label='primary')
    ax.plot(phases2, vrads2, label='secondary')
    ax.plot(phases1, np.zeros(num), 'g')
    ax.set_ylim((1.1*min(min(vrads1), min(vrads2)), 1.1*max(max(vrads1), max(vrads2))))
    ax.legend()


def plot_relative_orbit(ax, system):
    ax.clear()
    setup_relax(ax)
    num = 200
    ecc_anoms = np.linspace(0, 2 * np.pi, num)
    norths = system.relative.north_of_ecc(ecc_anoms)
    easts = system.relative.east_of_ecc(ecc_anoms)
    ax.plot(easts, norths)
    ax.plot([system.relative.east_of_true(-system.relative.omega),
             system.relative.east_of_true(-system.relative.omega + np.pi)],
            [system.relative.north_of_true(-system.relative.omega),
             system.relative.north_of_true(-system.relative.omega + np.pi)], 'k--',
            label='line of nodes')
    ax.plot([system.relative.east_of_ecc(0), system.relative.east_of_ecc(np.pi)],
            [system.relative.north_of_ecc(0), system.relative.north_of_ecc(np.pi)], 'k-',
            label='major axis')
    ax.plot([system.relative.east_of_ecc(0)], [system.relative.north_of_ecc(0)], marker='s', fillstyle='none')
    ax.set_xlim((1.1*max(easts), 1.1*min(easts)))
    ax.set_ylim((1.1*min(norths), 1.1*max(norths)))


def plot_data(rvax, relax, datadict, system):
    for key, data in datadict.items():
        if key == 'RV1':
            phases = system.primary.phase_of_hjd(data[:, 0])
            rvax.errorbar(phases, data[:, 1], yerr=data[:, 2], ls='', marker='o', ms=2)
        elif key == 'RV2':
            phases = system.secondary.phase_of_hjd(data[:, 0])
            rvax.errorbar(phases, data[:, 1], yerr=data[:, 2], ls='', marker='o', ms=2)
        elif key == 'AS':
            ellipses = [
                Ellipse((data[i, 1], data[i, 2]), data[i, 3], data[i, 4], 90 - data[i, 5], fill=False, color='r') for i
                in range(len(data[:, 0]))]
            for e in ellipses:
                relax.add_artist(e)


def make_slider(fig, rvax, relax, system):
    marker_point1, = rvax.plot(0, system.primary.radial_velocity_of_ecc_anom(0), 'ro')
    marker_point2, = rvax.plot(0, system.secondary.radial_velocity_of_ecc_anom(0), 'ro')
    marker_point3, = relax.plot(system.relative.east_of_ecc(0), system.relative.north_of_ecc(0), 'ro')

    # noinspection PyTypeChecker
    axphase = plt.axes([0.05, 0.05, 0.9, 0.03])
    sphase = Slider(axphase, 'phase', np.float64(0.), np.float64(1.), valinit=0, valstep=0.005)

    def update(dummy):  # dummy argument is required by matplotlib
        phase = sphase.val
        newvrad1, newtheta1, newecc_anom1 = system.primary.radial_velocity_of_phase(phase, getAngles=True)
        newvrad2, newtheta2, newecc_anom2 = system.secondary.radial_velocity_of_phase(phase, getAngles=True)
        newecc_anomr = system.relative.eccentric_anom_of_phase(phase)
        marker_point1.set_xdata(phase)
        marker_point1.set_ydata(newvrad1)
        marker_point2.set_xdata(phase)
        marker_point2.set_ydata(newvrad2)
        marker_point3.set_xdata(system.relative.east_of_ecc(newecc_anomr))
        marker_point3.set_ydata(system.relative.north_of_ecc(newecc_anomr))
        fig.canvas.draw()

    sphase.on_changed(update)
