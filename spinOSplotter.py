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


def make_plots(system: System, datadict: dict):
    print('Starting to make plots')
    plt.rc('text', usetex=True)
    fig = plt.figure()
    plt.subplots_adjust(bottom=0.15)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122, aspect=1)
    plot_rv_curves(ax1, system.primary, system.secondary)
    plot_relative_orbit(ax2, system.relative)
    ax2.set_xlim(ax2.get_xlim()[::-1])
    # make_slider(fig, ax1, ax2, system)
    plot_data(ax1, ax2, datadict, system)


def plot_rv_curves(ax, primary_orbit, secondary_orbit):
    num = 200
    ecc_anoms = np.linspace(0, 2 * np.pi, num)
    vrads1, vrads2 = np.zeros(num), np.zeros(num)
    phases1, phases2 = np.zeros(num), np.zeros(num)
    for i in range(num):
        vrads1[i] = primary_orbit.radial_velocity_of_ecc_anom(ecc_anoms[i])
        vrads2[i] = secondary_orbit.radial_velocity_of_ecc_anom(ecc_anoms[i])
        phases1[i] = primary_orbit.phase_of_ecc_anom(ecc_anoms[i])
        phases2[i] = secondary_orbit.phase_of_ecc_anom(ecc_anoms[i])
    ax.set_title('Radial velocity')
    ax.plot(phases1, vrads1, label='primary')
    ax.plot(phases2, vrads2, label='secondary')
    ax.plot(phases1, np.zeros(num), 'g')
    ax.set_xlabel(r'$phase$')
    ax.set_ylabel(r'$RV (km s^{-1})$')
    ax.grid()
    ax.legend()


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


def plot_relative_orbit(ax, orb):
    ax.set_title('Relative orbit')
    num = 200
    ecc_anoms = np.linspace(0, 2 * np.pi, num)
    norths = orb.north_of_ecc(ecc_anoms)
    easts = orb.east_of_ecc(ecc_anoms)
    ax.plot(easts, norths)
    ax.plot([orb.east_of_true(-orb.omega), orb.east_of_true(-orb.omega + np.pi)],
            [orb.north_of_true(-orb.omega), orb.north_of_true(-orb.omega + np.pi)], 'k--',
            label='line of nodes')
    ax.plot([orb.east_of_ecc(0), orb.east_of_ecc(np.pi)], [orb.north_of_ecc(0), orb.north_of_ecc(np.pi)], 'k-',
            label='major axis')
    ax.plot([orb.east_of_ecc(0)], [orb.north_of_ecc(0)], marker='s', fillstyle='none')
    ax.set_xlabel('East')
    ax.set_ylabel('North')
    ax.axhline(linestyle=':', color='black')
    ax.axvline(linestyle=':', color='black')


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
