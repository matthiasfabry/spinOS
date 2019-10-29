"""
This module plots a Radial velocity curve, as well as the corresponding orbit.
The module is used as follows:
$ python spinOSplotter.py a e i omega Omega T gamma P
with:
a       the semimajor axis of the orbit (in AU)
e       its eccentricity
i       its inclination (in degrees)
omega   its argument of periastron (in degrees)
Omega   its longitude of the ascending node (in degrees)
T       its time of periastron crossing (in HJD or days)
gamma   the systemic velocity (in km/s)
P       its period (in years)

This module is developed with scipy 1.3.1 and matplotlib 3.1.1


Author:
Matthias Fabry, Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
22 Oct 2019
"""
import plot_3d_orbit
import numpy as np
import matplotlib.pyplot as plt
import sys
from orbit import Orbit
from matplotlib.widgets import Slider
from matplotlib.patches import Ellipse


def make_plots(relative_orbit, primary_orbit, secondary_orbit, datadict: dict):
    print('Starting to make plots')
    plt.rc('text', usetex=True)
    fig = plt.figure()
    plt.subplots_adjust(bottom=0.15)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132, aspect=1)
    ax3 = fig.add_subplot(133, aspect=1)
    plot_rv_curves(ax1, primary_orbit, secondary_orbit)
    plot_orbit(ax2, primary_orbit, 'Absolute orbits')
    plot_orbit(ax2, secondary_orbit)
    ax2.set_xlim(ax2.get_xlim()[::-1])
    plot_orbit(ax3, relative_orbit, 'Relative orbit')
    ax3.set_xlim(ax3.get_xlim()[::-1])
    # make_slider(fig, ax1, ax2, ax3, relative_orbit, primary_orbit, secondary_orbit)
    plot_data(ax1, ax3, datadict, primary_orbit, secondary_orbit)


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


def plot_data(rvax, relax, datadict, primary_orbit: Orbit, secondary_orbit: Orbit):
    for key, data in datadict.items():
        if key == 'RV1':
            phases = primary_orbit.phase_of_hjd(data[:, 0])
            rvax.errorbar(phases, data[:, 1], yerr=data[:, 2], ls='', marker='o', ms=2)
        elif key == 'RV2':
            phases = secondary_orbit.phase_of_hjd(data[:, 0])
            rvax.errorbar(phases, data[:, 1], yerr=data[:, 2], ls='', marker='o', ms=2)
        elif key == 'AS':
            ellipses = [
                Ellipse((data[i, 1], data[i, 2]), data[i, 3], data[i, 4], 90 - data[i, 5], fill=False, color='r') for i
                in range(len(data[:, 0]))]
            for e in ellipses:
                relax.add_artist(e)


def plot_orbit(ax, orb, title: str = None):
    if title is not None:
        ax.set_title(title)
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


def make_slider(fig, rvax, absax, relax, relative_orbit: Orbit, primary_orbit: Orbit, secondary_orbit: Orbit):
    marker_point1, = rvax.plot(0, primary_orbit.radial_velocity_of_ecc_anom(0), 'ro')
    marker_point2, = rvax.plot(0, secondary_orbit.radial_velocity_of_ecc_anom(0), 'ro')
    marker_point3, = relax.plot(relative_orbit.east_of_ecc(0), relative_orbit.north_of_ecc(0), 'ro')
    marker_point4, = absax.plot(primary_orbit.east_of_ecc(0), primary_orbit.north_of_ecc(0), 'ro')
    marker_point5, = absax.plot(secondary_orbit.east_of_ecc(0), secondary_orbit.north_of_ecc(0), 'ro')

    # noinspection PyTypeChecker
    axphase = plt.axes([0.05, 0.05, 0.9, 0.03])
    sphase = Slider(axphase, 'phase', np.float64(0.), np.float64(1.), valinit=0, valstep=0.005)

    def update(dummy):  # dummy argument is required by matplotlib
        phase = sphase.val
        newvrad1, newtheta1, newecc_anom1 = primary_orbit.radial_velocity_of_phase(phase, getAngles=True)
        newvrad2, newtheta2, newecc_anom2 = secondary_orbit.radial_velocity_of_phase(phase, getAngles=True)
        newecc_anomr = relative_orbit.eccentric_anom_of_phase(phase)
        marker_point1.set_xdata(phase)
        marker_point1.set_ydata(newvrad1)
        marker_point2.set_xdata(phase)
        marker_point2.set_ydata(newvrad2)
        marker_point3.set_xdata(relative_orbit.east_of_ecc(newecc_anomr))
        marker_point3.set_ydata(relative_orbit.north_of_ecc(newecc_anomr))
        marker_point4.set_xdata(primary_orbit.east_of_ecc(newecc_anom1))
        marker_point4.set_ydata(primary_orbit.north_of_ecc(newecc_anom1))
        marker_point5.set_xdata(secondary_orbit.east_of_ecc(newecc_anom2))
        marker_point5.set_ydata(secondary_orbit.north_of_ecc(newecc_anom2))
        fig.canvas.draw()

    sphase.on_changed(update)
