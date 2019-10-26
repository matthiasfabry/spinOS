"""
This module plots a Radial velocity curve, as well as the corresponding orbit.
The module is used as follows:
$ python plot_RV_curve.py a e i omega Omega T gamma P
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
import orbit
from matplotlib.widgets import Slider

# read in arguments in standalone mode.
j = 1
if len(sys.argv) > 11:
    print('too many arguments!')
    exit()
try:
    a = np.float64(sys.argv[j])  # (mas)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    e = np.float64(sys.argv[j])
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    i = np.float64(sys.argv[j])  # (deg)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    omega = np.float64(sys.argv[j])  # (deg)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    Omega = np.float64(sys.argv[j])  # (deg)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    t0 = np.float64(sys.argv[j])  # (days)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    gamma = np.float64(sys.argv[j])  # (km/s)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    P = np.float64(sys.argv[j])  # (days)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    d = np.float64(sys.argv[j])  # (pc)
    j += 1
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()
try:
    do3d = bool(sys.argv[j] == 'True')
except IndexError:
    print("not enough arguments, was scanning argument", j)
    exit()

print("plotting orbit for:")
print("-------------------")
# noinspection PyUnboundLocalVariable
print("semimajor(mas) ", a, "\neccentricity ", e, "\ninclination(degrees) ", i, "\nomega(degrees) ", omega,
      "\nOmega(degrees) ", Omega, "\ntime of periastron(HJD) ", t0,
      "\nsystemic velocity(km/s)", gamma, "\nPeriod(days) ", P, "\nDistance(pc) ", d, "\nmaking 3d plot ", do3d)

orbit = orbit.RelativeOrbit(a, e, i, omega, Omega, t0, P, d)

# calculate RV curve

phases = np.linspace(0, 1, 500)
vrads = np.zeros(len(phases))
thetas = np.zeros(len(phases))
ecc_anoms = np.zeros(len(phases))
for i in range(len(vrads)):
    vrads[i], thetas[i], ecc_anoms[i] = orbit.radial_velocity(phases[i], getAngles=True)

norths = orbit.north_ecc(ecc_anoms)
easts = orbit.east_ecc(ecc_anoms)

# plot RV curve
plt.rc('text', usetex=True)
fig = plt.figure()
plt.subplots_adjust(bottom=0.25)
ax1 = fig.add_subplot(121)
ax1.plot(phases, vrads)
ax1.plot(phases, np.ones(len(phases)) * gamma, 'g')
ax1.set_xlabel(r'$phase$')
ax1.set_ylabel(r'$RV (km s^{-1})$')
ax1.grid()

# make 2D astrometric plot

ax3 = fig.add_subplot(122, aspect=1)
ax3.plot(easts, norths)
ax3.plot([orbit.east_true(-orbit.omega), orbit.east_true(-orbit.omega + np.pi)],
         [orbit.north_true(-orbit.omega), orbit.north_true(-orbit.omega + np.pi)],
         'k--', label='line of nodes')
ax3.plot([orbit.east_ecc(0), orbit.east_ecc(np. pi)],
         [orbit.north_ecc(0), orbit.north_ecc(np.pi)],
         'r', label='major axis')
ax3.plot([orbit.east_ecc(0)], [orbit.north_ecc(0)], marker='s', fillstyle='none')
ax3.set_xlabel('East')
ax3.set_xlim(ax3.get_xlim()[::-1])
ax3.set_ylabel('North')
ax3.axhline(linestyle=':', color='black')
ax3.axvline(linestyle=':', color='black')

# make a Slider
marker_point1, = ax1.plot(phases[0], vrads[0], 'ro')
marker_point3, = ax3.plot([easts[0]], [norths[0]], 'ro')

# noinspection PyTypeChecker
axphase = plt.axes([0.1, 0.1, 0.75, 0.03])
sphase = Slider(axphase, 'phase', 0, 1, valinit=0, valstep=0.01)


def update(dummy):  # dummy argument is required by matplotlib
    phase = sphase.val
    newvrad, newtheta, newecc_anom = orbit.radial_velocity(phase, getAngles=True)
    marker_point1.set_ydata(newvrad)
    marker_point1.set_xdata(phase)
    marker_point3.set_xdata(orbit.east_ecc(newecc_anom))
    marker_point3.set_ydata(orbit.north_ecc(newecc_anom))
    fig.canvas.draw_idle()


sphase.on_changed(update)

if do3d:
    plot_3d_orbit.plot_3d(orbit)

plt.show()
