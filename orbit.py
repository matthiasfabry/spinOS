"""
Module that defines the Orbit class and subclasses.

The motion of a star in a (Newtonian) binary is completely specified by:
    1) a, the semi major axis
    2) e, the eccentricity
    3) i, the inclination with respect to the plane of the sky
    4) omega, the argument of periastron
    5) Omega, the longitude of the ascending node
    6) t0, the epoch of periastron passage
    7) P, the period of the orbit
    8) gamma, the systemic velocity of the binary pair
    9) d, the distance to the focus of the orbit (assumed way larger than the dimensions of the orbit)
"""
import numpy as np
from numpy import float64
import constants as c
import scipy.optimize as spopt


def generate_orbits(pars):
    """
    generates orbit objects given the parameters.
    :param pars: list containing, in order:
                        - e        the eccentricity
                        - i        the inclination (deg)
                        - omega    the longitude of the periastron with respect to the ascending node (deg)
                        - Omega    the longitude of the ascending node of the seconday measured east of north (deg)
                        - t0       the time of periastron passage (hjd)
                        - K1       the semiamplitude of the radial velocity curve of the primary (km/s)
                        - K2       the semiamplitude of the radial velocity curve of the secondary (km/s)
                        - P        the period of the binary (days)
                        - gamma1   the (apparent) systemic velocity of the primary (km/s)
                        - gamma2   the (apparent) systemic velocity of the secondary (km/s)
                        - d        the distance (pc)
    :return: the relative, primary and secondary orbit
    """
    relative = RelativeOrbit(pars[5] + pars[6], pars[0], pars[1], pars[2], pars[3], pars[4], pars[7], pars[10])
    primary = AbsoluteOrbit(pars[5], pars[0], pars[1], pars[2] + 180, pars[3], pars[4], pars[7], pars[8], pars[10])
    secondary = AbsoluteOrbit(pars[6], pars[0], pars[1], pars[2], pars[3], pars[4], pars[7], pars[9], pars[10])
    return relative, primary, secondary


class Orbit:
    """
    Creates a general orbit object, storing all the orbital elemements as well as its period and systemic velocity.
    Within binary orbital solution finding however, use the subclasses instead to differenciate the absolute and
    the relative orbits.
    """

    def __init__(self, k, e, i, omega, Omega, t0, P, gamma, d):
        self.k = k
        self.e = e
        self.i_deg = i
        self.i = self.i_deg * np.pi / 180
        self.sini = np.sin(self.i)
        self.cosi = np.cos(self.i)
        self.omega_deg = omega
        self.omega = self.omega_deg * np.pi / 180
        self.sino = np.sin(self.omega)
        self.coso = np.cos(self.omega)
        self.Omega_deg = Omega
        self.Omega = self.Omega_deg * np.pi / 180
        self.sinO = np.sin(self.Omega)
        self.cosO = np.cos(self.Omega)
        self.t0 = t0
        self.P = P
        self.gamma = gamma
        self.d = d
        self.a = (self.P * self.k * np.sqrt(1 - self.e ** 2)) / (2 * np.pi * self.sini) * 86400 / (
                c.pc * self.d) * 180 * 3600 / np.pi * 1000
        self.thiele_A = self.a * (self.cosO * self.coso - self.sinO * self.sino * self.cosi)
        self.thiele_B = self.a * (self.sinO * self.coso + self.cosO * self.sino * self.cosi)
        self.thiele_F = self.a * (-self.cosO * self.sino - self.sinO * self.coso * self.cosi)
        self.thiele_G = self.a * (-self.sinO * self.sino + self.cosO * self.coso * self.cosi)

    def radial_velocity_of_phase(self, phase, getAngles: bool = False):
        E = self.eccentric_anom_of_phase(phase)
        return self.radial_velocity_of_ecc_anom(E, getAngles)

    def radial_velocity_of_ecc_anom(self, ecc_anom, getAngles: bool = False):
        theta = self.true_anomaly_of_ecc_anom(ecc_anom)
        if getAngles:
            return self.k * (np.cos(theta + self.omega) + self.e * self.coso) + self.gamma, theta, ecc_anom
        return self.k * (np.cos(theta + self.omega) + self.e * self.coso) + self.gamma

    def radial_velocity_of_hjd(self, hjd, getAngles: bool = False):
        phase = (hjd - self.t0) % self.P / self.P
        return self.radial_velocity_of_phase(phase, getAngles=getAngles)

    def eccentric_anom_of_phase(self, phase):
        def keplers_eq(ph):
            def kepler(ecc_an):
                return (ecc_an - self.e * np.sin(ecc_an)) - 2 * np.pi * ph

            # you might need these derivatives if you change root finding algorithm
            def kepler_der(ecc_an):
                return 1 - self.e * np.cos(ecc_an)

            def kepler_der_2(ecc_an):
                return self.e * np.sin(ecc_an)

            return kepler, kepler_der, kepler_der_2

        if phase.size == 1:
            result = spopt.root_scalar(keplers_eq(phase)[0], method='toms748', bracket=(0, 2 * np.pi)).root
        else:
            result = np.zeros(len(phase))
            for i in range(len(phase)):
                # current root finding algorithm is toms748, as it has the fastest convergence.
                # You can change this as needed.
                result[i] = spopt.root_scalar(keplers_eq(phase[i])[0], method='toms748', bracket=(0, 2 * np.pi)).root
        return result

    def phase_of_hjd(self, hjd):
        return (hjd - self.t0) % self.P / self.P

    def phase_of_ecc_anom(self, ecc_anom):
        return (ecc_anom - self.e * np.sin(ecc_anom)) / (2 * np.pi)

    def true_anomaly_of_hjd(self, hjd):
        return self.true_anomaly_of_ecc_anom(self.eccentric_anom_of_phase(self.phase_of_hjd(hjd)))

    def true_anomaly_of_phase(self, phase):
        return self.true_anomaly_of_ecc_anom(self.eccentric_anom_of_phase(phase))

    def ecc_anom_of_true_anom(self, theta):
        return 2 * np.arctan(np.sqrt((1 - self.e) / (1 + self.e)) * np.tan(theta / 2))

    def true_anomaly_of_ecc_anom(self, E):
        return 2 * np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(E / 2))

    def X(self, E):
        return np.cos(E) - self.e

    def Y(self, E):
        return np.sqrt(1 - self.e ** 2) * np.sin(E)

    def r(self, theta):
        return self.a * (1 - self.e ** 2) / (1 + self.e * np.cos(theta))

    def x(self, theta):
        return self.r(theta) * (self.cosO * np.cos(theta + self.omega)
                                - self.sinO * np.sin(theta + self.omega) * self.cosi)

    def y(self, theta):
        return self.r(theta) * (self.sinO * np.cos(theta + self.omega)
                                + self.cosO * np.sin(theta + self.omega) * self.cosi)

    def z(self, theta):
        return -self.r(theta) * np.sin(theta + self.omega) * self.sini

    def north_of_ecc(self, E):
        return self.thiele_A * self.X(E) + self.thiele_F * self.Y(E)

    def east_of_ecc(self, E):
        return self.thiele_B * self.X(E) + self.thiele_G * self.Y(E)

    def north_of_true(self, theta):
        return self.north_of_ecc(self.ecc_anom_of_true_anom(theta))

    def east_of_true(self, theta):
        return self.east_of_ecc(self.ecc_anom_of_true_anom(theta))

    def north_of_hjd(self, hjd):
        return self.north_of_ecc(self.eccentric_anom_of_phase(self.phase_of_hjd(hjd)))

    def east_of_hjd(self, hjd):
        return self.east_of_ecc(self.eccentric_anom_of_phase(self.phase_of_hjd(hjd)))


class AbsoluteOrbit(Orbit):
    def __init__(self, k, e, i, omega, Omega, t0, P, gamma, d):
        super().__init__(k, e, i, omega, Omega, t0, P, gamma, d)


class RelativeOrbit(Orbit):
    def __init__(self, k, e, i, omega, Omega, t0, P, d):
        super().__init__(k, e, i, omega, Omega, t0, P, 0, d)
