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
import constants as c
import scipy.optimize as spopt


class Orbit:
    """
    Creates a general orbit object, storing all the orbital elemements as well as its period and systemic velocity.
    Within binary orbital solution finding however, use the subclasses instead to differenciate the absolute and
    the relative orbits.
    """

    def __init__(self, a, e, i, omega, Omega, t0, P, gamma, d):
        self.a = a
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
        self.K = 2 * self.a / (1000 * 3600) * np.pi / 180 * self.d * c.pc * np.pi * self.sini / (
                    self.P * 86400 * (1 - self.e ** 2))
        self.thiele_A = self.a * (self.cosO * self.coso - self.sinO * self.sino * self.cosi)
        self.thiele_B = self.a * (self.sinO * self.coso + self.cosO * self.sino * self.cosi)
        self.thiele_F = self.a * (-self.cosO * self.sino - self.sinO * self.coso * self.cosi)
        self.thiele_G = self.a * (-self.sinO * self.sino + self.cosO * self.coso * self.cosi)

    def radial_velocity(self, phase, getAngles: bool = False):
        E = self.eccentric_anom_of_phase(phase)
        theta = self.true_anomaly_of_ecc_anom(E)
        if getAngles:
            return self.K * (np.cos(theta + self.omega) + self.e * self.coso) + self.gamma, theta, E
        return self.K * (np.cos(theta + self.omega) + self.e * self.coso) + self.gamma

    def eccentric_anom_of_phase(self, phase):
        def keplers_eq(ph):
            def kepler(ecc_an):
                return (ecc_an - self.e * np.sin(ecc_an)) - 2 * np.pi * (ph - self.t0 / self.P)

            # you might need these derivatives if you change root finding algorithm
            def kepler_der(ecc_an):
                return 1 - self.e * np.cos(ecc_an)

            def kepler_der_2(ecc_an):
                return self.e * np.sin(ecc_an)

            return kepler, kepler_der, kepler_der_2

        # current root finding algorithm is toms748, as it has the fastest convergence. You can change this as needed.
        return spopt.root_scalar(keplers_eq(phase)[0], method='toms748', bracket=(0, 2 * np.pi)).root

    def true_anomaly_of_phase(self, phase):
        E = self.eccentric_anom_of_phase(phase)
        return self.true_anomaly_of_ecc_anom(E)

    def ecc_anom_of_true_anom(self, theta):
        return 2 * np.arctan(np.sqrt((1 - self.e) / (1 + self.e)) * np.tan(theta / 2))

    def true_anomaly_of_ecc_anom(self, E):
        return 2 * np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(E / 2))

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

    def north_ecc(self, E):
        return self.thiele_A * self.X(E) + self.thiele_F * self.Y(E)

    def east_ecc(self, E):
        return self.thiele_B * self.X(E) + self.thiele_G * self.Y(E)

    def north_true(self, theta):
        return self.north_ecc(self.ecc_anom_of_true_anom(theta))

    def east_true(self, theta):
        return self.east_ecc(self.ecc_anom_of_true_anom(theta))

    def X(self, E):
        return np.cos(E) - self.e

    def Y(self, E):
        return np.sqrt(1 - self.e ** 2) * np.sin(E)


class PrimaryOrbit(Orbit):
    def __init__(self, k1, e, i, omega1, Omega1, t0, P, gamma1, d):
        super().__init__(None, e, i, omega1, Omega1, t0, P, gamma1, d)


class SecondaryOrbit(Orbit):
    def __init__(self, k2, e, i, omega2, Omega2, t0, P, gamma2, d):
        super().__init__(None, e, i, omega2, Omega2, t0, P, gamma2, d)


class RelativeOrbit(Orbit):
    def __init__(self, a, e, i, omega, Omega, t0, P, d):
        super().__init__(a, e, i, omega, Omega, t0, P, 0, d)
