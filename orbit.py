"""
Module that defines the System class, the Orbit class and its subclasses.

The motion of a star in a (Newtonian) binary is completely determined by the quantities:
    1) k, the semiamplitude of the observed radial velocity curve
    2) e, the eccentricity
    3) i, the inclination with respect to the plane of the sky
    4) omega, the argument of periastron
    5) Omega, the longitude of the ascending node
    6) t0, the epoch of periastron passage
    7) P, the period of the orbit
    8) gamma, the systemic velocity of the binary pair
    9) d, the distance to the focus of the orbit (assumed way larger than the dimensions of the orbit)
"""
import constants as c
import scipy.optimize as spopt
import numpy as np


class System:
    """
    generates orbit objects given the parameters.
    :param parameters: dictionary containing:
                        - e        the eccentricity
                        - i        the inclination (deg)
                        - omega    the longitude of the periastron with respect to the ascending node (deg)
                        - Omega    the longitude of the ascending node of the seconday measured east of north (deg)
                        - t0       the time of periastron passage (hjd)
                        - k1       the semiamplitude of the radial velocity curve of the primary (km/s)
                        - k2       the semiamplitude of the radial velocity curve of the secondary (km/s)
                        - p        the period of the binary (days)
                        - gamma1   the (apparent) systemic velocity of the primary (km/s)
                        - gamma2   the (apparent) systemic velocity of the secondary (km/s)
                        - d        the distance (pc)
    :return: the relative, primary and secondary orbit
    """

    def __init__(self, parameters: dict):
        self.e = parameters['e']
        self.i_deg = parameters['i']
        self.i = self.i_deg * np.pi / 180
        self.sini = np.sin(self.i)
        self.cosi = np.cos(self.i)
        self.Omega_deg = parameters['Omega']
        self.Omega = self.Omega_deg * np.pi / 180
        self.sinO = np.sin(self.Omega)
        self.cosO = np.cos(self.Omega)
        self.t0 = parameters['t0']
        self.p = parameters['p']
        self.d = parameters['d']
        self.primary = AbsoluteOrbit(self, parameters['k1'], parameters['omega'] + 180, parameters['gamma1'])
        self.secondary = AbsoluteOrbit(self, parameters['k2'], parameters['omega'], parameters['gamma2'])
        self.relative = RelativeOrbit(self, parameters['k1'] + parameters['k2'], parameters['omega'])


class Orbit:
    """
    Creates a general orbit object, storing all the orbital elemements as well as its period and systemic velocity.
    Within binary orbital solution finding however, use the subclasses instead to differenciate the absolute and
    the relative orbits.
    """

    def __init__(self, system, k, omega, gamma, ):
        self.system = system
        self.k = k
        self.omega_deg = omega
        self.omega = self.omega_deg * np.pi / 180
        self.sino = np.sin(self.omega)
        self.coso = np.cos(self.omega)
        self.gamma = gamma

    def radial_velocity_of_phase(self, phase, getAngles: bool = False):
        E = self.eccentric_anom_of_phase(phase)
        return self.radial_velocity_of_ecc_anom(E, getAngles)

    def radial_velocity_of_ecc_anom(self, ecc_anom, getAngles: bool = False):
        theta = self.true_anomaly_of_ecc_anom(ecc_anom)
        if getAngles:
            return self.k * (np.cos(theta + self.omega) + self.system.e * self.coso) + self.gamma, theta, ecc_anom
        return self.k * (np.cos(theta + self.omega) + self.system.e * self.coso) + self.gamma

    def radial_velocity_of_hjd(self, hjd, getAngles: bool = False):
        return self.radial_velocity_of_phase(self.phase_of_hjd(hjd), getAngles=getAngles)

    def eccentric_anom_of_phase(self, phase):
        def keplers_eq(ph):
            def kepler(ecc_an):
                return (ecc_an - self.system.e * np.sin(ecc_an)) - 2 * np.pi * ph

            return kepler

        if phase.size == 1:
            result = spopt.root_scalar(keplers_eq(phase), method='toms748', bracket=(0, 2 * np.pi)).root
        else:
            result = np.zeros(len(phase))
            for i in range(len(phase)):
                # current root finding algorithm is toms748, as it has the fastest convergence.
                # You can change this as needed.
                result[i] = spopt.root_scalar(keplers_eq(phase[i]), method='toms748', bracket=(0, 2 * np.pi)).root
        return result

    def phase_of_hjd(self, hjd):
        return (hjd - self.system.t0) % self.system.p / self.system.p

    def phase_of_ecc_anom(self, ecc_anom):
        return (ecc_anom - self.system.e * np.sin(ecc_anom)) / (2 * np.pi)

    def true_anomaly_of_hjd(self, hjd):
        return self.true_anomaly_of_ecc_anom(self.eccentric_anom_of_phase(self.phase_of_hjd(hjd)))

    def true_anomaly_of_phase(self, phase):
        return self.true_anomaly_of_ecc_anom(self.eccentric_anom_of_phase(phase))

    def ecc_anom_of_true_anom(self, theta):
        return 2 * np.arctan(np.sqrt((1 - self.system.e) / (1 + self.system.e)) * np.tan(theta / 2))

    def true_anomaly_of_ecc_anom(self, E):
        return 2 * np.arctan(np.sqrt((1 + self.system.e) / (1 - self.system.e)) * np.tan(E / 2))


class AbsoluteOrbit(Orbit):
    def __init__(self, system, k, omega, gamma):
        super().__init__(system, k, omega, gamma)


class RelativeOrbit(Orbit):
    def __init__(self, system, k, omega):
        super().__init__(system, k, omega, 0)
        self.a = (system.p * self.k * np.sqrt(1 - system.e ** 2)) / (2 * np.pi * system.sini) * 24 / (
                c.pc2km * system.d) * 180000 / np.pi  # mas separation of the relative
        self.thiele_A = self.a * (system.cosO * self.coso - system.sinO * self.sino * system.cosi)
        self.thiele_B = self.a * (system.sinO * self.coso + system.cosO * self.sino * system.cosi)
        self.thiele_F = self.a * (-system.cosO * self.sino - system.sinO * self.coso * system.cosi)
        self.thiele_G = self.a * (-system.sinO * self.sino + system.cosO * self.coso * system.cosi)

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

    def X(self, E):
        return np.cos(E) - self.system.e

    def Y(self, E):
        return np.sqrt(1 - self.system.e ** 2) * np.sin(E)

    def r(self, theta):
        return self.a * (1 - self.system.e ** 2) / (1 + self.system.e * np.cos(theta))

    def x(self, theta):
        return self.r(theta) * (self.system.cosO * np.cos(theta + self.omega)
                                - self.system.sinO * np.sin(theta + self.omega) * self.system.cosi)

    def y(self, theta):
        return self.r(theta) * (self.system.sinO * np.cos(theta + self.omega)
                                + self.system.cosO * np.sin(theta + self.omega) * self.system.cosi)

    def z(self, theta):
        return -self.r(theta) * np.sin(theta + self.omega) * self.system.sini
