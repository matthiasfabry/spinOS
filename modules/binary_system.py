"""
Copyright 2020-2024 Matthias Fabry
This file is part of spinOS.

spinOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

spinOS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with spinOS. If not, see <https://www.gnu.org/licenses/>.

Module that defines the System class, the Orbit class and its subclasses.
"""
import numpy as np
import scipy.optimize as spopt

import modules.constants as const


class BinarySystem:
    """
    The class that represents a binary system with its respective
    components. It is assumed that this binary has a
    distance to the observer that is way larger than the orbital separation.
    """

    def __init__(self, parameters: dict):
        """
        Creates a BinarySystem object, defining a binary system with the 11
        parameters supplied that fully determine the orbits:
                - e        the eccentricity
                - i        the inclination (deg)
                - omega    the longitude of the periastron with respect to
                the ascending node (deg)
                - Omega    the longitude of the ascending node of the
                secondary measured east of north (deg)
                - t0       the time of periastron passage (hjd)
                - p        the period of the binary (days)
                - d        the distance to the system(pc)
                - k1       the semi-amplitude of the radial velocity curve of
                the primary (km/s)
                - k2       the semi-amplitude of the radial velocity curve of
                the secondary (km/s)
                - gamma1   the (apparent) systemic velocity of the primary (
                km/s)
                - gamma2   the (apparent) systemic velocity of the secondary
                (km/s)
                - mt       the total dynamical mass of the system (Msun)
        :param parameters: dictionary containing the aforementioned parameters
        """
        if parameters['p'] == 0.0:
            raise ValueError('a binary system cannot have a period of zero days')
        self.e = parameters['e']
        self.i = parameters['i'] * const.DEG2RAD
        self.sini = np.sin(self.i)
        self.cosi = np.cos(self.i)
        self.Omega = parameters['Omega'] * const.DEG2RAD
        self.sinO = np.sin(self.Omega)
        self.cosO = np.cos(self.Omega)
        self.t0 = parameters['t0']
        self.p = parameters['p']  # day
        self.d = parameters['d'] * const.PC2KM  # km
        self.primary: AbsoluteOrbit = AbsoluteOrbit(self, abs(parameters['k1']),
                                                    parameters['omega'] - 180, parameters['gamma1'])
        self.secondary: AbsoluteOrbit = AbsoluteOrbit(self, abs(parameters['k2']),
                                                      parameters['omega'], parameters['gamma2'])
        self.ap_kk = (self.p * const.DAY2SEC) * (
                abs(parameters['k1']) + abs(parameters['k2'])) * np.sqrt(1 - self.e ** 2) / (
                             2 * np.pi * abs(self.sini))
        self.mt = parameters['mt']
        self.ap_mt = np.cbrt(
            parameters['mt'] * const.MSUN * const.G * (self.p * const.DAY2SEC) ** 2 / (
                    4 * np.pi ** 2))
        self.relative: RelativeOrbit = RelativeOrbit(self, self.ap_mt / self.d * const.RAD2MAS,
                                                     parameters['omega'])

    def extend_rvs_until_time(self, times, rvs, maxtime):
        out = np.copy(times)
        n = int((maxtime - times[0]) // self.p)
        for i in range(n):
            times += self.p
            out = np.concatenate((out, times))
        rvs = np.tile(rvs, n + 1)
        return out, rvs

    def semimajor_axis_from_RV(self):
        """
        Calculates the physical semi-major axis of the relative orbit.
        :return: semi-major axis (AU)
        """
        return self.ap_kk / const.AU2KM

    def semimajor_axis_from_distance(self):
        """
        Calculates the physical semi-major axis of the relative orbit from the distance and
        apparent size.
        :return: semi-major axis (AU)
        """
        return self.ap_mt / const.AU2KM

    def primary_mass(self):
        """
        Calculates the mass of the primary body of the system
        :return: mass of the primary (Solar Mass)
        """
        return np.power(1 - self.e ** 2, 1.5) * (
                self.primary.k + self.secondary.k) ** 2 * self.secondary.k * (
                self.p * const.DAY2SEC) / (2 * np.pi * const.G * abs(self.sini) ** 3) / const.MSUN

    def secondary_mass(self):
        """
        Calculates the mass of the secondary body of the system
        :return: mass of the secondary (in Solar Mass)
        """
        return np.power(1 - self.e ** 2, 1.5) * (
                self.primary.k + self.secondary.k) ** 2 * self.primary.k * (
                self.p * const.DAY2SEC) / (2 * np.pi * const.G * abs(self.sini) ** 3) / const.MSUN

    def total_mass(self):
        """
        Calculates the total mass of the system
        :return: mass of the system (in Solar Mass)
        """
        return np.power(1 - self.e ** 2, 1.5) * (self.primary.k + self.secondary.k) ** 3 * (
                self.p * const.DAY2SEC) / (2 * np.pi * const.G * abs(self.sini) ** 3) / const.MSUN

    def total_mass_from_distance(self):
        """
        Calculates the total dynamical mass of the system using the size of
        the apparent orbit.
        :return: m_total (Msun)
        """
        return 4 * np.pi ** 2 * np.power(self.d * self.relative.a * const.MAS2RAD, 3) / (
                const.G * (self.p * const.DAY2SEC) ** 2) / const.MSUN

    def phase_of_hjd(self, hjds):
        """
        Calculates the phase in the orbit given a Julian Date.
        :param hjds: Julian Date (day)
        :return: phase (rad)
        """
        return np.remainder((hjds - self.t0) / self.p, 1)

    def phase_of_ecc_anom(self, ecc_anoms):
        """
        Calculates the phase in the orbit given an eccentric anomaly.
        :param ecc_anoms: eccentric anomaly (rad)
        :return: phase (rad)
        """
        return (ecc_anoms - self.e * np.sin(ecc_anoms)) / (2 * np.pi)

    def ecc_anom_of_true_anom(self, theta):
        """
        Calculates the eccentric anomaly given a true anomaly
        :param theta: true anomaly (rad)
        :return: eccentric anomaly (rad)
        """
        return 2 * np.arctan(np.sqrt((1 - self.e) / (1 + self.e)) * np.tan(theta / 2))

    def true_anomaly_of_ecc_anom(self, E):
        """
        Calculates the true anomaly given an eccentric anomaly
        :param E: eccentric anomaly (rad)
        :return: true anomaly (rad)
        """
        return 2 * np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(E / 2))

    def ecc_anom_of_phase(self, phase):
        """
        Calculates the eccentric anomaly given a phase. This function is the hardest function to
        execute as it requires solving a transcendental equation: Kepler's Equation.
        :param phase: phase (rad)
        :return: eccentric anomaly (rad)
        """

        # define keplers equation as function of a phase
        def keplers_eq(ph):
            """
            wrapper that returns kepler's equation of a phase ph as a function object
            :param ph: phase
            :return: python function to find root of
            """

            # build a function object that should be zero for a certain eccentric anomaly
            def kepler(ecc_an):
                """
                defines kepler's equation in function of an eccentric anomaly and phase;
                if zero, the given eccentric anomaly corresponds to the phase ph of this orbit
                :param ecc_an: eccentric anomaly
                :return: float
                """
                return ecc_an - self.e * np.sin(ecc_an) - 2 * np.pi * ph

            return kepler

        # find the root of keplers_eq(phase), which by construction returns
        # a function for which the eccentric anomaly is the independent variable.
        # current root finding algorithm is toms748, as it has the best
        # convergence (2.7 bits per function evaluation)
        phase = np.remainder(phase, 1)
        try:
            _ = iter(phase)
        except TypeError:
            # 'not iterable'
            return spopt.root_scalar(keplers_eq(phase), method='toms748',
                                     bracket=(0, 2 * np.pi)).root
        else:
            # iterable
            Es = np.empty(len(phase))
            for i in range(len(phase)):
                Es[i] = spopt.root_scalar(keplers_eq(phase[i]), method='toms748',
                                          bracket=(0, 2 * np.pi)).root
            return Es

    def create_phase_extended_RV(self, rvdata, extension_range):
        """
        Creates a new RV dataset, where the phase folding of the data is
        extended outside of the (0, 1) interval by a given amount
        :param rvdata: an ndarray containing the hjds, RV measurements and
        errors
        :param extension_range: the phase amount to extend the folding with.
        :return: same dataset as supplied, only folded to phases (
        -extension_range,
        1+extension_range)
        """
        phases = self.phase_of_hjd(rvdata[:, 0])
        data = rvdata[:, 1]
        errors = rvdata[:, 2]
        left_extended_phases = phases[phases > (1 - extension_range)] - 1
        right_extended_phases = phases[phases < extension_range] + 1
        left_extended_data = data[phases > (1 - extension_range)]
        right_extended_data = data[phases < extension_range]
        left_extended_errors = errors[phases > (1 - extension_range)]
        right_extended_errors = errors[phases < extension_range]
        extended_phases = np.concatenate((left_extended_phases, phases, right_extended_phases))
        extended_data = np.concatenate((left_extended_data, data, right_extended_data))
        extended_errors = np.concatenate((left_extended_errors, errors, right_extended_errors))
        return extended_phases, extended_data, extended_errors


class Orbit:
    """
    Creates a general orbit object, storing all the orbital elements as
    well as its period and systemic velocity.
    Within binary orbital solution finding however, use the subclasses
    instead to differentiate the absolute and the relative orbits.

    The motion of a star in a (Newtonian) binary is completely determined by
    the quantities:
        -) e, the eccentricity
        -) i (rad),  the inclination with respect to the plane of the sky
        -) omega (rad), the argument of periastron
        -) Omega (rad), the longitude of the ascending node
        -) t0 (day), the epoch of periastron passage
        -) k (km/s), the semi-amplitude of the observed radial velocity curve
        -) p (day), the period of the orbit
        -) gamma (km/s), the systemic velocity of the binary pair
        -) d (pc), the distance to the focus of the orbit (assumed way
        larger than the dimensions
        of the orbit)

    An orbit inherits directly the quantities e, i, Omega, t0, p, d from the
    system it resides in.
    """

    def __init__(self, system, omega):
        self.system: BinarySystem = system
        self.omega = omega * const.DEG2RAD
        self.sino = np.sin(self.omega)
        self.coso = np.cos(self.omega)


class AbsoluteOrbit(Orbit):
    """
    An absolute orbit represents an orbit of either of the component
    masses separately. It is these orbits that
    determine the observed RV measurements
    """

    def __init__(self, system, k, omega, gamma):
        """
        An absolute orbit of a component body is an orbit with the systemic
        parameters as well as k, omega, gamma.
        :param system: system the component belongs to
        :param k: semiamplitude of its RV curve
        :param omega: its argument of periastron
        :param gamma: its peculiar velocity
        """
        self.k = k
        self.gamma = gamma
        super().__init__(system, omega)

    def radial_velocity_of_phase(self, phase, getAngles: bool = False):
        """
        Calculates the radial velocities of a component body given a list of
        phases
        :param phase: phase (rads) (cannot be an iterable)
        :param getAngles: Boolean (default=False) indicating to additionally
        return the
        corresponding true and eccentric
            anomalies
        :return: list of radial velocities (km/s) [optionally: list of true
        anomalies (rad) and
        list of eccentric
            anomalies (rad)]
        """
        E = self.system.ecc_anom_of_phase(phase)
        return self.radial_velocity_of_ecc_anom(E, getAngles)

    def radial_velocity_of_ecc_anom(self, ecc_anom, getAngles: bool = False):
        """
        Calculates the radial velocity of a component body given an
        eccentric anomaly
        :param ecc_anom: eccentric anomaly (rad)
        :param getAngles: Boolean (default=False) indicating to additionally
        return the eccentric
        anomaly
        :return: radial velocity (km/s) [optionally: eccentric anomaly (rad)]
        """
        if getAngles:
            return self.radial_velocity_of_true_anom(self.system.true_anomaly_of_ecc_anom(ecc_anom),
                                                     getAngles), ecc_anom
        return self.radial_velocity_of_true_anom(self.system.true_anomaly_of_ecc_anom(ecc_anom))

    def radial_velocity_of_true_anom(self, theta, getAngles: bool = False):
        """
        Calculates the radial velocity of a component body given a true anomaly
        :param theta: true anomaly (rad)
        :param getAngles: Boolean (default=False) indicating to additionally
        return the true anomaly
        :return: radial velocity (km/s) [optionally: true anomaly (rad)]
        """
        if getAngles:
            return self.k * (
                    np.cos(theta + self.omega) + self.system.e * self.coso) + self.gamma, theta
        return self.k * (np.cos(theta + self.omega) + self.system.e * self.coso) + self.gamma

    def radial_velocity_of_hjd(self, hjd, getAngles: bool = False):
        """
        Calculates the radial velocity of a component body given Julian dates
        :param hjd: julian dates
        :param getAngles: Boolean (default=False) indicating to additionally
        return the true anomaly
        :return: radial velocity (km/s) [optionally: true anomaly (rad)]
        """
        return self.radial_velocity_of_phase(self.system.phase_of_hjd(hjd), getAngles=getAngles)


class RelativeOrbit(Orbit):
    """
    A Relative orbit represents the relative orbit of the secondary with
    respect to the primary. It is this orbit that determine the observed astrometric measurements
    """

    def __init__(self, system, a, omega):
        super().__init__(system, omega)
        self.a = a
        self.thiele_A = self.a * (system.cosO * self.coso - system.sinO * self.sino * system.cosi)
        self.thiele_B = self.a * (system.sinO * self.coso + system.cosO * self.sino * system.cosi)
        self.thiele_F = self.a * (-system.cosO * self.sino - system.sinO * self.coso * system.cosi)
        self.thiele_G = self.a * (-system.sinO * self.sino + system.cosO * self.coso * system.cosi)

    def north_of_phase(self, ph):
        """
        Calculates the northward separations given a phase
        :param ph: phases
        :return: northward separations
        """
        return self.north_of_ecc(self.system.ecc_anom_of_phase(ph))

    def east_of_phase(self, ph):
        """
        Calculates the eastward separations given a phase
        :param ph: phases
        :return: eastward separations
        """
        return self.east_of_ecc(self.system.ecc_anom_of_phase(ph))

    def north_of_ecc(self, E):
        """
        Calculates the northward separations given an eccentric anomaly
        :param E: eccentric anomalies (rad)
        :return: northward separations
        """
        return self.thiele_A * self.X(E) + self.thiele_F * self.Y(E)

    def east_of_ecc(self, E):
        """
        Calculates the eastward separations given an eccentric anomaly
        :param E: eccentric anomaly (rad)
        :return: eastward separation
        """
        return self.thiele_B * self.X(E) + self.thiele_G * self.Y(E)

    def north_of_true(self, theta):
        """
        Calculates the northward separation given a true anomaly
        :param theta: true anomaly (rad)
        :return: northward separation
        """
        return self.north_of_ecc(self.system.ecc_anom_of_true_anom(theta))

    def east_of_true(self, theta):
        """
        Calculates the eastward separation given a true anomaly
        :param theta: eccentric anomaly (rad)
        :return: eastward separation
        """
        return self.east_of_ecc(self.system.ecc_anom_of_true_anom(theta))

    def north_of_hjd(self, hjd):
        """
        Calculates the northward separation given Julian dates
        :param hjd: julian date (days)
        :return: northward separations
        """
        return self.north_of_ecc(self.system.ecc_anom_of_phase(self.system.phase_of_hjd(hjd)))

    def east_of_hjd(self, hjd):
        """
        Calculates the eastward separation given Julian dates
        :param hjd: julian date (days)
        :return: eastward separations
        """
        return self.east_of_ecc(self.system.ecc_anom_of_phase(self.system.phase_of_hjd(hjd)))

    def separation_of_hjd(self, hjd):
        """
        Calculates the total separation (in mas) of the two stars at given Julian dates
        :param hjd: julian dates (days)
        :return: total separations (deg)
        """
        east = self.east_of_hjd(hjd)
        north = self.north_of_hjd(hjd)
        return np.sqrt(east ** 2 + north ** 2)

    def position_angle_of_hjd(self, hjd):
        """
        Calculates teh position angles (east of north) of the secondary wrt to the primary at given
        Julian dates
        :param hjd: julian dates (days)
        :return: position angles (deg)
        """
        east = self.east_of_hjd(hjd)
        north = self.north_of_hjd(hjd)
        return np.mod(np.arctan2(east, north) * const.RAD2DEG, 360)

    def X(self, E):
        """
        Calculates the elliptical rectangular coordinate X given an eccentric anomaly
        :param E: eccentric anomaly (rad)
        :return: elliptical rectangular coordinate X
        """
        return np.cos(E) - self.system.e

    def Y(self, E):
        """
        Calculates the elliptical rectangular coordinate Y given an eccentric anomaly
        :param E: eccentric anomaly (rad)
        :return: elliptical rectangular coordinate Y
        """
        return np.sqrt(1 - self.system.e ** 2) * np.sin(E)
