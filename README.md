presenting spinOS: the SPectroscopic and INterferometric Orbital Solution finder.

Goal:
    spinOS computes the best fit orbital solution given:
     1) radial and systemic velocity data for either or both components of a spectroscopic binary and/or
     2) astrometric data containing the separtions and error ellipses, and
     3) an initial guess of the parameters:
        - e        the eccentricity
        - i        the inclination (deg)
        - omega    the longitude of the periastron with respect to the ascending node of the secondary (deg)
        - Omega    the longitude of the ascending node of the seconday measured east of north (deg)
        - t0       the time of periastron passage (day)
        - k1       the semiamplitude of the radial velocity curve of the primary (km/s)
        - k2       the semiamplitude of the radial velocity curve of the secondary (km/s)
        - p        the period of the binary (day)
        - gamma1   the peculiar velocity of the primary (km/s)
        - gamma2   the peculiar velocity of the secondary (km/s)
        - d        the distance (pc)
        For each parameter, a tag True/False should be supplied to decide for the minimizer to vary this parameter. (So
        False means it will keep it fixed at the supplied value)

This package provides a commandline script interface as well as a GUI to easily visualize your data

Dependencies:
    python 3.7
    numpy 1.17.2
    scipy 1.3.1
    lmfit 0.9.14
    matplotlib 3.1.1
    emcee 3.0.0

Author:
    Matthias Fabry
    Instituut voor Sterrekunde, KU Leuven, Belgium

Date:
    21 Nov 2019

Version:
    1.5

Acknowledgements:
    This python3 implementation is heavily based on an earlier IDL implementation by Hugues Sana.
    We thank the authors of lmfit for the development of their package.
